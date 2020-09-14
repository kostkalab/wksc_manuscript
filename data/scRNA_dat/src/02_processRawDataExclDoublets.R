############################################################
## finally
############################################################
library(DropletUtils)
library(scater)
library(scran)
library(uwot)
library(pheatmap)
library(scds)

############################################################
# get counts from the molecule information file
############################################################
fname <- "../cr_out/outs/molecule_info.h5"
mol.info <- read10xMolInfo(fname)

## UMI counts before excluding empty droplets
cnts <- makeCountMatrix(gene = mol.info$data$gene,
                        cell = mol.info$data$cell,
                        all.genes = mol.info$genes)
dim(cnts) ## 27998 x 520454

############################################################
## (log)-total UMI vs. (log)-rank of each barcode
############################################################
br.out <- barcodeRanks(cnts)

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", main = "(log)-total UMI vs. (log)-rank")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", legend=c("knee", "inflection"), col=c("dodgerblue", "forestgreen"), lty=2)

############################################################
## detect & remove empty droplets with FDR 0.01
############################################################
set.seed(1000) ## for permutation testing to calculate p-values
get.ed = emptyDrops(cnts, lower = 100)
is.cell <- get.ed$FDR <= 0.01
sum(is.cell, na.rm=TRUE) ## 5887
table(Limited = get.ed$Limited, Significant = is.cell) ## no Limited==TRUE and Sig==FALSE

plot(get.ed$Total,
     -get.ed$LogProb,
     col = ifelse(is.cell, rgb(1,0,0,1/5), rgb(0,0,0,1/5)),
     xlab = "Total UMI count",
     ylab = "-Log Probability",pch=".")
## hmm..

## Exclude empty droplets
cnts = cnts[, which(is.cell)]; dim(cnts) ## 27998 x 5918

tmp = cnts
colnames(tmp) = NULL
rownames(tmp) = NULL

############################################################
## create SCE object
############################################################
sce = SingleCellExperiment(assays = list(counts = tmp))
rowData(sce)$ID = as.character(rownames(cnts))
colData(sce)$cell = as.character(colnames(cnts))

annot.df <- readRDS("../dat/annot.rds")
annot.df <- annot.df[rowData(sce)$ID,  2:ncol(annot.df)]

rowData(sce)$chr <- annot.df$chrom
rowData(sce)$symbol <- annot.df$gene_name
rowData(sce)$symbol = tolower(rowData(sce)$symbol)
rowData(sce)$feature <- annot.df$feature
rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$symbol)

## restrict chrs
sce = sce[rowData(sce)$chr %in% c(as.character(1:22),"X","Y","MT"),]; dim(sce)
## 27921 x 5918

############################################################
## QC metrics
############################################################
sce = calculateQCMetrics(sce, feature_controls = list(Mito = which(rowData(sce)$chr == "MT")))

par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
    xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

############################################################
## cell filtering -- mark cells by QC metrics
############################################################
## no total library filter + total no. of expressed features + mito content
low.nexprs <- sce$total_features_by_counts < 1000 ## colSums(counts(sce) > 0) < 1000
table(low.nexprs) ## 1428

lmed = median(log(sce$pct_counts_Mito))
lmad = mad(log(sce$pct_counts_Mito))
lcut = lmed + 3*lmad
cut  = exp(lcut) ; rm(lmed,lmad,lcut)
high.mito <- sce$pct_counts_Mito > cut
table(high.mito) ## 214

discard <- low.nexprs | high.mito ## 1515
colData(sce)$PassQC_dk <- !discard

## filter cells
sce <- sce[, sce$PassQC_dk]; dim(sce) ## 27921 x 4403

############################################################
## cell filtering -- mark putative doublets and filter
############################################################
#- Annotate doublet using co-expression based doublet scoring:
sce = cxds(sce)

#- Annotate doublet using binary classification based doublet scoring:
set.seed(1234) ## remember we need to set seed here!
sce = bcds(sce)

#- Combine both annotations into a hybrid annotation
sce = cxds_bcds_hybrid(sce)

#- Doublet scores are now available via colData:
CD  = colData(sce)
head(cbind(CD$cxds_score,CD$bcds_score, CD$hybrid_score))

frc = 0.05
numpos = round(nrow(CD) * frc)
possible_doublet = rep(FALSE, nrow(CD))
possible_doublet[order(CD$hybrid_score, decreasing = TRUE)[1:numpos]] = TRUE
sce$possible_doublet = possible_doublet
table(sce$possible_doublet) ## 220

## exclude doublets
inds <- !sce$possible_doublet
sce <- sce[, inds]; dim(sce) ## 27921 x 4183

############################################################
## gene filtering
############################################################
## genes with at least 3 counts in at least 3 cells
filtGenes <- function(obj){
    return(obj[(Matrix::rowSums(counts(obj) >= 3) >= 3), ])
}

sce <- filtGenes(sce) ; dim(sce)  ## 11155 x 4182

############################################################
## normalize and plot
############################################################
normalizeSCE <- function(obj){
    set.seed(1000)
    clusters = quickCluster(obj, method="igraph",
                            min.mean=0.1)
    obj = computeSumFactors(obj, min.mean = 0.1, cluster = clusters)
    obj = normalize(obj)

    plot(sizeFactors(obj),
         obj$total_counts/1e3,
         log = "xy",
         ylab = "Library size (thousands)",
         xlab = "Size factor")

    return(obj)
}

sce <- normalizeSCE(sce)

############################################################
## modeling mean-variance trend
############################################################
## endo-based with loess span set low
fit <- trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))

## decompose variance and get HVG by biol with dk's cutoffs
res <- decomposeVar(sce, fit = fit)
ii <- which((res$FDR < 0.01) & (pmax(0, res$bio)/res$total >= 0.25))
hvg  = sort(rownames(res)[ii]); length(hvg) ## 827; ## 818

## plot
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
points(fit$mean[hvg], fit$var[hvg], pch=16, col = "red")
m1 <- fit$mean[hvg]
v1 <- fit$var[hvg]
gg <- names(which(m1 >= 2 & v1 >= 3))
text(m1[gg] - 0.1, v1[gg] - 0.1, gg, col = "dodgerblue")

############################################################
## dimensionality reduction
############################################################
set.seed(1000)
sce <- denoisePCA(sce,
                  technical = res[, "tech"],
                  subset.row = hvg)

sce
ncol(reducedDim(sce, "PCA"))

set.seed(1000)
sce <- runTSNE(sce,
               use_dimred = "PCA",
               perplexity = 30)
sce
ncol(reducedDim(sce, "TSNE"))

tmp = umap(reducedDim(sce, "PCA"),
           verb = T,
           n_neighbors = 25)
reducedDim(sce, "UMAP") = tmp
sce

plotPCA(sce, ncomponents=3,
        colour_by="log10_total_features_by_counts")
plotTSNE(sce, colour_by="log10_total_features_by_counts")

############################################################
## shared snn graph based on top eigenvalues
############################################################
snn.gr       =  buildSNNGraph(sce, k=25, use.dimred="PCA")
clusters     =  igraph::cluster_walktrap(snn.gr, steps = 10)
sce$Cluster  = factor(factor(igraph::cut_at(clusters, 10)))

plotTSNE(sce, colour_by="Cluster")

## checking modularity per wf
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)
pheatmap(log.ratio,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100))

plotTSNE(sce, colour_by="cited1")
plotTSNE(sce, colour_by="six2")

############################################################
## save processed object
############################################################
saveRDS(sce, file = paste0("../dat/sce.rds"))

## procData_20190429.rds ## mix-up; was before we excluded doublets at this stage itself
## procData_20190424.rds
## procData_20190416.rds
############################################################
############################################################
