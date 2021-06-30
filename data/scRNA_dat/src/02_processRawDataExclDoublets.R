############################################################
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
fname <- "../cr_out/molecule_info.h5"
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

############################################################
## detect & remove empty droplets with FDR 0.01
############################################################
set.seed(1000) ## for permutation testing to calculate p-values
get.ed = emptyDrops(cnts, lower = 100)
is.cell <- get.ed$FDR <= 0.01
sum(is.cell, na.rm=TRUE) 

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
per.cell = perCellQCMetrics(sce, subsets = list(Mito = which(rowData(sce)$chr == "MT")))
colData(sce) <- cbind(colData(sce), per.cell)

par(mfrow=c(1,3))
hist(log10(sce$sum), breaks=20, col="grey80",
     xlab="Log-total UMI count")
hist(log10(sce$detected), breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
hist(sce$subsets_Mito_percent, breaks=20, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

############################################################
## cell filtering -- mark cells by QC metrics
############################################################
## no total library filter + total no. of expressed features + mito content
low.nexprs <- sce$detected < 1000 
table(low.nexprs) ## 1428

lmed = median(log(sce$subsets_Mito_percent))
lmad = mad(log(sce$subsets_Mito_percent))
lcut = lmed + 3*lmad
cut  = exp(lcut) ; rm(lmed,lmad,lcut)
high.mito <- sce$subsets_Mito_percent > cut
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
saveRDS(sce, file = "../dat/sce_noGF.rds")

############################################################
## gene filtering
############################################################
## genes with at least 3 counts in at least 3 cells
filtGenes <- function(obj){
    return(obj[(Matrix::rowSums(counts(obj) >= 3) >= 3), ])
}

sce <- filtGenes(sce) ; dim(sce)  ## 11155 x 4183

############################################################
## normalize and plot
############################################################
normalizeSCE <- function(obj){
    set.seed(1000)
    clusters = quickCluster(obj, method="igraph",
                            min.mean=0.1)
    obj = computeSumFactors(obj, min.mean = 0.1, cluster = clusters)
    obj = logNormCounts(obj)
    return(obj)
}

sce <- normalizeSCE(sce)

############################################################
## hvg
############################################################
fit <- trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
res <- decomposeVar(sce, fit = fit)
ii <- which((res$FDR < 0.01) & (pmax(0, res$bio)/res$total >= 0.25))
hvg  = sort(rownames(res)[ii]); length(hvg)  ## 826

############################################################
## dimensionality reduction
############################################################
set.seed(1000)
sce <- denoisePCA(sce,
                  technical = res[, "tech"],
                  subset.row = hvg)
set.seed(1000)
sce <- runTSNE(sce,
               dimred = "PCA",
               perplexity = 30)

############################################################
## shared snn graph based on top eigenvalues
############################################################
snn.gr       =  buildSNNGraph(sce, k=25, use.dimred="PCA")
clusters     =  igraph::cluster_walktrap(snn.gr, steps = 10)
sce$Cluster  = factor(factor(igraph::cut_at(clusters, 10)))

############################################################
## save processed object
############################################################
saveRDS(sce, file = paste0("../dat/sce.rds"))

############################################################
############################################################
