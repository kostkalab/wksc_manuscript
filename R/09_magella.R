library(data.table)
library(dplyr)
library(scater)
library(scran)
library(AUCell)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(stringr)

annot.df <- readRDS("../data/scRNA_dat/dat/annot.rds")


#===============================================
#- MAGELLA et al. data
#===============================================
## download: GSM2796989	(GSE104396); syn11001759
a <- fread("../data/external/magella/GSM2796989_E14_5_WT_kidney_10X_matrix.txt.gz")
counts = as(as.matrix(a[, 2:ncol(a)]), "dgCMatrix")
genes = a$UID; rm(a)

a3 <- fread("../data/external/magella/syn11027925/knnClassified16populations/groups.10X-log2-NearestNeighbor.txt", stringsAsFactors = F, head = F)
dim(a3) ## 2161 x 3
a3.df = data.frame(a3)
rownames(a3.df) = a3$V1

a4 <- fread("../data/external/magella/syn11027925/knnClassified16populations/exp.10X-log2-NearestNeighbor.txt", stringsAsFactors = F, head = T)
dim(a4) ## 27998 x 2162

## get processed obj for magella 10X data; 
sce = SingleCellExperiment(assays = list(counts = counts))
rowData(sce)$ID = annot.df$ensembl_gene_id
rowData(sce)$symbol = genes
colData(sce)$cell = as.character(colnames(counts))
rowData(sce)$chr <- annot.df$chrom
rowData(sce)$symbol <- annot.df$gene_name
rowData(sce)$symbol = tolower(rowData(sce)$symbol)
rowData(sce)$feature <- annot.df$feature
rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$symbol)

## restrict chrs and cells
sce = sce[rowData(sce)$chr %in% c(as.character(1:22),"X","Y","MT"),]; dim(sce)
counts = as.matrix(assay(sce, "counts"))
rownames(counts) = rownames(sce)
sce = sce[, colnames(a4)[2:ncol(a4)]]; dim(sce) ## 27921 x 2161
## working with knn class results; add cell cluster label from knn class
pD = colData(sce)

pD$knn_cluster = paste0("c", a3.df[rownames(pD), "V2"])
tmp = as.character(pD$knn_cluster)
tmp[tmp == "c1"] = "medullary_collecting_duct"
tmp[tmp == "c2"] = "collecting_duct"
tmp[tmp == "c3"] = "ureteric_tip"
tmp[tmp == "c4"] = "loop_of_henle"
tmp[tmp == "c5"] = "distal_comma_shaped_body"
tmp[tmp == "c6"] = "podocytes"
tmp[tmp == "c7"] = "mid_s_shaped_body"
tmp[tmp == "c8"] = "proximal_tubular"
tmp[tmp == "c9"] = "pretubular_aggregate"
tmp[tmp %in% c("c10", "c11", "c12")] = "cap_mesenchyme"
tmp[tmp == "c13"] = "endothelium"
tmp[tmp == "c14"] = "nephrogenic_zone_stroma"
tmp[tmp == "c15"] = "cortical_stroma"
tmp[tmp == "c16"] = "medullary_stroma"
pD$knn_celltype = factor(tmp, levels = unique(tmp))
colData(sce) = pD

#- rm very low expressed genes in this subset of cells
#=====================================================
filtGenes <- function(obj){
    return(obj[(Matrix::rowSums(counts(obj) >= 3) >= 3), ])
}
sce <- filtGenes(sce); dim(sce) ## 7223 x 2161

normalizeSCE <- function(obj){
    set.seed(1000)
    clusters = quickCluster(obj,
                            method = "igraph",
                            min.mean = 0.1)
    obj = computeSumFactors(obj, min.mean = 0.1, cluster = clusters)
    obj = logNormCounts(obj)
    return(obj)
}

sce <- normalizeSCE(sce)
## get hvg
mgv               = modelGeneVar(sce)
ii                = mgv$FDR <= .1
rowData(sce)$hvg = ii

## dimensionality reduction
set.seed(1000)
sce <- denoisePCA(sce,
                  technical = mgv$tech,
                  subset.row = rowData(sce)$hvg)
set.seed(1000)
sce <- runTSNE(sce,
               dimred = "PCA")
set.seed(1000)
sce <- runUMAP(sce,
               dimred = "PCA",
               n_neighbors = 25)

#- transfer labels from wk to magella
#---------------------------------
sce_mag = sce
pD = colData(sce_mag)
emat = logcounts(sce_mag)
cells_rankings <- AUCell_buildRankings(emat)

## whole kidney
sce_wk = readRDS("../results/sce_fully-annotated.rds")
fm = findMarkers(sce_wk,
                 sce_wk$cluster_it,
                 pval.type = "all",
                 direction = "up",
                 lfc = log2(1.5))

gl = lapply(fm, function(x) rownames(x[1:20,]))
auc <- AUCell_calcAUC(gl, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
tmp = t(assay(auc))
colnames(tmp) = paste0(gsub("/", "_",colnames(tmp)), "_wktop20_aucsc")
pD = cbind(pD, tmp)

## annotate cell type by max score
dat = assay(auc)
metadata(sce_mag)$wk20_auc = dat
datNew = (dat - rowMeans(dat))/apply(dat, 1, sd)
new_annot = apply(datNew, 2, function(x){
    return(rownames(datNew)[which.max(x)])
})
pD$wktop20_annot = new_annot

## np lineage cell types 
sce_np = readRDS("../results/sce_np_fully-annotated.rds")
fm = findMarkers(sce_np,
                 sce_np$cluster_tme,
                 pval.type = "all",
                 direction = "up",
                 lfc = log2(1.5))
gl = lapply(fm, function(x) rownames(x[1:20,]))
auc <- AUCell_calcAUC(gl, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
tmp = t(assay(auc))
colnames(tmp) = paste0(gsub("/", "_",colnames(tmp)), "_nptop20_aucsc")
pD = cbind(pD, tmp)
## annotate cell type by max score
dat = assay(auc)
metadata(sce_mag)$np20_auc = dat
datNew = (dat - rowMeans(dat))/apply(dat, 1, sd)
new_annot = apply(datNew, 2, function(x){
    return(rownames(datNew)[which.max(x)])
})
pD$nptop20_annot = new_annot
colData(sce_mag) = pD

saveRDS(sce_mag, file = "../results/sce_magella.rds")


#=========================
#- Magella et al data
#=========================
#
#- Supplemental Figure 6A

sce = readRDS("../results/sce_magella.rds")
tab1 = table(sce_mag$knn_celltype, sce_mag$wktop20_annot)
## percentage of orig annotations that got assigned to a cell type
tt1 = tab1
tt1 = round(tt1/rowSums(tt1), 2)
mat <- matrix(as.numeric(tt1), nrow = nrow(tt1), ncol = ncol(tt1))
rownames(mat) = rownames(tt1)
colnames(mat) = colnames(tt1)
ord = c(7,6,2,9,8,1,4,5,3,10)

p1 <- pheatmap(mat[,ord],
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               fontsize = 8,
               fontsize_row = 10,
               fontsize_col = 10,
               color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
               border = "white",
               cellwidth = 12,
               cellheight = 12)

pdf("../figures/supFig6a_magella-comp.pdf", width = 8, height = 8)
p1
dev.off()


#=========================
# birc5, np
#=========================

#- Supplemental Figure 6B
sce_mag = readRDS("../results/sce_magella.rds")
look_cls = c("loop_of_henle",
             "distal_comma_shaped_body",
             "podocytes",
             "mid_s_shaped_body",
             "proximal_tubular",
             "pretubular_aggregate",
             "cap_mesenchyme")

sce_mag_np = sce_mag[, sce_mag$knn_celltype %in% look_cls]
tcls <- droplevels(sce_mag_np$knn_celltype[sce_mag_np$knn_celltype %in% look_cls])
tcls = forcats::fct_relevel(tcls,look_cls)
sce_mag_np$knn_celltype = tcls
dim(sce_mag_np) ## 7223 x 692
sce = sce_mag_np

sce = sce[rowSums(counts(sce)>=2) >= 3, ]
mgv               = modelGeneVar(sce)
ii                = mgv$FDR <= .1
set.seed(1000)
sce <- denoisePCA(sce,
	technical = mgv$tech,
	subset.row = ii)
set.seed(1000)
sce <- runTSNE(sce,
	 dimred = "PCA")

#- Birc5 and different types of NP cells:
p1 <- plotTSNE(sce, col = "knn_celltype")
p2 <- plotTSNE(sce, col = "birc5") +
    guides(fill = guide_colorbar(title = "Birc5", label = FALSE, ticks = FALSE))

sfu <- function(x) {
    x[x >= quantile(x,.99)] = quantile(x,.99)
    x = x - min(x)
    x = x/max(x)
    x
}

sce$self_renew <- sce$'self-renew_nptop20_aucsc'  %>% sfu() 
p3 <- plotTSNE(sce,col="self_renew") +
    guides(fill = guide_colorbar(title="self renew",label=FALSE, ticks=FALSE))

#plot_grid(p1, p2, p3, nrow = 1)

p_leg   <- plot_grid(p1,p2,p3, nrow=1)
p_noLeg <- plot_grid(p1 + theme(legend.position = "none"),
                     p2 + theme(legend.position = "none") ,
                     p3 + theme(legend.position = "none"),
                     nrow = 1)
save_plot(p_leg,
          file="../figures/supFig6b_magella_tsne-np-birc5_leg.pdf", base_height=4, base_aspect_ratio=3)
save_plot(p_noLeg,
          file="../figures/supFig6b_magella_tsne-np-birc5_noLeg.pdf", base_height=4, base_aspect_ratio=3)


#=========================
# aucell scores
#=========================

#- Supplemental Figure 6C
sce$'self_renew'      <- sce$'self-renew_nptop20_aucsc'         %>% sfu() 
sce$'primed'          <- sce$'primed_nptop20_aucsc'             %>% sfu()
sce$'differentiating' <- sce$'differentiating_nptop20_aucsc'    %>% sfu()
sce$'podocyte'        <- sce$'m_podocyte_nptop20_aucsc'         %>% sfu()
sce$'tub_prox'        <- sce$'m_proximal_tubular_nptop20_aucsc' %>% sfu()
sce$'tub_dist'        <- sce$'m_distal_tubular_nptop20_aucsc'   %>% sfu()

p4 <- plotTSNE(sce,col="podocyte") +
    guides(fill = guide_colorbar(title="podocyte",label=FALSE, ticks=FALSE))
p5 <- plotTSNE(sce,col="tub_prox") +
    guides(fill = guide_colorbar(title="tub_prox",label=FALSE, ticks=FALSE))
p6 <- plotTSNE(sce,col="tub_dist") +
    guides(fill = guide_colorbar(title="tub_dist",label=FALSE, ticks=FALSE))

p_leg   <- plot_grid(p4,p5,p6, nrow=1)
p_noLeg <- plot_grid(p4 + theme(legend.position = "none"),
                     p5 + theme(legend.position = "none"),
                     p6 + theme(legend.position = "none"),
                     nrow=1)

save_plot(p_leg,
          file="../figures/supFig6c_magella_tsne-np-celltypes_leg.pdf", base_height=4, base_aspect_ratio=3)
save_plot(p_noLeg,
          file="../figures/supFig6c_magella_tsne-np-celltypes_noLeg.pdf", base_height=4, base_aspect_ratio=3)

#- end

