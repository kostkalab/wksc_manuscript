library(scater)
library(scran)
library(uwot)
library(reshape2)
library(slingshot)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(dendsort)
library(scds)
library(data.table)
library(SAVER)
library(ggsci)
library(RColorBrewer)
library(colorspace)

#- our SCE
sce = readRDS("../results/sce_annotated.rds")


#- Reduce cells in sce (NP-type-cells only)
#==========================================

sce$selected = sce$Cluster %in% c("nephron-progenitor",
                                  "mixed/differentiating",
                                  "tubular",
                                  "podocytes")

sce2         = sce[, sce$selected]; dim(sce2) ## 11181 x 1272

#- rm very low expressed genes in this subset of cells
#=====================================================

filtGenes <- function(obj){
    return(obj[(Matrix::rowSums(counts(obj) >= 3) >= 3), ])
}
sce2 <- filtGenes(sce2); dim(sce2) ## 9611 x 1273

#- pull out cells and keep old projections
reducedDim(sce2, "PCA_old")  = reducedDim(sce, "PCA")[sce$selected, ]
reducedDim(sce2, "TSNE_old") = reducedDim(sce, "TSNE")[sce$selected, ]


#- Find variable genes (again)
#=============================

mgv               = modelGeneVar(sce2)
ii                = mgv$FDR <= .1
rowData(sce2)$hvg = ii

#- RE-calculate projections
#==========================

set.seed(1000)
sce2 <- denoisePCA(sce2,
                   technical = mgv$tech,
                   subset.row = rowData(sce2)$hvg)
set.seed(1000)
sce2 <- runTSNE(sce2,
                use_dimred = "PCA")
set.seed(1000)
sce2 <- runUMAP(sce2,
                dimred = "PCA")

saveRDS(sce2,file="../results/sce_np.rds")


#- end
