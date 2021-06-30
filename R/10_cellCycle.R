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
library(data.table)
library(ggsci)
library(colorspace)
library(openxlsx)
library(stringr)
library(pheatmap)

library(gridExtra)
library(cowplot)

source("./utilities_func.R")

library(Seurat)
library(stringr)

#===============================================
#- cell cycle effects
#===============================================
## run once
## s.genes <- cc.genes$s.genes
## g2m.genes <- cc.genes$g2m.genes
## ## seurat github; https://github.com/satijalab/seurat/issues/2493
## convertHumanGeneList <- function(x){
##     require("biomaRt")
##     human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
##     mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
##     genesV2 = getLDS(attributes = c("hgnc_symbol"),
##                      filters = "hgnc_symbol",
##                      values = x ,
##                      mart = human,
##                      attributesL = c("mgi_symbol"),
##                      martL = mouse,
##                      uniqueRows=T)
##     humanx <- unique(genesV2[, 2])
##     return(humanx)
## }
## s.genes <- convertHumanGeneList(cc.genes$s.genes)
## g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
## mus_cc.genes = list("s.genes" = s.genes,
##                     "g2m.genes" = g2m.genes)
## saveRDS(mus_cc.genes, file = "../data/external/mus_cc_genes.rds")

mus_cc.genes <- readRDS(file = "../data/external/mus_cc_genes.rds")
s.genes <- tolower(mus_cc.genes$s.genes)
g2m.genes <- tolower(mus_cc.genes$g2m.genes)

## cell cycle scoring with Seurat
sce = readRDS("../data/scRNA_dat/dat/sce_noGF.rds")
colnames(sce) <- sce$cell
counts = assay(sce, "counts")
seu <- CreateSeuratObject(counts = counts,
                          meta.data = as.data.frame(colData(sce)))

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst")
seu <- ScaleData(seu, features = rownames(seu))
seu@assays$RNA@meta.features = as.data.frame(rowData(sce))
seu <- CellCycleScoring(seu,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = F)
md = seu@meta.data
tmp = md[, c("S.Score", "G2M.Score", "Phase")]
colnames(tmp) = paste0("seucc_", c("S_score", "G2M_score", "phase"))

colData(sce) = cbind(colData(sce), tmp)
saveRDS(sce,"../data/scRNA_dat/dat/sce_noGFCC.rds")

sce0 = sce

## CC annotation to our final objects
sce = readRDS("../results/sce_fully-annotated.rds")
colData(sce) = cbind(colData(sce), colData(sce0)[colnames(sce),
                                                 c("seucc_S_score",
                                                   "seucc_G2M_score",
                                                   "seucc_phase")])
saveRDS(sce, "../results/sce_fully-annotated_CC.rds")

sce = readRDS("../results/sce_np_fully-annotated.rds")
colData(sce) = cbind(colData(sce),
                     colData(sce0)[colnames(sce), c("seucc_S_score",
                                                    "seucc_G2M_score",
                                                    "seucc_phase")])
saveRDS(sce, "../results/sce_np_fully-annotated_CC.rds")

#===============================================
# remove cell cycle effects - np focus
#===============================================
sce  <- readRDS("../results/sce_np_fully-annotated_CC.rds")
rD = rowData(sce)
hvg = rownames(rD)[which(rD$hvg)]

seu = as.Seurat(sce)
seu@assays$RNA@meta.features = as.data.frame(rD)
seu@assays$RNA@var.features = gsub("_", "-", hvg)

# regress on cell cycle scores
seu1 <- ScaleData(seu,
                  vars.to.regress = c("seucc_S_score", "seucc_G2M_score"),
                  features = rownames(seu))
set.seed(1000)
seu1 <- RunPCA(seu1,
               features = VariableFeatures(seu))
set.seed(1000)
seu1 <- RunTSNE(seu1,
                perplexity = 30)
saveRDS(seu1, "../results/seu_np_fully-annotated_CCReg1.rds")


#============================
#- self-renew vs primed after regressing out cell cycle effects
#============================
#
#- Supplemental Table 11

sce <- readRDS("../results/sce_np_fully-annotated_CC.rds")
cls = metadata(sce)$Cluster_colors
ii   = sce$cluster_tme %in% c("self-renew","primed")
sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)

lfc = log2(1.2)
fdrCut = 0.2
pD = colData(sce2)
scores = as.matrix(pD[, c("seucc_S_score", "seucc_G2M_score")])
design <- model.matrix(~ scores)
mkrs2 = findMarkers(sce2,
                    sce2$cluster_tme,
                    design = design,
                    pval.type="all",
                    direction="any",
                    lfc=lfc,
                    assay.type = "logcounts_SAVER")
mkrs = mkrs2
ii       = mkrs[[1]]$FDR < fdrCut
gg       = rownames(mkrs[[1]])[ii]
gg.sr.pr = gg
gg.sr.pr.up = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed < 0]
gg.sr.pr.dn = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed > 0]
saveRDS(mkrs, file="../results/np_self-vs-primed-DEG_CCScores_lfc1.2.rds")

dat = as.data.frame(mkrs[[1]][ii,])
dat$gene = rownames(dat)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=dat, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/supTab11_np_self-vs-primed-DEG_CCScores_lfc1.2.xlsx')


#- Supplemental Table 12

sce <- readRDS("../results/sce_np_fully-annotated_CC.rds")
cls = metadata(sce)$Cluster_colors
ii   = sce$cluster_tme %in% c("self-renew","primed")
sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)

## excl CC associated genes
diff <- getVarianceExplained(sce2, DataFrame(sce2$seucc_phase))
discard <- diff > 5
summary(discard) ## 518
discard.g = rownames(discard)[discard]

sce2_subset = sce2[(!rownames(sce2) %in% discard.g), ]

lfc = log2(1.2)
fdrCut = 0.2
mkrs31 = findMarkers(sce2_subset,
                     sce2_subset$cluster_tme,
                     pval.type="all",
                     direction="any",
                     lfc=lfc,
                     assay.type = "logcounts_SAVER")
mkrs = mkrs31
ii       = mkrs[[1]]$FDR < fdrCut
gg       = rownames(mkrs[[1]])[ii]
gg.sr.pr = gg
gg.sr.pr.up = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed < 0]
gg.sr.pr.dn = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed > 0]
saveRDS(mkrs, file="../results/np_self-vs-primed-DEG_CCRmPhaseGenes.rds")

dat = as.data.frame(mkrs[[1]][ii,])
dat$gene = rownames(dat)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=dat, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/suppTab12_np_self-vs-primed-DEG_CCRmPhaseGenes.xlsx')


#============================
#- birc5 correlating genes in distal tubular and ub/cd
#============================
#
#- Supplemental Table 13

sce_wk <- readRDS("../results/sce_fully-annotated_CC.rds")
sce_np <- readRDS("../results/sce_np_fully-annotated_CC.rds")

n1 = colnames(sce_np)[sce_np$cluster_tme %in% c("i_distal_tubular", "m_distal_tubular")]
n2 = colnames(sce_wk)[sce_wk$cluster_it == "ureteric_bud/collecting_duct"]
sce = sce_wk[,c(n1,n2)]

mat = logcounts(sce)

set.seed(100)
var.cor <- correlatePairs(sce)
sig.cor <- var.cor$FDR <= 0.05

var.cor1 = var.cor[which(var.cor$gene1 == "birc5"), ]
var.cor1 = var.cor1[order(var.cor1$rho, decreasing = T), ]
sig.var.cor1 = var.cor1[which(var.cor1$FDR <= 0.05), ]
sig.var.cor1 

## writing out
df = data.frame(sig.var.cor1)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=df, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/supTab13_birc5CorGenes_dist_tub_ub.xlsx')

#=========================
#- Cell cycle effects
#=========================

#- Supplemental Figure 7

#======================
#- TSNE of NP CELLCYCLE assignments
#======================
#

sce <- readRDS("../results/sce_np_fully-annotated_CC.rds")
cls = metadata(sce)$Cluster_colors
pD = colData(sce)

dat = data.frame(reducedDim(sce,"TSNE"))
colnames(dat) = c("tSNE1","tSNE2")
dat$cluster = sce$cluster_tme
dat = cbind(dat, as.data.frame(pD[, c("seucc_phase", "seucc_S_score", "seucc_G2M_score")]))

#- median cell for each cluster
cfu <- function(cname,dat){
  tmp = dat[dat$cluster==cname,]
  t1  = tmp[,"tSNE1"] ; m1 = median(t1)
  t2  = tmp[,"tSNE2"] ; m2 = median(t2)
  ii  = which.min((t1-m1)^2 + (t2-m2)^2)
  return(c(t1[ii],t2[ii]))
}
labdat = sapply(levels(dat$cluster),cfu,dat)
labdat = data.frame(t(labdat))
labdat$label = rownames(labdat)
colnames(labdat) = c("tSNE1","tSNE2","cluster")


cls_tmp = cls[names(cls) %in% dat$cluster]
ggplot(dat,aes(x=tSNE1,y=tSNE2,col=cluster)) +
    geom_point(alpha=.6,size=1) +
    scale_color_manual(values = cls_tmp) +
    geom_text_repel(  data = subset(labdat,tSNE1<2.5),
                    nudge_x = -23  -subset(labdat,tSNE1<2.5)$tSNE1,
                    aes(label=cluster,fill=cluster),
                    alpha=.5,
                    segment.size=1/2,
                    segment.alpha=1/3,
                    fontface = 'bold',
                    color = 'black',
                    direction="y",
                    hjust=1) +
  geom_text_repel(  data = subset(labdat,tSNE1>2.5),
                    nudge_x = 22  -subset(labdat,tSNE1>2.5)$tSNE1,
                    aes(label=cluster,fill=cluster),
                    alpha=.5,
                    segment.size=1/2,
                    segment.alpha=1/3,
                    fontface = 'bold',
                    color = 'black',
                    direction="y",
                    hjust=0) +
  scale_x_continuous( breaks = c(-20,-10,0,10,20),
                      limits=c(-35,35)) +
  theme_minimal() +
  theme( panel.background = element_blank(), panel.border = element_blank(),
         panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
         legend.position= "none")

p2 <- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=seucc_phase)) +
    geom_point(alpha=.6,size=.8) 
p4 <- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=seucc_G2M_score)) +
    geom_point(alpha=.6,size=.8) 
p5 <- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=seucc_S_score)) +
    geom_point(alpha=.6,size=.8) 


pdf(file="../figures/supFig7_wk_np_CCPhase-tsne.pdf",width=6,height=4.5)
p2
dev.off()
pdf(file="../figures/supFig7_wk_np_CCG2MScore-tsne.pdf",width=6,height=4.5)
p4
dev.off()
pdf(file="../figures/supFig7_wk_np_CCSScore-tsne.pdf",width=6,height=4.5)
p5
dev.off()


#======================
#- TSNE after regressing out CC effects
#======================
#

seu <- readRDS("../results/seu_np_fully-annotated_CCReg1.rds")

p1 <- DimPlot(seu, reduction = "tsne", group.by = "cluster_tme") +
    ggtitle("NP cells after regressing out CC effects using classification scores (seurat)")
p2 <- DimPlot(seu, reduction = "tsne", group.by = "seucc_phase") +
    ggtitle("NP cells after regressing out CC effects using classification scores (seurat)")
pp = plot_grid(p1, p2, ncol = 2)
figFile = paste0("../figures/supFig7_wk_np_CCScoresReg1-tsne.pdf")
save_plot(filename = figFile,
          pp,
          ncol = 2,
          nrow = 1,
          base_width = 7,
          base_height = 7)


#- end
