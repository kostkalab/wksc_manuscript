library(scater)
library(scran)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(dendsort)
library(ggsci)
library(RColorBrewer)
library(colorspace)
library(ComplexHeatmap)
library(slingshot)
library(cowplot)

#=========================
#- WK marker genes heatmap
#=========================
#
#- Supplemental Figure 1

sce = readRDS("../results/sce_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors

#- NOTE: the talbe have "some" not "all" for pval.type.
fm = findMarkers( sce,
                  sce$cluster_it,
                  pval.type="all",
                  direction="any",
                  lfc=log2(1.5))

gl = lapply(fm,function(x) rownames(x[1:15,]))
glu = unique(unlist(gl))

nce = 45
inds = tapply(1:ncol(sce),sce$cluster_it,sample,nce)
inds = (unlist(inds))
clus = sce$cluster_it[inds]
mat  = as.matrix(logcounts(sce)[glu,inds])

dd     = mat - rowMeans(mat)
dd     = dd / apply(dd,1,sd)

source("./utilities_func.R")
ord = 1:ncol(mat)

dst = dist((dd))
hc  = hclust(dst,method="ward.D2")
dend = as.dendrogram(hc)
dend = dendextend::seriate_dendrogram( dend,dst)

colfun = colorRampPalette((diverge_hsv(n = 7)))(100)
dd[dd>3] = 3 ; dd[dd < -3] = -3

hm <- Heatmap(dd[,ord],name="expression",col=colfun,cluster_columns=FALSE,
column_labels=rep("",ncol(dd)), cluster_rows=dend,row_split=10,
column_split= clus,use_raster=TRUE,row_names_gp = gpar(fontsize = 8),
show_row_dend=FALSE,row_title="",column_title_rot=90)

pdf("../figures/wk_celltype-DEG-heatmap.pdf",width=6,height=16)
draw(hm)
dev.off()

#=========================
#- NP lineages TSNE
#=========================
#
#- Supplemental Figure 2


sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
tt = cbind(sce$slingPseudotime_1,sce$slingPseudotime_2,sce$slingPseudotime_3)
sce$slingNL = rowSums(!is.na(tt))

sds = SlingshotDataSet(sce)

sce$Cluster6 = forcats::fct_relevel(sce$Cluster6,c("self-renew","primed",
"differentiating","i_tubular_dist","m_tubular_dist","i_tubular_prox",
"m_tubular_prox","i_podocytes" ,"m_podocytes"))

pl = plotTSNE(sce,col="cluster_tme",shape_by='slingNL') +
      scale_color_manual(values = cls) + scale_fill_manual(values=cls)

d1 = data.frame(sds@curves[[1]]$s)
colnames(d1) = c("x","y")
d2 = data.frame(sds@curves[[2]]$s)
colnames(d2) = c("x","y")
d3 = data.frame(sds@curves[[3]]$s)
colnames(d3) = c("x","y")

pl = pl + geom_line(aes(x=x,y=y),data=d1,size=2,alpha=.5) +
     geom_line(aes(x=x,y=y),data=d2,size=2,alpha=.5) +
     geom_line(aes(x=x,y=y),data=d3,size=2,alpha=.5) +
     xlab("tSNE1") + ylab("tSNE2")

save_plot("../figures/np_pseudotime-lineages.pdf",pl,base_height = 3.71, base_asp=1.618)

#=========================
#- NP marker genes heatmap
#=========================
#
#- Supplemental Figures 3 and 4

sce = readRDS("../results/sce_np_fully-annotated.rds")
fm = findMarkers( sce,
                  sce$cluster_tme,
                  pval.type="any",
                  direction="any",
                  lfc=log2(1.5),            #- this is base2 because log-counts are; so OK
                  full.stats = TRUE)

hmfu <- function(sce, gens, look_cls, nrowsep=5){

  mat = logcounts(sce[gns,sce$cluster_tme %in% look_cls])
  tcls = droplevels(sce$cluster_tme[sce$cluster_tme %in% look_cls])
  tcls = forcats::fct_relevel(tcls,look_cls)
  ncells = 45 #min(table(tcls))
  inds = tapply(1:ncol(mat),tcls,sample,ncells)
  inds = (unlist(inds))
  clus = tcls[inds]
  mat  = as.matrix(mat[,inds])

  dd     = mat - rowMeans(mat)
  dd     = dd / apply(dd,1,sd)

  source("./utilities_func.R")
  ord = 1:ncol(mat)

  dst = dist((dd))
  hc  = hclust(dst,method="ward.D2")
  dend = as.dendrogram(hc)
  dend = dendextend::seriate_dendrogram( dend,dst)

  colfun = colorRampPalette((diverge_hsv(n = 7)))(100)
  dd[dd>3] = 3 ; dd[dd < -3] = -3

  rownames(dd) <- sub("_ENSMUSG.*","",rownames(dd))

  hm <- Heatmap(dd[,ord],name="expression",col=colfun,cluster_columns=FALSE,
  column_labels=rep("",ncol(dd)), cluster_rows=dend, row_split=nrowsep,
  column_split= clus,use_raster=TRUE,row_names_gp = gpar(fontsize = 8),
  show_row_dend=FALSE,row_title="", column_title_rot=90)

  return(hm)
}

look_cls = c("differentiating","i_podocyte","m_podocyte","i_proximal_tubular",
              "m_proximal_tubular","i_distal_tubular","m_distal_tubular")

tst = sapply(look_cls,function(x) fm[[x]][fm[[x]]$Top<=10,])
gns = sort(unique(unlist(lapply(tst,rownames))))

hm <- hmfu(sce,gns,look_cls,nrowsep=6)

pdf("../figures/np_celltype-DEG-heatmap.pdf",width=6,height=15) #- FS 3
draw(hm)
dev.off()

#- TUBULAR DISTAL VS PROXIMAL (SF 4A)
#============================

gns      <- rownames(fm[["m_distal_tubular"]])[which(fm[["m_distal_tubular"]][,"stats.m_proximal_tubular"][,"log.FDR"] < log(0.1))]
gns      <- gns[1:100]
look_cls <- c("i_proximal_tubular","m_proximal_tubular","i_distal_tubular","m_distal_tubular")
hm       <- hmfu(sce,gns,look_cls,nrowsep=4)

pdf("../figures/np_celltype-tp-vs-td-DEG-heatmap.pdf",width=6,height=12)
draw(hm)
dev.off()

#- TUBULAR DISTAL MATURE vs IMMATURE (SF 4B)
#====================================

gns <- rownames(fm[["m_distal_tubular"]])[which(fm[["m_distal_tubular"]][,"stats.i_distal_tubular"][,"log.FDR"] < log(0.1))]
look_cls = c("i_distal_tubular","m_distal_tubular")
hm <- hmfu(sce,gns,look_cls,nrowsep=2)
pdf("../figures/np_celltype-mtd-vs-itd-DEG-heatmap.pdf",width=6,height=6)
draw(hm)
dev.off()

#- TUBULAR POXIMAL MATURE vs IMMATURE (SF 4C)
#====================================

gns <- rownames(fm[["m_proximal_tubular"]])[which(fm[["m_proximal_tubular"]][,"stats.i_proximal_tubular"][,"log.FDR"] < log(0.1))]
look_cls = c("i_proximal_tubular","m_proximal_tubular")
hm <- hmfu(sce,gns,look_cls,nrowsep=2)
pdf("../figures/np_celltype-mtp-vs-itp-DEG-heatmap.pdf",width=6,height=6)
draw(hm)
dev.off()

#- PODOCYTES MATURE vs IMMATURE (SF 4D)
#==============================

gns <- rownames(fm[["m_podocyte"]])[which(fm[["m_podocyte"]][,"stats.i_podocyte"][,"log.FDR"] < log(0.1))]
gns <- gns[1:100]
look_cls = c("i_podocyte","m_podocyte")
hm <- hmfu(sce,gns,look_cls,nrowsep=2)
pdf("../figures/np_celltype-mpod-vs-ipod-DEG-heatmap.pdf",width=6,height=12)
draw(hm)
dev.off()

#=======================
#- BIRC5 TUB and UB
#- Supplemental Figure 5
#=======================

sce_wk <- readRDS("../results/sce_fully-annotated.rds")
sce_np <- readRDS("../results/sce_np_fully-annotated.rds")

n1 = colnames(sce_np)
n2 = colnames(sce_wk)[sce_wk$cluster_it == "ureteric_bud/collecting_duct"]
sce = sce_wk[,c(n1,n2)]

mgv               = modelGeneVar(sce)
ii                = mgv$FDR <= .1
rowData(sce)$hvg = ii

set.seed(1000)
sce <- denoisePCA(sce,
                  technical = mgv$tech,
                  subset.row = rownames(sce)[which(rowData(sce)$hvg)])

set.seed(1000)
sce <- runTSNE(sce,
              use_dimred = "PCA")

cls = metadata(sce)$Cluster_colors

dat = data.frame(reducedDim(sce,"TSNE"))
colnames(dat) = c("tSNE1","tSNE2")
dat$cluster = sce$cluster_tme
dat$birc5 = logcounts(sce)["birc5",]
cls_tmp = cls[names(cls) %in% dat$cluster]

p1 <- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=cluster)) +
  geom_point(alpha=.8,size=.8) +
  scale_x_continuous( breaks = c(-20,-10,0,10,20),
                      limits=c(-30,30)) +
  theme_minimal() +
  theme( panel.background = element_blank(),
         panel.border = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
  scale_color_manual(values = cls_tmp)


p2<- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=birc5)) +
  geom_point(alpha=.8,size=.8) +
  scale_x_continuous( breaks = c(-20,-10,0,10,20),
                      limits=c(-30,30)) +
  theme_minimal() +
  theme( panel.background = element_blank(),
         panel.border = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
         scale_color_gradient(low="steelblue",high="firebrick")


 p3 <- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=cluster)) +
   geom_point(alpha=.8,size=.8) +
   scale_x_continuous( breaks = c(-20,-10,0,10,20),
                       limits=c(-30,30)) +
   theme_minimal() +
   theme( panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position= "none")+
   scale_color_manual(values = cls_tmp)


 p4<- ggplot(dat,aes(x=tSNE1,y=tSNE2,col=birc5)) +
   geom_point(alpha=.8,size=.8) +
   scale_x_continuous( breaks = c(-20,-10,0,10,20),
                       limits=c(-30,30)) +
   theme_minimal() +
   theme( panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position= "none")+
          scale_color_gradient(low="steelblue",high="firebrick")



save_plot("../figures/np-ub_tsne-clust.pdf",p1)
save_plot("../figures/np-ub_tsne-birc5.pdf",p2)

save_plot("../figures/np-ub_tsne-noleg-clust.pdf",p3)
save_plot("../figures/np-ub_tsne-noleg-birc5.pdf",p4)
