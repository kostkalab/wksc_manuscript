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
library(ggsci)
library(RColorBrewer)
library(colorspace)
library(openxlsx)
library(stringr)
library(pheatmap)

#======================
#- TSNE of WK CELLTYPES
#======================
#
#- FIGURE 1B

sce = readRDS("../results/sce_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
dat = data.frame(reducedDim(sce,"TSNE"))
colnames(dat) = c("tSNE1","tSNE2")
dat$cluster = sce$cluster_it

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

pdf(file="../figures/fig1b_wk_celltype-tsne.pdf",width=8,height=4.5)

cls = metadata(sce)$Cluster_colors
cls_tmp = cls[names(cls) %in% dat$cluster]
ggplot(dat,aes(x=tSNE1,y=tSNE2,col=cluster)) +
  geom_point(alpha=.6,size=.8) +
  scale_color_manual(values = cls_tmp) +
  geom_text_repel(  data = subset(labdat,tSNE1<0),
                    nudge_x = -45  -subset(labdat,tSNE1<0)$tSNE1,
                    aes(label=cluster,fill=cluster),
                    alpha=.5,
                    segment.size=1/2,
                    segment.alpha=1/3,
                    fontface = 'bold',
                    color = 'black',
                    direction="y",
                    hjust=1) +
  geom_text_repel(  data = subset(labdat,tSNE1>0),
                    nudge_x = 50  -subset(labdat,tSNE1>0)$tSNE1,
                    aes(label=cluster,fill=cluster),
                    alpha=.5,
                    segment.size=1/2,
                    segment.alpha=1/3,
                    fontface = 'bold',
                    color = 'black',
                    direction="y",
                    hjust=0) +
  scale_x_continuous( breaks = c(-40,-30,-20,-10,0,10,20,30,40),
                      limits=c(-45,50)) +
  theme_minimal() +
  theme( panel.background = element_blank(), panel.border = element_blank(),
         panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
         legend.position= "none")
dev.off()

#======================
#- TSNE of NP CELLTYPES
#======================
#
#- FIGURE 2A

sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
dat = data.frame(reducedDim(sce,"TSNE"))
colnames(dat) = c("tSNE1","tSNE2")
dat$cluster = sce$cluster_tme

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

pdf(file="../figures/fig2a_np_celltype-tsne.pdf",width=8,height=4.5)

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
dev.off()


#=========================
#- VIOLINS of WK CELLTYPES
#=========================
#
#- FIGURE 1C

sce = readRDS("../results/sce_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
gens = tolower(c("Cited1","Six2"  , #- NPs
                 "Lhx1", "Pax8",    #- mixed / differentiating
                 "Podxl" ,"Nphs1" , #- podocytes
                 "Pdzk1", "Slc34a1" ,
                 "Tmem52b", "Shd",
                 "Calb1" ,"Gata3" , #- ureteric bud / collecting duct (UB/CD)
                 "Col1a1","meis1" , #- stroma 
                 "Aldh1a2", 
                 "Dlk1", "Postn",   #- cortical / mesangial
                 "Cldn11", 
                 "Emcn"  ,"Kdr"   , #- endothelial
                 "Cd52"  ,"Fcer1g"))


df            = data.frame(t(as.matrix(logcounts(sce)[gens,])))
colnames(df)  = str_to_title(colnames(df))
df$cluster    = sce$cluster_it
dfm           = reshape2::melt(df)
colnames(dfm) = c("cluster","gene","expression")


pdf(file="../figures/fig1c_wk_celltype-violins.pdf",width=14.5,height=5)
cls_tmp = cls[names(cls) %in% df$cluster]
ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
   geom_violin(scale="width") +
   coord_flip() +
   facet_wrap(.~gene,nrow=1) +
   theme_minimal() +
   theme(axis.text.x = element_blank()     ,
         panel.spacing.x = unit(0.5, "lines"),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         legend.position = "none",
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text.x = element_text(angle = 45, hjust = 0, vjust=0)
         )+
         scale_color_manual(values = cls_tmp) +
         scale_fill_manual(values=cls_tmp) +
         guides(fill = guide_legend(reverse=TRUE)) +
         guides(color = guide_legend(reverse=TRUE))
dev.off()

#=========================
#- VIOLINS of NP CELLTYPES
#=========================
#
#- FIGURE 2B

sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
gens = tolower(c("cited1" ,  "six2"   , "pax8",
                 "smc2"   ,  "smc4"   ,
                 "cdca3"  ,  "birc5"  ,
                 "lhx1"  ,
                 "podxl" ,"nphs1",
                 "pdzk1"  ,  "slc34a1",
                 "tmem52b",  "shd"))

df            = data.frame(t(as.matrix(logcounts(sce)[gens,])))
colnames(df)  = str_to_title(colnames(df))
df$cluster    = sce$cluster_tme
dfm           = reshape2::melt(df)
colnames(dfm) = c("cluster","gene","expression")

pw = 14/22*14.5
ph = 9/11*5
pdf(file="../figures/fig2b_np_celltype-violins.pdf",width=pw,height=ph)
cls_tmp = cls[names(cls) %in% df$cluster]
ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
   geom_violin(scale="width") +
   coord_flip() +
   facet_wrap(.~gene,nrow=1) +
   theme_minimal() +
   theme(axis.text.x = element_blank()     ,
         panel.spacing.x = unit(0.5, "lines"),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         legend.position = "none",
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text.x = element_text(angle = 45, hjust = 0, vjust=0)
         )+
         scale_color_manual(values = cls_tmp) +
         scale_fill_manual(values=cls_tmp) +
         guides(fill = guide_legend(reverse=TRUE)) +
         guides(color = guide_legend(reverse=TRUE))
dev.off()

#===========================
#- VIOLINS of POD  CELLTYPES
#===========================
#
#- FIGURE 4B

sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors

gens = tolower(c("cdkn1c","gas1","wt1","robo2","sparc","fxyd6",
"mafb","cldn5","plat","slc9a3r2","tcf21","synpo","nphs2", "magi2",
"pdlim2","pard3b","ptpro"))


df            = data.frame(t(as.matrix(logcounts(sce)[gens,])))
colnames(df)  = str_to_title(colnames(df))
df$cluster    = sce$cluster_tme
dfm           = reshape2::melt(df)
colnames(dfm) = c("cluster","gene","expression")

pw = 17/22*14.5
ph = 9/11*5

pdf(file="../figures/fig4b_np_podocyte-violins.pdf",width=pw,height=ph)
cls_tmp = cls[names(cls) %in% df$cluster]
ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
  geom_violin(scale="width") +
  coord_flip() +
  facet_wrap(.~gene,nrow=1) +
  theme_minimal() +
  theme(axis.text.x = element_blank()     ,
        panel.spacing.x = unit(0.5, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(angle = 45, hjust = 0, vjust=0)
        )+
        scale_color_manual(values = cls_tmp) +
        scale_fill_manual(values=cls_tmp) +
        guides(fill = guide_legend(reverse=TRUE)) +
        guides(color = guide_legend(reverse=TRUE))
dev.off()

#===========================
#- VIOLINS of TUB  CELLTYPES
#===========================
#
#- FIGURE 4A

sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors

gens = tolower(c( "cldn6","cldn7" ,"ly6a",
"ttc36","calml4","aldob","slc34a1","mep1a",
"tmem150a","lrp2","tmem174","cdkn1a",
"fxyd2","ccnd1", "tuba1a","mal","tmem213","slc12a1",
"neat1","foxq1","lgals3", "cdkl1"))


df            = data.frame(t(as.matrix(logcounts(sce)[gens,])))
colnames(df)  = str_to_title(colnames(df))
df$cluster    = sce$cluster_tme
dfm           = reshape2::melt(df)
colnames(dfm) = c("cluster","gene","expression")

pw = 22/22*14.5
ph = 9/11*5

pdf(file="../figures/fig4a_np_tubular-violins.pdf",width=pw,height=ph)
cls_tmp = cls[names(cls) %in% df$cluster]
ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
  geom_violin(scale="width") +
  coord_flip() +
  facet_wrap(.~gene,nrow=1) +
  theme_minimal() +
  theme(axis.text.x = element_blank()     ,
        panel.spacing.x = unit(0.5, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(angle = 45, hjust = 0, vjust=0)
        )+
        scale_color_manual(values = cls_tmp) +
        scale_fill_manual(values=cls_tmp) +
        guides(fill = guide_legend(reverse=TRUE)) +
        guides(color = guide_legend(reverse=TRUE))
dev.off()


#========================
#- CROSS LINEAGE MARKDERS
#========================
#
#- FIGURE 5

sce = readRDS("../results/sce_fully-annotated.rds")

#- "known" marker genes for lineages
gens = c("six2","crym","cited1","gdnf","meis1","foxd1","crabp1","aldh1a2")

#- plotting
dat               = data.matrix(logcounts(sce)[gens,])
rownames(dat) = str_to_title(rownames(dat))
colnames(dat)     = sce$cell
coldat            = sce$cluster_it
coldat            = as.data.frame(coldat)
rownames(coldat)  = colnames(dat)
colnames(coldat)  = "cluster"
levels(coldat$cluster) = levels(sce$cluster_it)[levels(sce$cluster_it) %in% levels(coldat$cluster)]

cls               = metadata(sce)$Cluster_colors
cls_tmp           = cls[names(cls) %in% coldat$cluster]
cc                = list(cluster = cls_tmp)
cc$cluster        = cc$cluster[levels(coldat$cluster)]

pal   = colorRampPalette(rev(sequential_hcl(n = 7,"Gray")))(100)
hgaps = cumsum(table(coldat$cluster)[levels(coldat$cluster[order(coldat$cluster)])])
hgaps = hgaps[seq_len(length(hgaps)-1)]


pheatmap(dat[,order(coldat$cluster)],
                scale="none",
                cluster_cols=FALSE,
                gaps_col=hgaps,
                color =pal,
                clustering_method="ward.D",
                annotation_col=coldat,
                annotation_colors=cc,
                cutree_row=8,
                labels_col=rep("",ncol(dat)),
                fontsize_row=9,
                cellheight=10,
                treeheight_row=0)

dev.print(pdf,width=9,height=4.5,file="../figures/fig5_xl_np-stromal_heatmap.pdf")





#- "known" marker genes for lineages
gens = c("six2","crym","cited1","gdnf", "meis1", "foxd1","crabp1","aldh1a2", "dlk1", "postn", "cldn11", "col1a1")

#- plotting
dat               = data.matrix(logcounts(sce)[gens,])
rownames(dat)     = str_to_title(rownames(dat))
colnames(dat)     = sce$cell
coldat            = sce$cluster_it
coldat            = as.data.frame(coldat)
rownames(coldat)  = colnames(dat)
colnames(coldat)  = "cluster"
levels(coldat$cluster) = levels(sce$cluster_it)[levels(sce$cluster_it) %in% levels(coldat$cluster)]

cls               = metadata(sce)$Cluster_colors
cls_tmp           = cls[names(cls) %in% coldat$cluster]
cc                = list(cluster = cls_tmp)
cc$cluster        = cc$cluster[levels(coldat$cluster)]

pal   = colorRampPalette(rev(sequential_hcl(n = 7,"Gray")))(100)
hgaps = cumsum(table(coldat$cluster)[levels(coldat$cluster[order(coldat$cluster)])])
hgaps = hgaps[seq_len(length(hgaps)-1)]


pheatmap(dat[,order(coldat$cluster)],
                scale="none",
                cluster_cols=FALSE,
                gaps_col=hgaps,
                color =pal,
                clustering_method="ward.D",
                annotation_col=coldat,
                annotation_colors=cc,
                cutree_row=12,
                labels_col=rep("",ncol(dat)),
                fontsize_row=9,
                cellheight=10,
                treeheight_row=0)

dev.print(pdf,width=9,height=4.5,file="../figures/fig5_xl_np-stromal_heatmap_Rev.pdf")


#- end
