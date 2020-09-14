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
library(openxlsx)

library(limma)
library(org.Mm.eg.db)
library(GO.db)

#===============================================
#- DIFFERENTIALLY EXPRESSED GENES (supp tables )
#===============================================

sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors


#- DEGs: self-renew vs primed
#----------------------------

ii   = sce$cluster_tme %in% c("self-renew","primed")
sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)

mkrs = findMarkers( sce2,
                    sce2$cluster_tme,
                    pval.type="all",
                    direction="any",
                    lfc=log2(1.5),
                    assay.type = "logcounts_SAVER")

ii       = mkrs[[1]]$FDR < 0.1
gg       = rownames(mkrs[[1]])[ii]
gg.sr.pr = gg
gg.sr.pr.up = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed < 0]
gg.sr.pr.dn = gg.sr.pr[mkrs[[1]][ii,]$logFC.primed > 0]


saveRDS(mkrs, file="../results/np_self-vs-primed-DEG.rds")
dat = as.data.frame(mkrs[[1]][ii,])
dat$gene = rownames(dat)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=dat, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/np_self-vs-primed-DEG.xlsx')


#- DEGs: primed vs differentiating
#---------------------------------

ii   = sce$cluster_tme %in% c("primed","differentiating")
sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)

mkrs = findMarkers( sce2,
                    sce2$cluster_tme,
                    pval.type="all",
                    direction="any",
                    lfc=log2(1.5),
                  assay.type = "logcounts_SAVER")

ii       = mkrs[[1]]$FDR < 0.1
gg       = rownames(mkrs[[1]])[ii]
gg.pr.dr = gg
gg.pr.dr.up = gg.pr.dr[mkrs[[1]][ii,]$logFC.differentiating < 0]
gg.pr.dr.dn = gg.pr.dr[mkrs[[1]][ii,]$logFC.differentiating > 0]


saveRDS(mkrs, file="../results/np_primed-vs-differentiating-DEG.rds")
dat = as.data.frame(mkrs[[1]][ii,])
dat$gene = rownames(dat)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=dat, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/np_primed-vs-differentiating-DEG.xlsx')



#==================================================
#- PLOT HEATMAP: 100 random cells from each cluster
#==================================================
#
#- FIGURE 3A

#- Get cells/data and make small sce
#-----------------------------------
set.seed(23432)
gg  = unique(c(gg.sr.pr,gg.pr.dr))
ii  = sort(c(   sample(which(sce$cluster_tme=="self-renew"),100),
                sample(which(sce$cluster_tme=="primed"),100),
                sample(which(sce$cluster_tme=="differentiating"),100)))

sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)

dat  = data.matrix(logcounts(sce2)[gg,])
colnames(dat) = sce2$cell


#- Scale and Saturate
#--------------------
colfun = colorRampPalette((diverge_hsv(n = 7)))(100)
dd     = dat - rowMeans(dat)
dd     = dd / apply(dd,1,sd)
dd[dd>3] = 3 ; dd[dd < -3] = -3

#- choose genes to annotate
#--------------------------
gens = c("birc5","smc4","smc2","six2","cited1","cdca3","lhx1","pax8","crym","meis2","gas1","eya1","rspo1","fxyd2","epcam")
ha = rowAnnotation(foo = anno_mark(at =match(gens,rownames(dat)), labels = gens))

#- Column order in each cluster
#------------------------------
source("./utilities_func.R")
tmp  = cluster_by_fac(t(dat),sce2$cluster_tme,method="ward.D2")
ord  = rev(tmp$ord[tmp$hres$order])


hm = Heatmap(dd[,ord],name="expression",col=colfun,cluster_columns=FALSE,
column_labels=rep("",ncol(dd)), row_labels=rep("",nrow(dd)),
right_annotation = ha, split=3, column_split =sce2$cluster_tme[ord],
clustering_method_rows = "ward.D2",show_row_dend = FALSE,
column_title_gp = gpar(fontsize = 12),row_title_gp = gpar(fontsize = 0))

pdf(width=8,height=6,file="../figures/np_np-clusters-deg-heatmap.pdf")
draw(hm)
dev.off()


#===============
#- GO enrichment
#===============

ensg2eg.lst <- as.list(org.Mm.egENSEMBL2EG)

make_topgo <- function(gens){
  ensg = rowData(sce[gens,])$ID
  entz <- unlist(ensg2eg.lst[ensg])
  go <- goana(entz, species = "Mm")
  topgo <- topGO(go, ontology = "BP", number = 20)
}

tgo.sr.pr = make_topgo(gg.sr.pr)
tgo.pr.dr = make_topgo(gg.pr.dr)


plot_tgo <- function(tgo){
#=========================

  df = tgo
  df <- df[, c(1, 4, 5)]
  colnames(df) = c("TermDesc", "Count", "Pval")
  df$TermDesc <- factor(df$TermDesc, levels = df$TermDesc)

  g <- ggplot(  data = df,
                aes(x = TermDesc, y = -log10(Pval))) +
                geom_bar(stat = "identity",fill=gray(1/2),width=0.75) +
                scale_x_discrete(limits = rev(levels(df$TermDesc))) +
                geom_text(aes(label = df$Count,
                       y = -log10(df$Pval)+1/2),
                       position = position_dodge(0),
                       vjust = 0.5,size=3,color=gray(1/2)) +
                coord_flip() +
                theme_minimal() +
                theme(  plot.title = element_text(hjust = 1),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.y = element_blank(),
                        legend.position = "none",
                        text = element_text(size=14,family = 'Roboto')
                      ) +
                ylab("-log10(p-value)") +
                xlab("")
    return(g)
  }


x11(width=8,height=1.75,type="cairo")
plot_tgo(tgo.pr.dr[1:15,])
dev.print(png,file="../figures/np_np-clusters-go-pr-vs-diff.png", width=8*600,height=3*600,res=600)
plot_tgo(tgo.sr.pr[1:15,])
dev.print(png,file="../figures/np_np-clusters-go-sr-vs-pr.png",width=8*600,height=3*600,res=600)
dev.off(); dev.off()

#- tables

make_topgo <- function(gens){
  ensg = rowData(sce[gens,])$ID
  entz <- unlist(ensg2eg.lst[ensg])
  go <- goana(entz, species = "Mm")
  topgo <- topGO(go, ontology = "BP",number = 1000L)
}

r1 = make_topgo(gg.sr.pr)
r2 = make_topgo(gg.pr.dr)
r1 = r1[r1$N>9,]
r2 = r2[r2$N>9,]
r1 = r1[r1$P.DE<.01,]
r2 = r2[r2$P.DE<.01,]

saveRDS(r1,file="../results/np_np-clusters-go-sr-vs-pr.rds")
saveRDS(r2,file="../results/np_np-clusters-go-pr-vs-diff.rds")

r1$Term.ID = rownames(r1)
r2$Term.ID = rownames(r2

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "self-renew-vs-primed")
addWorksheet(wb = wb, sheetName = "primed-vs-differentiating")
writeData(wb=wb,x=r1, sheet="self-renew-vs-primed",rowNames=FALSE)
writeData(wb=wb,x=r2, sheet="primed-vs-differentiating",rowNames=FALSE)
saveWorkbook(wb,'../results/np_np-clusters-go-combined.xlsx')



#- Pseudotime plots of rspo1, birc5 and ccnd1
#=============================================

sce             <- readRDS("../results/sce_np_fully-annotated.rds")
cls             <- metadata(sce)$Cluster_colors
atme            <- apply(cbind(sce$slingPseudotime_1,sce$slingPseudotime_2,sce$slingPseudotime_3),1,mean,na.rm=T)
sce$slingAvgTme <- atme

df1 = data.frame(expression= assay(sce,"logcounts_SAVER")["rspo1",], cluster = sce$cluster_tme, pseudotime = sce$slingAvgTme)
df1$gene = "rspo1"
df2 = data.frame(expression= assay(sce,"logcounts_SAVER")["birc5",], cluster = sce$cluster_tme, pseudotime = sce$slingAvgTme)
df2$gene = "birc5"
df3 = data.frame(expression= assay(sce,"logcounts_SAVER")["ccnd1",], cluster = sce$cluster_tme, pseudotime = sce$slingAvgTme)
df3$gene = "ccnd1"
df = rbind(df1,df2,df3)


pl = ggplot(df,aes(x=pseudotime,y=expression,col=cluster)) +
    geom_point(size=.7) + scale_color_manual(values = cls) +
    facet_wrap(~gene,ncol=1,scale="free_y") +
  theme_minimal() +
  theme( panel.background = element_blank(),
         panel.border = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
         guides(colour = guide_legend(override.aes = list(size=3)))

         save_plot("../figures/np_pseudotime_birc5-rspo1-ccnd1.pdf",pl,base_heigh = 6, base_asp=.9)
