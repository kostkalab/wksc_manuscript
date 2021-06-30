#============================
# misc plots 
#============================

#============================
#- Pseudotime plots for TFs; R1.17
#============================
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
library(data.table)
library(stringr)


sce             <- readRDS("../results/sce_np_fully-annotated.rds")
cls             <- metadata(sce)$Cluster_colors
atme            <- apply(cbind(sce$slingPseudotime_1,
                               sce$slingPseudotime_2,
                               sce$slingPseudotime_3),1,
                         mean,na.rm=T)
sce$slingAvgTme <- atme

makeDF <- function(sce, gene){
    df = data.frame(expression = assay(sce, "logcounts_SAVER")[gene, ],
                    cluster = sce$cluster_tme,
                    pseudotime = sce$slingAvgTme)
    df$gene = str_to_title(gene)
    return(df)
}

tfs.lst = list("podocytes" = c("cers6", "foxc2", "lmx1b", "lrrfip1", "mafb", "tcf21", "wt1", "zbtb7c"),
               "distal_tubular" = c("emx1", "foxq1", "gata3", "hoxd8", "mecom", "sim1", "tfap2b"),
               "proximal_tubular" = c("hif3a", "hnf4a", "maf"))

lapply(seq_along(tfs.lst), function(ind){
    x = tfs.lst[[ind]]
    type = names(tfs.lst)[ind]
    x = x[x %in% rownames(sce)]
    
    df.lst = lapply(seq_along(x), function(ii){
        this_g = x[ii]
        df = makeDF(sce, this_g)
        return(df)
    })
    df = rbindlist(df.lst)
    pl = ggplot(df,aes(x=pseudotime,y=expression,col=cluster)) +
        geom_point(size=.7) + scale_color_manual(values = cls) +
        facet_wrap(~gene,ncol=3,scale="free_y") +
        theme_minimal() +
        theme( panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        guides(colour = guide_legend(override.aes = list(size=3)))
    figFile = paste0("../figures/fig-R1.17_np_pseudotime_", type, "-tfs.pdf")
    save_plot(figFile,
              pl,
              base_width = 8,
              base_height = 10)
})



#============================
## birc5 cell cycle phase; R1.18
#============================
sce_wk <- readRDS("../results/sce_fully-annotated_CC.rds")

sce = sce_wk
sce$Birc5 = logcounts(sce)["birc5", ]
p1 <- plotColData(sce,
                  y = "Birc5",
                  x = "seucc_phase",
                  colour_by = "seucc_phase")

figFile = paste0("../figures/fig-R1.18_birc5vsCCPhase_violinplot.pdf")
save_plot(filename = figFile,
          p1,
          base_width = 7,
          base_height = 7)



#============================
## deg primed vs self up/down ; R3.2
#============================
sce = readRDS("../results/sce_np_fully-annotated.rds")
cls = metadata(sce)$Cluster_colors
ii   = sce$cluster_tme %in% c("self-renew","primed")
sce2 = sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)
rownames(sce2) = str_to_title(rownames(sce2))

mwilcox <- function(gene){
    x1 = assay(sce2,"logcounts_SAVER")[gene, sce2$cluster_tme == 'primed']
    x2 = assay(sce2,"logcounts_SAVER")[gene, sce2$cluster_tme != 'primed']
    res = wilcox.test(x1,x2)
}

down_genes = str_to_title(c("six2", "cited1", "rspo1", "crym", "meis2", "gas1"))
wes <- sapply(down_genes, mwilcox)
anno <- data.frame(Feature = names(wes[3,]), lab = paste("p =", format(unlist(wes[3,]),digits=3)))
anno$x = 1.2
anno$y = 0.5
pe <- plotExpression(sce2, features= down_genes, x="cluster_tme", exprs_values="logcounts_SAVER")
p  <- pe +
    geom_text(data = anno, aes(x = x,  y = y, label = lab)) + ylab("Expression") + xlab("NP cell type")
save_plot(p, file="../figures/supFig8_deg_primed-vs-self_down.pdf",base_width=3.5, base_height=7)


up_genes = str_to_title(c("pax8", "lhx1", "wnt4"))
wes <- sapply(up_genes, mwilcox)
anno <- data.frame(Feature = names(wes[3,]), lab = paste("p =", format(unlist(wes[3,]),digits=3)))
anno$x = 1.2
anno$y = 2.1
pe <- plotExpression(sce2, features= up_genes, x="cluster_tme", exprs_values="logcounts_SAVER")
p  <- pe +
    geom_text(data = anno, aes(x = x,  y = y, label = lab)) +
    ylab("Expression") +
    xlab("NP cell type") +
    facet_wrap(vars(Feature), ncol=1)
save_plot(p, file="../figures/supFig8_deg_primed-vs-self_up.pdf",base_width=3.5/2*1.1, base_height=7)


#-end
