
library(tidyverse)
library(scater)
library(scran)
library(ggplot2)
library(reshape2)

#- NP CELLS
#==========

#- neat 1 expression
sce_np <- readRDS("../results/sce_np_fully-annotated.rds")

p1_np     <- plotTSNE(sce_np, col="neat1") + theme(legend.position = "none")
tmp       <- plotTSNE(sce_np, col="neat1") 
p1_np_leg <- cowplot::get_legend(tmp) %>% ggpubr::as_ggplot()

#- clusters
cls       <- metadata(sce_np)$Cluster_colors
cls       <- cls[sce_np$cluster_tme %>% as.character() %>% unique()]
p2_np     <- plotTSNE(sce_np, col="cluster_tme")  + theme(legend.position = "none")
p2_np     <- p2_np + scale_color_manual(values=cls) 
tmp       <- plotTSNE(sce_np, col="cluster_tme")
tmp       <- tmp  + scale_color_manual(values=cls)
p2_np_leg <- cowplot::get_legend(tmp) %>% ggpubr::as_ggplot()

#- violins
df            <- data.frame(t(as.matrix(logcounts(sce_np)["neat1",,drop=FALSE])))
df$cluster    <- sce_np$cluster_tme
colnames(df)  <- str_to_title(colnames(df))
dfm           <- reshape2::melt(df)
colnames(dfm) <- c("cluster","gene","expression")

p3_np <- ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
   		geom_violin(scale="width") +
   		facet_wrap(.~gene,nrow=1) +
   		theme_minimal() +
   		scale_color_manual(values=cls) +
   		scale_fill_manual(values=cls) +
   		theme(  legend.position = "none",
           		axis.text.x = element_text(angle = 90, hjust=1))

cowplot::save_plot(p1_np, file="../figures/neat1_p1_np.pdf",base_height=3)
cowplot::save_plot(p1_np_leg, file="../figures/neat1_p1_np_leg.pdf",base_height=3)
cowplot::save_plot(p2_np, file="../figures/neat1_p2_np.pdf",base_height=3)
cowplot::save_plot(p2_np_leg, file="../figures/neat1_p2_np_leg.pdf",base_height=3)
cowplot::save_plot(p3_np, file="../figures/neat1_p3_np.pdf",base_height=3)

#- WHOLE KIDNEY
#==============

#- neat 1 expression
sce_wk <- readRDS("../results/sce_fully-annotated.rds")

#- annotate proximal vs. distal tubular cells
cc = sce_wk$cluster_tme %>% as.character()
cc[ cc %in% c("self-renew","primed")] = 'nephron-progenitor'
cc[ cc %in% c("i_distal_tubular","m_distal_tubular")] = 'distal_tubular'
cc[ cc %in% c("i_proximal_tubular","m_proximal_tubular")] = 'proximal_tubular'
cc[ cc %in% c("i_podocyte","m_podocyte")] = 'podocyte'
sce_wk$cluster = cc

p1_wk     <- plotTSNE(sce_wk, col="neat1") + theme(legend.position = "none")
tmp       <- plotTSNE(sce_wk, col="neat1")
p1_wk_leg <- cowplot::get_legend(tmp) %>% ggpubr::as_ggplot()

#- clusters
cls       <- metadata(sce_wk)$Cluster_colors
cls       <- cls[sce_wk$cluster %>% as.character() %>% unique()]
p2_wk     <- plotTSNE(sce_wk, col="cluster")  + theme(legend.position = "none")
p2_wk     <- p2_wk + scale_color_manual(values=cls)
tmp       <- plotTSNE(sce_wk, col="cluster")
tmp       <- tmp  + scale_color_manual(values=cls)
p2_wk_leg <- cowplot::get_legend(tmp) %>% ggpubr::as_ggplot()

#- violins
df            <- data.frame(t(as.matrix(logcounts(sce_wk)["neat1",,drop=FALSE])))
df$cluster    <- sce_wk$cluster
colnames(df)  <- str_to_title(colnames(df))
dfm           <- reshape2::melt(df)
colnames(dfm) <- c("cluster","gene","expression")

p3_wk <- ggplot(data=dfm, aes(x=cluster, y=expression, col=cluster, fill=cluster)) +
                geom_violin(scale="width") +
                facet_wrap(.~gene,nrow=1) +
                theme_minimal() +
                scale_color_manual(values=cls) +
                scale_fill_manual(values=cls) +
                theme(  legend.position = "none",
                        axis.text.x = element_text(angle = 90, hjust=1))

cowplot::save_plot(p1_wk, file="../figures/neat1_p1_wk.pdf",base_height=3)
cowplot::save_plot(p1_wk_leg, file="../figures/neat1_p1_wk_leg.pdf",base_height=3)
cowplot::save_plot(p2_wk, file="../figures/neat1_p2_wk.pdf",base_height=3)
cowplot::save_plot(p2_wk_leg, file="../figures/neat1_p2_wk_leg.pdf",base_height=3)
cowplot::save_plot(p3_wk, file="../figures/neat1_p3_wk.pdf",base_height=3)
