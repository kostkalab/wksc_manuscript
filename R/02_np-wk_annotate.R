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
library(SAVER)

#- Annotate NP celltypes
sce = readRDS("../results/sce_np.rds")
#=====================================
#- CLUSTER THE NP-LINEAGE CELLS ALONE
#=====================================

snn.gr       = buildSNNGraph(sce,k=25, use.dimred="PCA")
clusters     = igraph::cluster_walktrap(snn.gr, steps = 10)
sce$Cluster2 = factor(igraph::cut_at(clusters, 6))

plotTSNE(sce, colour_by="Cluster2")

cl = as.character(sce$Cluster2)
cl[cl == "1"] = "differentiating"
cl[cl == "2"] = "podocytes"
cl[cl == "3"] = "primed"
cl[cl == "4"] = "tubular_prox"
cl[cl == "5"] = "tubular_dist"
cl[cl == "6"] = "self-renew"
sce$Cluster2 = factor(cl, levels = c("self-renew",
                                     "primed",
                                     "differentiating",
                                     "tubular_prox",
                                     "tubular_dist",
                                     "podocytes"))

#- Cluster2 colors are still needed
cl   = metadata(sce)$Cluster_colors
lv   = levels(sce$Cluster2)
new  = lv[! lv %in% names(cl)]
nold = length(levels(sce$Cluster))
names(cl)[(nold+1):(nold+length(new))] = new

metadata(sce)$Cluster_colors  = cl
saveRDS(sce,file="../results/sce_np_annotated.rds")

#===============================================
#- UPDATE THE WK CLUSTERS WITH TUB DIST vs. PROX
#===============================================

sce_np    = readRDS("../results/sce_np_annotated.rds")
sce       = readRDS("../results/sce_annotated.rds")
cl        = sce$Cluster
cl        = as.character(cl)
names(cl) = sce$cell

colnames(sce_np) = sce_np$cell
ids = sce$cell[sce$Cluster %in% c("mixed/differentiating","tubular")]
cl[ids] = as.character(sce_np[,ids]$Cluster2)
cl[cl %in% c("primed","self-renew")] = "nephron-progenitor"
cl[cl %in% c("differentiating")] = "mixed/differentiating"
cl[cl %in% c("tubular_prox")] = "tubular/proximal"
cl[cl %in% c("tubular_dist")] = "tubular/distal"
cl = factor(cl,levels=c("nephron-progenitor","stromal/cortical","stromal/mesangial","stromal/medullary","ureteric_bud/collecting_duct","endothelial","blood","podocytes","tubular/distal","tubular/proximal","mixed/differentiating"))
sce$Cluster_fine = cl

metadata(sce)$Cluster_colors["tubular/distal"]   = metadata(sce_np)$Cluster_colors["tubular_dist"]
metadata(sce)$Cluster_colors["tubular/proximal"] = metadata(sce_np)$Cluster_colors["tubular_prox"]

saveRDS(sce,file="../results/sce_annotated_fine.rds")

#=========================================
#- PSEUDOTIME ANNOTAION: MATURE / IMMATURE
#=========================================

sce = readRDS("../results/sce_np_annotated.rds")
sce = slingshot(  sce,
                  clusterLabels = 'Cluster2',
                  reducedDim = 'TSNE',
                  start.clus="self-renew")

#- SLINGSHOT LINEAGE ANNOTATIONS
################################

sds = SlingshotDataSet(sce)

#- cells that are going to be podocytes
t1      <- slingPseudotime(sds, na=FALSE)[,1]
t2      <- slingPseudotime(sds, na=FALSE)[,2]
t3      <- slingPseudotime(sds, na=FALSE)[,3]

w1      <- slingCurveWeights(sds)[,1]
w2      <- slingCurveWeights(sds)[,2]
w3      <- slingCurveWeights(sds)[,3]

ind_2 <- (w2 >= .8) & (w1 < .2) & (w3 < .2)
ind_1 <- (w1 >= .8) & (w2 < .2) & (w3 < .2)
ind_3 <- (w3 >= .8) & (w1 < .2) & (w2 < .2)

ind_12 <- (w3 < .2) & (!ind_1) & (!ind_2)
ind_13 <- (w2 < .2) & (!ind_1) & (!ind_3)
ind_23 <- (w1 < .2) & (!ind_2) & (!ind_3)

ind_123 <- (!ind_1) & (!ind_2) & (!ind_3) & (!ind_12) & (!ind_23) & (!ind_13)

#- TIME SMOOTHING / function
#- predict cluster label by pseudotime
mke_tme_clus_prob <- function(tme,clua=sce$Cluster4,nms) {
#======================================================

  clua = clua[!is.na(tme)]
  clua = droplevels(clua)
  levs = levels(clua)[sort(unique(as.numeric(clua)))]
  nms  =  nms[!is.na(tme)]
  tme  = tme[!is.na(tme)]
  ncla = length(unique(clua))
  tdat = data.frame(t=tme,y=as.integer(clua)-1)

  tres = nnet::multinom(y~t,tdat)
  prd  = tres$fitted.values

  colnames(prd) = levs
  dprd = data.frame(prd)
  dprd$tme = tme
  dprd$cluster = clua
  rownames(dprd) = nms

  return(list(res=tres,prd=dprd, tme=tme))
}

#- LINEAGE-COMMITTED PODOCYTES
#==============================

table(ind_1,sce$Cluster2=="podocytes")
cl                                      = as.character(sce$Cluster2)
cl[cl=="podocytes"]                     = "m_podocytes"
cl[ind_1 & (sce$Cluster2!="podocytes")] = "i_podocytes"
cl                                      = as.factor(cl)

#- smooth via time
#-----------------
#- minimum podocyte time for i_podocytes (and m_podocytes)
tmp      = mke_tme_clus_prob(t1,cl,nms=colData(sce)$cell)
tmin.ind = min(which(colnames(tmp$prd)[1:7][apply(tmp$prd[order(t1),][,1:7],1,which.max)] == "i_podocytes"))
tmin     = sort(t1)[tmin.ind]
#- for t>=tmin now use the time-based predictions; don't touch t < tmin
i1       = (t1>= tmin) & (ind_1)
p1       = c("i_podocytes","m_podocytes")[apply(tmp$prd[i1,c("i_podocytes","m_podocytes")],1,which.max)]
cl       = as.character(cl)
cl[i1]   = p1
#- for t < tmin, replace podocyte classification with next-best:
i2       = (!i1) & (cl %in% c("i_podocytes","m_podocytes")) & (!ind_2) & (!ind_3)
p2       = c("primed","self-renew","differentiating")[apply(tmp$prd[i2,c("primed","self.renew","differentiating")],1,which.max)]
cl[i2]   = p2
sce$Cluster3 = as.factor(cl)
TMIN_POD = tmin

#- LINEAGE-COMMITTED DISTAL TUBULAR CELLS
#=========================================
table(ind_2,sce$Cluster3=="tubular_dist")
cl                                      = as.character(sce$Cluster3)
cl[cl=="tubular_dist"]                  = "m_tubular_dist"
cl[ind_2 & (sce$Cluster2!="tubular_dist")] = "i_tubular_dist"
cl                                      = as.factor(cl)

#- smooth via time
#-----------------
#- minimum time
tmp      = mke_tme_clus_prob(t2,cl,nms=colData(sce)$cell)
tmin.ind = min(which(colnames(tmp$prd)[1:7][apply(tmp$prd[order(t2),][,1:7],1,which.max)] == "i_tubular_dist"))
tmin     = sort(t2)[tmin.ind]
#- for t>=tmin now use the time-based predictions; don't touch t < tmin
i1       = (t2>= tmin) & (ind_2)
p1       = c("i_tubular_dist","m_tubular_dist")[apply(tmp$prd[i1,c("i_tubular_dist","m_tubular_dist")],1,which.max)]
cl       = as.character(cl)
cl[i1]   = p1
#- for t < tmin, replace podocyte classification with next-best:
i2       = (!i1) & (cl %in% c("i_tubular_dist","m_tubular_dist")) & (!ind_1) & (!ind_3)
p2       = c("primed","self-renew","differentiating")[apply(tmp$prd[i2,c("primed","self.renew","differentiating")],1,which.max)]
cl[i2]   = p2

sce$Cluster4 = as.factor(cl)
TMIN_DT = tmin

#- LINEAGE-COMMITTED PROXIMAL TUBULAR CELLS
#=========================================
#- here we do a bit of re-clustering, because otherwise almost no immature cells
snn.gr       =  buildSNNGraph(sce, use.dimred="PCA")
clusters     =  igraph::cluster_walktrap(snn.gr, steps = 10)
sce$ClusterX = factor(igraph::cut_at(clusters, 7))
cl = as.character(sce$Cluster4)
cl[cl=="tubular_prox"] = "m_tubular_prox"
cl[(sce$ClusterX == "4")] = "i_tubular_prox"
cl = as.factor(cl)
#- time-smoothing
tmp = mke_tme_clus_prob(t3,cl,nms=colData(sce)$cell)
tmin.ind = min(which(colnames(tmp$prd)[1:7][apply(tmp$prd[order(t3),][,1:7],1,which.max)] == "i_tubular_prox"))
tmin     = sort(t3)[tmin.ind]
i1       = (t3>= tmin) & (ind_3)
p1       = c("i_tubular_prox","m_tubular_prox")[apply(tmp$prd[i1,c("i_tubular_prox","m_tubular_prox")],1,which.max)]
cl       = as.character(cl)
cl[i1]   = p1
i2       = (!i1) & (cl %in% c("i_tubular_prox","m_tubular_prox")) & (!ind_2) & (!ind_1)
p2       = c("primed","self-renew","differentiating")[apply(tmp$prd[i2,c("primed","self.renew","differentiating")],1,which.max)]
cl[i2]   = p2

sce$Cluster5 = as.factor(cl)
TMIN_PT = tmin


########################
#- UPDATE CLUSTER COLORS
########################

flatten_alpha <- function(x,al,bg=rep(255,3)) {
  cl = col2rgb(x)
  cl = (1-al) * bg + al*cl
  rgb(cl[1],cl[2],cl[3],maxColorValue = 255)
}

md  = metadata(sce)
cls = md$Cluster_colors
names(cls) = str_replace(names(cls),"tubular_dist","m_tubular_dist")
names(cls) = str_replace(names(cls),"tubular_prox","m_tubular_prox")
names(cls) = str_replace(names(cls),"podocytes","m_podocytes")
cls = c(cls, flatten_alpha(alpha(cls["m_tubular_dist"],.5),.5)); names(cls)[length(cls)] = "i_tubular_dist"
cls = c(cls, flatten_alpha(alpha(cls["m_tubular_prox"],.5),.5)); names(cls)[length(cls)] = "i_tubular_prox"
cls = c(cls, flatten_alpha(alpha(cls["m_podocytes"]   ,.5),.5)); names(cls)[length(cls)] = "i_podocytes"

metadata(sce)$Cluster_colors = cls
saveRDS(sce, file="../results/sce_np_annotated_time-clustered.rds")

#=============
#- RUN SAVER
#=============

if( !file.exists("../results/sce_np_annotated_time-clustered_saver.rds") ){
   sce <- readRDS("../results/sce_np_annotated_time-clustered.rds")
    sv = saver(counts(sce),ncores=1)
    assay(sce,"logcounts_SAVER") = log2(1+sv$estimate)
    metadata(sce)$result_SAVER   = sv
    saveRDS(sce, file="../results/sce_np_annotated_time-clustered_saver.rds")
} else {
  sce = readRDS("../results/sce_np_annotated_time-clustered_saver.rds")
}

#===================================
#- CLEAN UP FILE AND ANNOTAION CHAOS
#===================================

sce_wk           <- readRDS("../results/sce_annotated_fine.rds")
sce_np           <- readRDS("../results/sce_np_annotated_time-clustered_saver.rds")
colnames(sce_wk) <- sce_wk$cell
colnames(sce_np) <- sce_np$cell

sce_np$cluster     <- as.character(sce_np$Cluster)
sce_np$cluster_it  <- as.character(sce_np$Cluster2)
sce_np$cluster_tme <- as.character(sce_np$Cluster5)
colData(sce_np)    <- colData(sce_np)[, -grep("Cluster*", colnames(colData(sce_np)))]

sce_wk$cluster     <- as.character(sce_wk$Cluster)
sce_wk$cluster_it  <- as.character(sce_wk$Cluster_fine)
sce_wk$cluster_tme <- as.character(sce_wk$Cluster_fine)

ids                       <- colnames(sce_np)
sce_wk[, ids]$cluster_it  <- sce_np$cluster_it
sce_wk[, ids]$cluster_tme <- sce_np$cluster_tme

sce_wk$cluster_it[sce_wk$cluster_it == "self-renew"]      <- "nephron-progenitor"
sce_wk$cluster_it[sce_wk$cluster_it == "primed"]          <- "nephron-progenitor"
sce_wk$cluster_it[sce_wk$cluster_it == "differentiating"] <- "mixed/differentiating"

#- make factors with "correctly" ordered levels
#----------------------------------------------

mk_new_labs <- function(labs){
  labs <- gsub("blood",             "hematopoietic",    labs)
  labs <- gsub("stromal/medullary", "medullary_stroma", labs)
  labs <- gsub("stromal/cortical",  "cortical_stroma",  labs)
  labs <- gsub("stromal/mesangial", "mesangial_stroma", labs)
  labs <- gsub("tubular_prox",      "proximal_tubular", labs)
  labs <- gsub("tubular/proximal",  "proximal_tubular", labs)
  labs <- gsub("tubular_dist",      "distal_tubular",   labs)
  labs <- gsub("tubular/distal",    "distal_tubular",   labs)
  labs <- gsub("podocytes",         "podocyte",         labs)
  labs <- gsub("i_podocytes",       "i_podocyte",       labs)
  labs <- gsub("m_podocytes",       "m_podocyte",   labs)
  return(labs)
}

all_ordered = c("self-renew", "primed", "nephron-progenitor",
"differentiating", "mixed/differentiating", "i_podocyte", "m_podocyte",
"podocyte", "i_proximal_tubular", "m_proximal_tubular", "proximal_tubular",
"i_distal_tubular", "m_distal_tubular", "distal_tubular", "tubular",
"ureteric_bud/collecting_duct", "cortical_stroma", "mesangial_stroma",
"medullary_stroma", "endothelial", "hematopoietic")

mk_fac <- function(labs){
  tmp <- mk_new_labs(labs)
  tmp <- factor(tmp, levels = all_ordered[all_ordered %in% tmp])
  return(tmp)
}

sce_wk$cluster     <- mk_fac(sce_wk$cluster)
sce_wk$cluster_it  <- mk_fac(sce_wk$cluster_it)
sce_wk$cluster_tme <- mk_fac(sce_wk$cluster_tme)

sce_np$cluster     <- mk_fac(sce_np$cluster)
sce_np$cluster_it  <- mk_fac(sce_np$cluster_it)
sce_np$cluster_tme <- mk_fac(sce_np$cluster_tme)

cls_wk <- metadata(sce_wk)$Cluster_colors
cls_np <- metadata(sce_np)$Cluster_colors

cnms          <- names(cls_wk)
names(cls_wk) <- mk_new_labs(cnms)
cls_wk        <- cls_wk[names(cls_wk) %in% all_ordered]
cnms          <- names(cls_np)
names(cls_np) <- mk_new_labs(cnms)
cls_np        <- cls_np[names(cls_np) %in% all_ordered]
cls           <- c(cls_wk , cls_np[!(names(cls_np) %in% names(cls_wk))])

metadata(sce_wk)$Cluster_colors = cls
metadata(sce_np)$Cluster_colors = cls

saveRDS(sce_wk, "../results/sce_fully-annotated.rds")
saveRDS(sce_np, "../results/sce_np_fully-annotated.rds")


#- END
