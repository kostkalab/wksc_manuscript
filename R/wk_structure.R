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
library(openxlsx)
library(ComplexHeatmap)

#- helper functions
source("./utilities_func.R")

#- our SCE
sce = readRDS("../data/scRNA_dat/dat/sce.rds")


#- consistent cluster-color assignment
#=====================================

cls = pal_d3(palette="category20")(15) #- the 15 is post-hoc
names(cls)     = 1:15

#=================================
#- Identify cellt types in kidney
#================================

gens = tolower(c("Cited1","Six2"  , #- NPs
                 "Col1a1","meis1" , #- stroma
                 "Aldh1a2",
                 "Dlk1", "Postn", #- cortical / nephrogenic
                 "Cldn11",
                 "Calb1" ,"Gata3" , #- ureteric bud / collecting duct (UB/CD)
                 "emcn"  ,"kdr"   , #- endothelial
                 "cd52"  ,"fcer1g", #- blood
                 "podxl" ,"nphs1" , #- podocytes
                 "fxyd2" ,"hnf4a" , #- tubular cells (prox, dist. loop of henle)
                 "lhx1"  ,"pax8"))  #- mixed / differentiating

cl = as.character(sce$Cluster)
cl[cl == "1"]  = "ureteric_bud/collecting_duct"
cl[cl == "2"]  = "nephron-progenitor"
cl[cl == "3"]  = "stromal/medullary"
cl[cl == "4"]  = "mixed/differentiating"
cl[cl == "5"]  = "stromal/cortical"
cl[cl == "6"]  = "blood"
cl[cl == "7"]  = "tubular"
cl[cl == "8"]  = "podocytes"
cl[cl == "9"]  = "stromal/mesangial"
cl[cl == "10"] = "endothelial"
sce$Cluster = as.factor(cl)
sce$Cluster = factor(sce$Cluster,levels=c("nephron-progenitor","stromal/cortical","stromal/mesangial","stromal/medullary","ureteric_bud/collecting_duct","endothelial","blood","podocytes","tubular","mixed/differentiating"))


#====================================
#- save the first 10 cluster colors
names(cls)[1:10] = levels(sce$Cluster)
metadata(sce)$Cluster_colors = cls
#=====================================

saveRDS(sce, file="../results/sce_annotated.rds")

################################################################
# That's it for here.                                          #
# Tubular cluster is going to split in the np-related analysis #
################################################################
