library(scater)
library(scran)
library(openxlsx)
library(colorspace)
library(stringr)

#===================================
#- UNEXPECTED GENE EXPRESSION TABLES
#===================================
#
#- Tables 1 and 2

#- TABLE 1
#---------
sce = readRDS("../results/sce_fully-annotated.rds")

fet_fu <- function(g1, g2, mat, subset = NULL, thresh = 0) {
#-------------------------------------------------
#- FET for association between expression of the two genes
#- if subset != NULL only the indicated cells are taken
  if (is.null(subset)) {
    subset <- TRUE
  }
  tmat       <- mat[c(g1,g2),subset,drop=FALSE]
  f1         <- factor(tmat[1,]>thresh,levels=c("TRUE","FALSE"))
  f2         <- factor(tmat[2,]>thresh,levels=c("TRUE","FALSE"))
  tmp        <- fisher.test(table(f1,f2))
  expr.geom  <- mean(sqrt(mat[g1,subset]*mat[g2,subset]))
  res        <- c(tmp$estimate, tmp$p.value, expr.geom)
  names(res) <- c("estimate","pval","expr")
  return(res)
}

gens.1 <- sort(c("col1a1","meis1"))
gens.2 <- sort(c("six2","cited1","crym","gdnf"))
cluss  <- sort(c("cortical_stroma","mesangial_stroma","medullary_stroma"))

clus <- cluss[1]

RES = NULL
for(cl in 1:length(cluss)){
  clus = cluss[cl]
  for(i in 1:length(gens.1)){
    gen1 = gens.1[i]
    #- gene expression in clus:
    num.pos  = sum((sce$cluster_it==clus) & (logcounts(sce)[gen1,]>0))
    frac.pos = round(num.pos / sum(sce$cluster_it==clus),3)

    #- co-expression of gene.2 with gene.1 in clus
    for(j in 1:length(gens.2)){
      message(gen1)
      gen2 = gens.2[j]
      num.pos.g2   = sum((sce$cluster_it==clus) & (logcounts(sce)[gen2,]>0))
      frac.pos.g2 = round(num.pos.g2 / sum(sce$cluster_it==clus),3)
      num.pos.both = sum((sce$cluster_it==clus) & (logcounts(sce)[gen1,]>0) & (logcounts(sce)[gen2,]>0))
      frac.pos.both = round(num.pos.both/sum(sce$cluster_it==clus),3)
      frac.pos.both.exp = frac.pos *  frac.pos.g2
      res  = fet_fu(gen1,gen2,mat=logcounts(sce),subset = sce$cluster_it==clus)
      resl = data.frame(cluster = clus, gene1 = gen1, gene2 = gen2,num_gene1_pos = num.pos, frac_gene1_pos = frac.pos, num_gene2_pos = num.pos.g2, frac_gene2_pos= frac.pos.g2, num_gene1_gene2 = num.pos.both,OR = res[1], p_value = res[2])
      RES = rbind(RES,data.frame(resl))
    }
  }
}
RES = RES[RES$cluster == "cortical_stroma",]
RES$gene1 = str_to_title(RES$gene1)
RES$gene2 = str_to_title(RES$gene2)
saveRDS(RES,"../results/xl_np-stroma_coexpression.rds")
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=RES, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/tab1_xl_np-stroma_coexpression.xlsx')

#- TABLE 2
#---------
gens <- c("calb1","mal","mecom","wfdc2")
bmat <- counts(sce)[gens, ] > 0
r1   <- apply(bmat, 1, function(x) tapply(x, sce$cluster_tme, mean))
r2   <- round(r1,3)*100

tab  <- r2[c("i_proximal_tubular", "m_proximal_tubular", "i_distal_tubular",
             "m_distal_tubular", "ureteric_bud/collecting_duct"),]
colnames(tab) = str_to_title(colnames(tab))
saveRDS(tab,"../results/xl_tub-ub_expression.rds")
tab = cbind(rownames(tab), tab)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=tab, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/tab2_xl_tub-ub_expression.xlsx')

#============================
#- DEGs between WK cell-types
#============================
#
#- Supplemental Table 1

sce = readRDS("../results/sce_fully-annotated.rds")

fm = findMarkers( sce,
                  sce$cluster_it,
                  pval.type="some",
                  direction="any",
                  lfc=log2(1.5))

afu <- function(ind){
  res         = fm[[ind]]
  res         = res[res$FDR <= .05,]
  res$cluster = names(fm)[ind]
  return(res)
}

restab = sapply(1:length(fm),afu)
names(restab) = names(fm)
saveRDS(restab,"../results/wk_celltype-DEG.rds")

wb <- createWorkbook()
names(restab) = sub("/","-",names(restab))
for(i in 1:length(restab)){
  sname = names(restab)[i]
  message(sname)
  df    = as(restab[[i]],"data.frame")
  df$gene = str_to_title(rownames(df))
  addWorksheet(wb = wb, sheetName = sname)
  writeData(wb=wb,sheet=sname,x=df,rowNames=FALSE)
  conditionalFormatting(wb, sname, cols=3:(ncol(df)-2), rows=1:(nrow(df)+1),
          style = diverging_hcl(5, palette = "Blue-Red 3")[2:4],
          rule = c(-2,0,2),
          type = "colourScale")
}
saveWorkbook(wb,'../results/supTab1_wk_celltype-DEG.xlsx')

#============================
#- DEGs between NP cell-types
#============================
#
#- Supplemental Table 2

sce = readRDS("../results/sce_np_fully-annotated.rds")

fm = findMarkers( sce,
                  sce$cluster_tme,
                  pval.type="some",
                  direction="any",
                  lfc=log2(1.5))

afu <- function(ind){
  res         = fm[[ind]]
  res         = res[res$FDR <= .05,]
  res$cluster = names(fm)[ind]
  return(res)
}

restab = sapply(1:length(fm),afu)
names(restab) = names(fm)
saveRDS(restab,"../results/np_celltype-DEG.rds")

wb <- createWorkbook()
names(restab) = sub("/","-",names(restab))
for(i in 1:length(restab)){
  sname = names(restab)[i]
  message(sname)
  df    = as(restab[[i]],"data.frame")
  df$gene = str_to_title(rownames(df))
  addWorksheet(wb = wb, sheetName = sname)
  writeData(wb=wb,sheet=sname,x=df,rowNames=FALSE)
  conditionalFormatting(wb, sname, cols=3:(ncol(df)-2), rows=1:(nrow(df)+1),
          style = diverging_hcl(5, palette = "Blue-Red 3")[2:4],
          rule = c(-2,0,2),
          type = "colourScale")
}
saveWorkbook(wb,'../results/supTab2_np_celltype-DEG.xlsx')

#=================================
#- MATURE vs. IMMATURE comparisons
#=================================
#
#- Supplemental Table

sce = readRDS("../results/sce_np_fully-annotated.rds")

fm = findMarkers( sce,
                  sce$cluster_tme,
                  pval.type = "any",
                  direction = "any",
                  lfc = log2(1.5),
                  full.stats = TRUE)

DAT           <- NULL
ind           <- which(fm[["m_podocyte"]][,"stats.i_podocyte"][,"log.FDR"] < log(0.1))
dat           <- fm[["m_podocyte"]][,"stats.i_podocyte"][ind,]
dat$cluster_1 <- "m_podocyte"
dat$cluster_2 <- "i_podocyte"
DAT           <- rbind(DAT,dat)

ind           <- which(fm[["m_distal_tubular"]][,"stats.i_distal_tubular"][,"log.FDR"] < log(0.1))
dat           <- fm[["m_distal_tubular"]][,"stats.i_distal_tubular"][ind,]
dat$cluster_1 <- "m_distal_tubular"
dat$cluster_2 <- "i_distal_tubular"
DAT           <- rbind(DAT,dat)

ind           <- which(fm[["m_proximal_tubular"]][,"stats.i_proximal_tubular"][,"log.FDR"] < log(0.1))
dat           <- fm[["m_proximal_tubular"]][,"stats.i_proximal_tubular"][ind,]
dat$cluster_1 <- "m_proximal_tubular"
dat$cluster_2 <- "i_proximal_tubular"
DAT           <- rbind(DAT,dat)

rownames(DAT) = str_to_title(rownames(DAT))
saveRDS(DAT,"../results/np_celltype-mature-vs-immature-DEG.rds")

DAT$gene = rownames(DAT)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=DAT, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/supTab9_np_celltype-mature-vs-immatrue-DEG.xlsx')


DAT           <- NULL
ind           <- which(fm[["m_proximal_tubular"]][,"stats.m_distal_tubular"][,"log.FDR"] < log(0.1))
dat           <- fm[["m_proximal_tubular"]][,"stats.m_distal_tubular"][ind,]
dat$cluster_1 <- "m_proximal_tubular"
dat$cluster_2 <- "m_distal_tubular"
DAT           <- rbind(DAT,dat)
rownames(DAT) = str_to_title(rownames(DAT))


DAT$gene = rownames(DAT)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=DAT, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/supTab10_np_celltype-tub-dist-vs-prox-DEG.xlsx')



#- end
