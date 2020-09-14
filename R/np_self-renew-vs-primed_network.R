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
library("SCENIC")
library(mgcv)
library(stringr)
library(foreach)

#=========================
#- SCENIC regulon analysis
#=========================

#- set up SCENIC / bend paths
dburl    <- "https://resources-mirror.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
org      <- "mgi" # or hgnc, or dmel
dbDir    <- "../data/external/rcistargetdat/" # RcisTarget databases location
filename <- paste(dbDir,"/",basename(dburl),sep="")
if(!file.exists(filename)) download.file(dburl, destfile <- filename)

myDatasetTitle="SCENIC_analysis" # choose a name for your analysis
data(defaultDbNames)
dbs <- "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"

scenicOptions <- initializeScenic(  org=org, dbDir=dbDir, dbs=dbs,
                                    datasetTitle=myDatasetTitle, nCores=1)

scenicOptions@fileNames$int = gsub("int/","../results/tmp/SCENIC/int/",scenicOptions@fileNames$int)

scenicOptions@fileNames$output = gsub("output/","../results/tmp/SCENIC/output/",scenicOptions@fileNames$output)
scenicOptions@settings$seed <- 123


#- FIND GENES TO SELECT -> PSEUDOTIME ASSOCIATED GENES
#=====================================================

#- 1. Shorten sce to have only self-renew, primed and differentiating cells
#--------------------------------------------------------------------------

sce  <- readRDS("../results/sce_np_fully-annotated.rds")
ii   <- sce$cluster_tme %in% c("self-renew","primed","differentiating")
sce2 <- sce[,ii] ; sce2$cluster_tme = droplevels(sce2$cluster_tme)


#- 2. Use time-smoothing to collect self-renew and primed cells mainly
#-------------------------------------------------------------------

mke_tme_clus_prob <- function(tme,clua,nms) {
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

atme = apply(cbind(sce2$slingPseudotime_1,sce2$slingPseudotime_2,sce2$slingPseudotime_3),1,mean,na.rm=T)
tmeA = mke_tme_clus_prob(atme,sce2$cluster_tme,nms=colData(sce2)$cell)
sce2$slingAvgTme = atme

#- find a cutoff time for "early" cells

rfu <- function(tme){
  predict(tmeA$res,newdat = data.frame(t=tme),type="prob")[2]
}

#- misses 20 primed cells and includes 28 differentiating - pretty good.
mxPr          = max(atme[sce2$cluster_tme=="primed"])
mnPr          = atme[which.max(sapply(atme,rfu))]
tcut          = uniroot(function(x) rfu(x)-.5,lower=mnPr,upper=mxPr)
tcut          = tcut$root
sce3          = sce2[,sce2$slingAvgTme < tcut]
sce3$cluster_tme = droplevels(sce3$cluster_tme)
sce3          = sce3[rowSums(logcounts(sce2)>0) > floor(ncol(sce3)*.05),]

#- 3. Calcluate pseudotime-associated genes (across selected cells)
#------------------------------------------------------------------

#- associated genes
#==================
gam_gene_assoc <- function(gname,tsce){
#=======================================
  message(gname)

  #ex  = logcounts(tsce)[gname,]
  ex  = assay(tsce,"logcounts_SAVER")[gname,]
  tm  = tsce$slingAvgTme
  sf  = colMeans(assay(tsce,"logcounts_SAVER"))
  dat = data.frame(tm=tm,ex=ex,sf=sf)
  res = gam(ex ~ sf + s(tm))
  pvl = summary(res)$s.table[1,"p-value"]
  #pfu = function(tm,sf)  predict( res,
  #                    newdata= data.frame(tm=tm,sf=sf))
  return(list(res=res,pval=pvl))#,predfu=pfu))
}

#
mgv = modelGeneVar(sce3,assay.type="logcounts_SAVER")
ind = mgv$FDR < 0.1

gam.res = lapply(rownames(sce3)[ind],gam_gene_assoc,sce3)
names(gam.res) = rownames(sce3)[ind]
gam.res = gam.res[unlist(lapply(gam.res,function(x) ((x$res$mgcv.conv$fully.converged))))]
pvs = unlist(lapply(gam.res,function(x) x$pval))

mm = assay(sce3,"logcounts_SAVER")[names(gam.res),]
iii = (pvs*length(pvs) < 0.01)

crs    = apply(mm,1,cor,sce3$slingAvgTme,method="spearman")
mds = apply(mm,1,mad)
up.ind = (crs > .4)  & iii & (mds > median(mds))
dn.ind = (crs < -.4)  & iii & (mds > median(mds))
tgens  = names(pvs[iii])
up.gens = names(pvs[up.ind])
dn.gens = names(pvs[dn.ind])

dat1 = data.frame(gene=dn.gens,direction="down")
dat2 = data.frame(gene=up.gens,direction="up")
dat  = rbind(dat1,dat2)

saveRDS(dat,file="../results/np_self-renew-vs-primed-pseudotime_upDown.rds")

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "increasing_with_time")
addWorksheet(wb = wb, sheetName = "decreasing_with_time")
writeData(wb=wb,x=dat1, sheet="decreasing_with_time",rowNames=FALSE)
writeData(wb=wb,x=dat2, sheet="increasing_with_time",rowNames=FALSE)
saveWorkbook(wb,'../results/np_self-renew-primed-pseudotime-upDown.xlsx')



#- UP  AND DOWN GO ANALYSIS


library(GSEABase)
gsc <- getBroadSets("../data/external/MSigDB/msigdb_v7.0.xml")
nuniv = length(unique(unlist(geneIds(gsc))))

fetfu <- function(set1,set2,nuniv){
#----------------------------------
  n1 = length(set1)
  n2 = length(set2)
  ov = sum(set1 %in% set2)
  tb = rbind( c( ov,    n1-ov         ),
              c( n2-ov, nuniv-n1-n2+ov))
  return(fisher.test(tb))
}

fets_dn <- sapply(gsc,
                  function(x) fetfu(dn.gens, str_to_lower(geneIds(x)),nuniv ))

fets_up <- sapply(gsc,
                  function(x) fetfu(up.gens, str_to_lower(geneIds(x)),nuniv ))

mk_tab <- function(fets,gens){
  pvs  = unlist(fets["p.value",])
  pvsa = p.adjust(pvs)
  ord  = order(pvs,decreasing=FALSE)
  nms  = sapply(gsc[ord][pvsa[ord]<.01],setName)
  cat  = sapply(gsc[ord][pvsa[ord]<.01],function(x)
                bcCategory(collectionType(x)))
  scat = sapply(gsc[ord][pvsa[ord]<.01],function(x)
                bcSubCategory(collectionType(x)))
  gns  = sapply(gsc[ord][pvsa[ord]<.01],function(x)
                paste(sort(intersect(str_to_lower(geneIds(x)),gens)),collapse=";"))
  desc = sapply(gsc[ord][pvsa[ord]<.01],function(x) description(x))
  pp   = format(pvs[ord][pvsa[ord]<.01],digits=4,scitentific=TRUE)
  ppa  = format(pvsa[ord][pvsa[ord]<.01],digits=4,scitentific=TRUE)

  tbl = cbind(nms,cat,scat,pp,ppa,gns)
  colnames(tbl) =c("term_name","term_category","term_subcategory","p_value","adjusted_p_value","genes")

  return(tbl)
}

tbl = mk_tab(fets_dn,dn.gens)
i1 = tbl[,"term_category"] == "h"
i2  = tbl[,"term_subcategory"]  == "BP"
i2[is.na(i2)] = FALSE
res.dn = tbl[i1|i2,]

tbl = mk_tab(fets_up,up.gens)
i1 = tbl[,"term_category"] == "h"
i2  = tbl[,"term_subcategory"]  == "BP"
i2[is.na(i2)] = FALSE
res.up = tbl[i1|i2,]

saveRDS(res.dn,file="../results/np_self-renew-primed-pseudotime_go-dn.rds")
saveRDS(res.up,file="../results/np_self-renew-primed-pseudotime_go-up.rds")

res.up = data.frame(res.up)
res.dn = data.frame(res.dn)

res.dn$p_value          = as.numeric(res.dn$p_value)
res.dn$adjusted_p_value = as.numeric(res.dn$adjusted_p_value)
res.up$p_value          = as.numeric(res.up$p_value)
res.up$adjusted_p_value = as.numeric(res.up$adjusted_p_value)

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "increasing_with_time")
addWorksheet(wb = wb, sheetName = "decreasing_with_time")
writeData(wb=wb,x=res.dn, sheet="decreasing_with_time",rowNames=FALSE)
writeData(wb=wb,x=res.up, sheet="increasing_with_time",rowNames=FALSE)
saveWorkbook(wb,'../results/np_self-renew-primed-pseudotime-go.xlsx')

#- PERFORM SCENIC ANALYSIS
#=========================

exprMat = mm[sort(tgens),]
rownames(exprMat) = str_to_title(rownames(exprMat))
colnames(exprMat) = colnames(sce3)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)


aucRegulons = readRDS("../results/tmp/SCENIC/int/3.1_regulons_forAUCell.Rds")
aucRes = readRDS("../results/tmp/SCENIC/int/3.4_regulonAUC.Rds")

#- use only extended wherer possible....; create list with modules.
mat = assay(aucRes)

TF     <-  sort(unique(gsub("_.*", "", gsub(" ", "_", rownames(mat)))))
TF.ind <-  sapply(TF, grep, rownames(mat))
afu <- function(ids){
  if(length(ids)==1) return(ids)
  ii = grep("_extended",rownames(mat)[ids])
  return(ids[ii])
}
inds = unlist(lapply(TF.ind,afu))
 df = melt(aucRegulons[rownames(mat)][inds])
colnames(df) <- c("gene","module")

saveRDS(df,file="../results/np_pseudotime-gene-set-modules.rds")
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Sheet1")
writeData(wb=wb,x=df, sheet="Sheet1",rowNames=FALSE)
saveWorkbook(wb,'../results/np_pseudotime-gene-set-modules.xlsx')



ca = HeatmapAnnotation( cluster=as.character(sce3$cluster_tme),
                        pseudotime=sce3$slingAvgTme,
                        col = list(cluster=metadata(sce3)$Cluster_colors,
                                   pseudotime= circlize::colorRamp2(c(0, 10, 20), c("white", "gray", "black"))) )
colfun = colorRampPalette((diverge_hsv(n = 7)))(100)
pmat = mat[inds,]
pmat = t(scale(t(pmat)))
pmat[pmat>2] = 2
pmat[pmat< -2] = -2
rownames(pmat) = gsub("_.*", "", gsub(" ", "_", rownames(pmat)))
hm =  Heatmap(t(scale(t(pmat)))[,order(sce3$slingAvgTme)],cluster_columns=FALSE,
column_labels=rep("",ncol(mat)),show_row_dend = FALSE,split=2, col=colfun,
clustering_method_rows = "ward.D2",top_annotation=ca[order(sce3$slingAvgTme)],
name="module score")

pdf("../figures/np_pseudotime-gene-set-heatmap.pdf",height=6,width=14)
draw(hm)
dev.off()
