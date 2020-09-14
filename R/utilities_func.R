
cluster_by_fac <- function(dat,fac,...){
#=======================================

  ord = order(fac)
  dat = dat[ord,]

  lfu = function(lev,fac,dat,...){
  #---------------------------
    ii  = fac == lev
    res = hclust(dist(dat[ii,]),...)
  }

  trs = lapply(levels(fac),lfu,fac,dat,...)
  trs = lapply(trs, as.dendrogram)
  res = trs[[1]]
  for(i in 2:length(trs)){
    res = merge(res,trs[[i]])
  }
  hres = as.hclust(res)
  return(list(hres=hres,ord=ord))

}
