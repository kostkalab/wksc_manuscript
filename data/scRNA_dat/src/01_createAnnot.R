############################################################
## get gene symbols/ensg id 
############################################################
library(data.table)
library(stringi)
gtf <- fread(cmd = 'grep -v "#" /data/opt/bio/cellRanger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf')
tmp <- data.table(t(stri_list2matrix(lapply(strsplit(gtf$V9, split = ";"), function(x){
    gene_id = gsub("gene_id \\\"|\\\"| ", "", x[grep("gene_id", x)])
    gene_name = gsub("gene_name \\\"|\\\"| ", "", x[grep("gene_name", x)])
    gene_biotype = gsub("gene_biotype \\\"|\\\"| ", "", x[grep("gene_biotype", x)])
    return(c(gene_id, gene_name, gene_biotype))
}))))
setnames(tmp, c("ensembl_gene_id", "gene_name", "gene_biotype"))
tmp1 = cbind(gtf[, c(1,4:5, 7, 3,2)], tmp)
setnames(tmp1, colnames(tmp1)[1:6],  c( "chrom", "start", "end", "strand", "feature", "source"))
tmp2 <- cbind(tmp, tmp1[,1:6])
tmp2 = tmp2[feature == "gene"] ## 27998 x 27998
annot.df <- DataFrame(tmp2)
rownames(annot.df) <- annot.df$ensembl_gene_id
saveRDS(annot.df, file = paste0("../dat/annot.rds"))

############################################################
############################################################
############################################################
