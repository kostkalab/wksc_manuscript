library(SingleCellExperiment)

sce = readRDS("../results/sce_fully-annotated.rds")
cts = counts(sce)
lnc = logcounts(sce)
md  = data.frame(barcode=sce$cell,cluster=as.character(sce$cluster_it),cluster_tme=as.character(sce$cluster_tme))

Matrix::writeMM(cts,"../results/GEO_counts.mtx")
Matrix::writeMM(lnc,"../results/GEO_log-norm-counts.mtx")
readr::write_csv(md,"../results/GEO_metadata.csv")
readr::write_csv(data.frame(genes=rownames(sce)),"../results/GEO_genes.txt")
readr::write_csv(data.frame(barcodes=sce$cell),  "../results/GEO_barcodes.txt")

#-  make the archive:
#-  cd ../results
#-  tar -czvf ./GEO_processed-data.tgz ./GEO_counts.mtx ./GEO_log-norm-counts.mtx ./GEO_genes.txt ./GEO_barcodes.txt ./GEO_metadata.csv
#-  ln -s ./GEO_processed-data.tgz Single1_processed.tgz
