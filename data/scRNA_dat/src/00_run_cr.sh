#!/bin/bash

CELLRANGER=/data/opt/bio/cellRanger/cellranger-2.2.0/cellranger

nice -19 $CELLRANGER count 	--id=Single1 								\
				--transcriptome=/data/opt/bio/cellRanger/refdata-cellranger-mm10-1.2.0 	\
				--fastqs=/projects/wksc/data/raw 					\
				--sample=Single1 							\
				--force-cells=10000 							\
				--chemistry=SC3Pv2 							\
				--localcores=12								\


