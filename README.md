

###Single cell RNA sequencing reveals differential cell cycle activity in key cell populations during nephrogenesis

<p align=center>
Abha S. Bais, Débora M. Cerqueira, Andrew  Clugston, Jacqueline Ho and Dennis Kostka
</p>

This repository contains code and instruction to reproduce analyses carried out in the context of the manuscript ***Single cell RNA sequencing reveals differential cell cycle activity in key cell populations during nephrogenesis***. Pleas follow [this link](TBD) to access the manuscript on bioRxiv. Data generated is available on the Gene Expression Omnibus (GEO) under [TBD](TBD).

The overall directory structure is the following:
```
.
├── data
├── figures
├── R
├── readme.md (this file)
└── results
```

* ```data``` contains new and existing/annotation data used for the project.
* ```figures``` contains plots/figures output from ```R``` scripts
* ```results``` contains other results like gene lists and result-objects derived from input data


The other directories are self-explanatory, ```results``` contains all results we report except the figures, which are in the ```figures``` directory.


**Single cell experssion data**

Single cell expression data is located in the data/scRNA_dat subdirectory:

```
./data/scRNA_dat/
├── cr_out -> [CELLRANGER_OUTPUT_DIR]
├── dat
│   ├── annot.rds
│   ├── sce.rds
│   └── sce.rds_md5sum
├── raw -> [RAW_FASTQ_DIR]
└── src
    ├── 00_run_cr.sh
    ├── 01_createAnnot.R
    └── 02_processRawDataExclDoublets.R
```
Scripts in this directory are provided for completeness; data is available from GEO code (path variables etc.) will have to be adjusted to fit other environments. Processed and filtered expression data (without annotations) will be generated in the dat subdirectory in a serialized single cell experiment file,  ```./data/scRNA_dat/dat/sce.rds```. 

**Iterative clustering reveals kidney structure**

In this step of the analysis we identify major cell types of the kidney and derive single cell experiment objects with the corresponding annotations. We also generate figures and tables reported in the manuscript.


*A. Whole Kidney Strucutre*

First we annotate whole kidney structure. The file ```wk_structure.R``` contains analyses to stratify cells and produce annotated object.

```
$ cd R && R --vanilla < ./wk_structure.R &> wk_structure_out.txt
$ cd ..
```

Will generate the file(s)

```
- ./results/sce_annotated.rds
```


*B. Nephron progenitor related cell-types and lineages*

This step will group nephron-progenitor lineages and also create the whole-kidney results.

 ```
$ cd R
$ R --vanilla < ./wk_np-select.R    &> wk-np-select_out.txt
$ R --vanilla < ./np-wk_annotate.R  &> np-wk_annotate_out.txt
$ cd ..
```

Results generated:

```
- ./results/sce_fully-annotated.rds    #- sce of all cells with annotation
- ./results/sce_np_fully-annotated.rds #- sce of NP-lineage cells with annotation
```

These single cell experiments have cluster annotations in their ```colData``` slots:

```
- $cluster      #- clustering of all kidney (10 clusters)
- $cluster_it   #- contains NP-lineage clustering results (11 clusters)
- $cluster_tme  #- mature vs. immature subclusters for NP-lineage cells
```

*C. Cell-type-/clustering-/annotation-related figures*

In this step we generate (most of) the plots/figures related to comparing cell-types.

```
$ cd R
$ R --vanilla < ./np-wk_plots.R &> np-wk_plots_out.txt
$ cd ..
```

The following figures are generated:

```
- ./figures/wk_celltype-tsne.pdf      #- Figure 1B
- ./figures/wk_celltype-violins.pdf   #- Figure 1C
- ./figures/np_celltype-tsne.pdf      #- Figure 2A
- ./figures/np_celltype-violins.pdf   #- Figure 2B
- ./figures/np_podocyte-violins.pdf   #- Figure 4A
- ./figures/np_tubular-violins.pdf    #- Figure 4B
- ./figures/xl_np-stromal_heatmap.pdf #- Figure 5
```

*D. Cell-type-/clustering-/annotation-related Tables*

Here we produce tables related to the previous steps:

```
$ cd R
$ R --vanilla < ./np-wk_tables.R &> np-wk_tables_out.txt
$ cd ..
```

The following results are generated:

```
- ./results/xl_np-stromal_coexpression.rds            #- Table 1
- ./results/xl_np-stromal_coexpression.xlsx           #- Table 1
- ./results/xl_tub-ub_expression.rds                  #- Table 2
- ./results/xl_tub-ub_expression.xlsx                 #- Table 2
- ./results/wk_celltype-DEG.rds                       #- Supplemental Table 1
- ./results/wk_celltype-DEG.xlsx                      #- Supplemental Table 1
- ./results/np_celltype-DEG.rds                       #- Supplemental Table 2
- ./results/np_celltype-DEG.xlsx                      #- Supplemental Table 2
- ./results/np_celltype-mature-vs-immature-DEG.rds    #- Supplemental Table 9
- ./results/np_celltype-mature-vs-immature-DEG.xlsx   #- Supplmeental Table 9
- ./results/np_celltype-tub-dist-vs-prox-DEG.rds      #- Supplemental Table 10
- ./results/np_celltype-tub-dist-vs-prox-DEG.xlsx     #- Supplemental Table 10
```

*E. Cell-type-/clustering-/annotation-related Supplemental Figures*

Here we focus on the supplemental figures included in the manuscript:

```
$ cd R
$ R --vanilla < ./np-wk_supplemental-plots.R &> np-wk_supplemental-plots_out.txt
$ cd ..
```
Produces:

```
- ./figures/wk_celltype-DEG-heatmap.pdf              #- Supplemental Figure 1
- ./figures/np_pseudotime-lineages.pdf               #- Supplemental Figure 2
- ./figures/np_celltype-DEG-heatmap.pdf              #- Supplemental Figure 3
- ./figures/np_celltype-tp-vs-td-DEG-heatmap.pdf     #- Supplemental Figure 5#A
- ./figures/np_celltype-mtd-vs-itd-DEG-heatmap.pdf   #- Supplemental Figure 5#B
- ./figures/np_celltype-mtp-vs-itp-DEG-heatmap.pdf   #- Supplemental Figure 5#C
- ./figures/np_celltype-mpod-vs-ipod-DEG-heatmap.pdf #- Supplemental Figure 5A
- ./figures/np-ub_tsne-clust.pdf                     #- Supplemental Figure 4A
- ./figures/np-ub_tsne-birc5.pdf                     #- Supplemental Figure 4B
- ./figures/np-ub_tsne-noleg-clust.pdf               #- Supplemental Figure 4A (no legend)
- ./figures/np-ub_tsne-noleg-birc5.pdf               #- Supplemental Figure 4A (no legend)
```

**Comparison of self-renewing and primed nephron progenitor cells**

Here we specifically focus on comparing nephron progenitor cell (sub)types "self-renew" and "primed".

*A. Differential Expression*

In this analysis we report differentially expressed genes (between self-renew vs primed nephron progenitor cells, and between primed nephron progenitors and differentiating cells) and we report enriched Gene Ontology terms.

```
$ cd R
$ R --vanilla < np_self-renew-vs-primed.R &> \
  np_self-renew-vs-primed_out.txt
$ cd ..
```

Produces the following tables in the ```./results``` direcotry

```
- np_self-vs-primed-DEG.rds              #- Supplemental Table 3
- np_self-vs-primed-DEG.xlsx             #- Supplemental Table 3
- np_primed-vs-differentiating-DEG.rds   #- Supplemental Table 4
- np_primed-vs-differentiating-DEG.xlsx  #- Supplemental Table 4
- np_np-clusters-go-combined.xlsx        #- Supplemental Table 5
- np_np-clusters-go-sr-vs-pr.rds         #- Supplemental Table 5
- np_np-clusters-go-pr-vs-diff.rds       #- Supplemental Table 5
```

... and these figures:

```
- ./figures/np_np-clusters-deg-heatmap.pdf       #- Figure 3A
- ./figures/np_np-clusters-go-sr-vs-pr.png       #- Figure 3B
- ./figures/np_np-clusters-go-pr-vs-diff.png     #- Figure 3C
- ./figures/np_pseudotime_birc5-rspo1-ccnd1.pdf" #- Figure 3E
```

*B. Regulatory modules*

Here we use pseudotime analysis and SCENIC to identify regulatory modules in the context of nephron progenitor "maturation" between self-renew and primed.

Running

```
$ mkdir -p ../results/tmp/SCENIC/int/
$ mkdir -p ../results/tmp/SCENIC/output/
$ cd R
$ R --vanilla < ./np_self-renew-vs-primed_network.R &> \
  np_self-renew-vs-primed_network_out.txt
$ cd ..
```

Generates
```
- ./figures/np_pseudotime-gene-set-heatmap.pdf          #- Figure 3g
- ./results/np_pseudotime-gene-set_regulons.rds         #- Supplemental Table 8
- ./results/np_pseudotime-gene-set_regulons.xlsx        #- Supplemental Table 8
- ./results/np_self-renew-primed-pseudotime-upDown.rds  #- Supplemental Table 6
- ./results/np_self-renew-primed-pseudotime-upDown.xlsx #- Supplemental Table 6
- ./results/np_self-renew-primed-pseudotime-go.xlsx     #- Supplemental Table 7
- ./results/np_self-renew-primed-pseudotime_go-dn.rds   #- Supplemental Table 7
- ./results/np_self-renew-primed-pseudotime_go-up.rds   #- Supplemental Table 7
```
