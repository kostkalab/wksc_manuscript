

###Single cell RNA sequencing reveals differential cell cycle activity in key cell populations during nephrogenesis

<p align=center>
Abha S. Bais, Débora M. Cerqueira, Andrew  Clugston, Andrew J. Bodnar,
Jacqueline Ho and Dennis Kostka
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


**Single cell expression data**

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

First we annotate whole kidney structure. The file ```00_wk_structure.R``` contains analyses to stratify cells and produce annotated object.

```
$ cd R && R --vanilla < ./00_wk_structure.R &> 00_wk_structure_out.txt
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
$ R --vanilla < ./01_wk_np-select.R    &> 01_wk-np-select_out.txt
$ R --vanilla < ./02_np-wk_annotate.R  &> 02_np-wk_annotate_out.txt
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
$ R --vanilla < ./03_np-wk_plots.R &> 03_np-wk_plots_out.txt
$ cd ..
```

The following figures are generated:

```
- ./figures/fig1b_wk_celltype-tsne.pdf      #- Figure 1B
- ./figures/fig1c_wk_celltype-violins.pdf   #- Figure 1C
- ./figures/fig2a_np_celltype-tsne.pdf      #- Figure 2A
- ./figures/fig2b_np_celltype-violins.pdf   #- Figure 2B
- ./figures/fig4a_np_tubular-violins.pdf    #- Figure 4A
- ./figures/fig4b_np_podocyte-violins.pdf   #- Figure 4B
- ./figures/fig5_xl_np-stromal_heatmap.pdf #- Figure 5
```

*D. Cell-type-/clustering-/annotation-related Tables*

Here we produce tables related to the previous steps:

```
$ cd R
$ R --vanilla < ./04_np-wk_tables.R &> 04_np-wk_tables_out.txt
$ cd ..
```

The following results are generated:

```
- ./results/xl_np-stromal_coexpression.rds                    #- Table 1
- ./results/tab1_xl_np-stromal_coexpression.xlsx              #- Table 1
- ./results/xl_tub-ub_expression.rds                          #- Table 2
- ./results/tab2_xl_tub-ub_expression.xlsx                    #- Table 2
- ./results/wk_celltype-DEG.rds                               #- Supplemental Table 1
- ./results/supTab1_wk_celltype-DEG.xlsx                      #- Supplemental Table 1
- ./results/np_celltype-DEG.rds                               #- Supplemental Table 2
- ./results/supTab2_np_celltype-DEG.xlsx                      #- Supplemental Table 2
- ./results/np_celltype-mature-vs-immature-DEG.rds            #- Supplemental Table 9
- ./results/supTab9_np_celltype-mature-vs-immature-DEG.xlsx   #- Supplmeental Table 9
- ./results/np_celltype-tub-dist-vs-prox-DEG.rds              #- Supplemental Table 10
- ./results/supTab10_np_celltype-tub-dist-vs-prox-DEG.xlsx    #- Supplemental Table 10
```


*E. Cell-type-/clustering-/annotation-related Supplemental Figures*

Here we focus on the supplemental figures included in the manuscript:

```
$ cd R
$ R --vanilla < ./05_np-wk_supplemental-plots.R &> 05_np-wk_supplemental-plots_out.txt
$ cd ..
```
Produces:

```
- ./figures/supFig1_wk_celltype-DEG-heatmap.pdf               #- Supplemental Figure 1
- ./figures/supFig2_np_pseudotime-lineages.pdf                #- Supplemental Figure 2
- ./figures/supFig3_np_celltype-DEG-heatmap.pdf               #- Supplemental Figure 3
- ./figures/supFig4a_np_celltype-tp-vs-td-DEG-heatmap.pdf     #- Supplemental Figure 4#a
- ./figures/supFig4b_np_celltype-mtd-vs-itd-DEG-heatmap.pdf   #- Supplemental Figure 4#b
- ./figures/supFig4c_np_celltype-mtp-vs-itp-DEG-heatmap.pdf   #- Supplemental Figure 4#c
- ./figures/supFig4d_np_celltype-mpod-vs-ipod-DEG-heatmap.pdf #- Supplemental Figure 4d
- ./figures/supFig5a_np-ub_tsne-clust.pdf                     #- Supplemental Figure 5a
- ./figures/supFig5b_np-ub_tsne-birc5.pdf                     #- Supplemental Figure 5b
- ./figures/supFig5a_np-ub_tsne-noleg-clust.pdf               #- Supplemental Figure 5a (no legend)
- ./figures/supFig5b_np-ub_tsne-noleg-birc5.pdf               #- Supplemental Figure 5b (no legend)
```



**Comparison of self-renewing and primed nephron progenitor cells**

Here we specifically focus on comparing nephron progenitor cell (sub)types "self-renew" and "primed".

*A. Differential Expression*

In this analysis we report differentially expressed genes (between self-renew vs primed nephron progenitor cells, and between primed nephron progenitors and differentiating cells) and we report enriched Gene Ontology terms.

```
$ cd R
$ R --vanilla < 06_np_self-renew-vs-primed.R &> \
  06_np_self-renew-vs-primed_out.txt
$ cd ..
```

Produces the following tables in the ```./results``` direcotry

```
- np_self-vs-primed-DEG.rds                      #- Supplemental Table 3
- supTab3_np_self-vs-primed-DEG.xlsx             #- Supplemental Table 3
- np_primed-vs-differentiating-DEG.rds           #- Supplemental Table 4
- supTab4_np_primed-vs-differentiating-DEG.xlsx  #- Supplemental Table 4
- supTab5_np_np-clusters-go-combined.xlsx        #- Supplemental Table 5
- np_np-clusters-go-sr-vs-pr.rds                 #- Supplemental Table 5
- np_np-clusters-go-pr-vs-diff.rds               #- Supplemental Table 5
```

... and these figures:

```
- ./figures/fig3a_np_np-clusters-deg-heatmap.pdf       #- Figure 3A
- ./figures/fig3b_np_np-clusters-go-sr-vs-pr.png       #- Figure 3B
- ./figures/fig3c_np_np-clusters-go-pr-vs-diff.png     #- Figure 3C
- ./figures/fig3e_np_pseudotime_birc5-rspo1-ccnd1.pdf" #- Figure 3E
```

*B. Regulatory modules*

Here we use pseudotime analysis and SCENIC to identify regulatory modules in the context of nephron progenitor "maturation" between self-renew and primed.

Running

```
$ mkdir -p ../results/tmp/SCENIC/int/
$ mkdir -p ../results/tmp/SCENIC/output/
$ cd R
$ R --vanilla < ./07_np_self-renew-vs-primed_network.R &> \
  07_np_self-renew-vs-primed_network_out.txt
$ cd ..
```

Generates
```
- ./figures/fig3g_np_pseudotime-gene-set-heatmap.pdf            #- Figure 3g
- ./results/np_pseudotime-gene-set_regulons.rds                 #- Supplemental Table 8
- ./results/supTab8_np_pseudotime-gene-set_regulons.xlsx        #- Supplemental Table 8
- ./results/np_self-renew-primed-pseudotime-upDown.rds          #- Supplemental Table 6
- ./results/supTab6_np_self-renew-primed-pseudotime-upDown.xlsx #- Supplemental Table 6
- ./results/supTab7_np_self-renew-primed-pseudotime-go.xlsx     #- Supplemental Table 7
- ./results/np_self-renew-primed-pseudotime_go-dn.rds   #- Supplemental Table 7
- ./results/np_self-renew-primed-pseudotime_go-up.rds   #- Supplemental Table 7
```

*C. Miscellaneous*
This includes analyses related to miscellaneous points including Magella et al. dataset and cell cycle effects.
Note for Magella et al analyses, data needs to be downloaded from GSE104396 and syn11027925, and placed under ```../data/external/magella/```

```
$ cd R
$ R --vanilla < 09_magella.R &> \
  09_magella_out.txt
$ R --vanilla < 10_cellCycle.R &> \
  10_cellCycle_out.txt
$ R --vanilla < 11_misc.R &> \
  11_misc_out.txt
$ cd ..
```

```
- ./figures/supFig6a_magella-comp.pdf                         #- Supplemental Figure 6a
- ./figures/supFig6b_magella_tsne-np-birc5_leg.pdf            #- Supplemental Figure 6b
- ./figures/supFig6b_magella_tsne-np-birc5_noLeg.pdf          #- Supplemental Figure 6b (no legend)
- ./figures/supFig6c_magella_tsne-np-celltypes_leg.pdf        #- Supplemental Figure 6c
- ./figures/supFig6c_magella_tsne-np-celltypes_noLeg.pdf      #- Supplemental Figure 6c (no legend)
- ./figures/supFig7_wk_np_CCPhase-tsne.pdf                    #- Supplemental Figure 7 (bottom row) 
- ./figures/supFig7_wk_np_CCG2MScore-tsne.pdf                 #- Supplemental Figure 7 (bottom row)
- ./figures/supFig7_wk_np_CCSScore-tsne.pdf                   #- Supplemental Figure 7 (bottom row)
-./figures/supFig7_wk_np_CCScoresReg1-tsne.pdf                #- Supplemental Figure 7 (top right)
- ./figures/supFig8_deg_primed-vs-self_down.pdf                        #- Supplemental Figure 8
- ./figures/supFig8_deg_primed-vs-self_up.pdf                          #- Supplemental Figure 8
- ./results/supTab11_np_self-vs-primed-DEG_CCScores_lfc1.2.xlsx        #- Supplemental Table 11
- ./results/suppTab12_np_self-vs-primed-DEG_CCRmPhaseGenes_lfc1.2.xlsx #- Supplemental Table 12
- ./results/supTab13_birc5CorGenes_dist_tub_ub.xlsx                    #- Supplemetal Table 13
```

pdf("../figures/supFig6a_magella-comp.pdf", width = 8, height = 8)
