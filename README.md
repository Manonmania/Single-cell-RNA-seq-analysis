![Build Status](https://gitlab.com/pages/plain-html/badges/master/build.svg)

---

## Single Cell RNA-seq data analysis tutorial

---

This repository contains scripts and documentation associated with the analyses of Single Cell RNA-seq data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- Doublets detection
- Setup Seurat Object
- Quality control
- Normalizing the data
- Identify variable features
- Scaling the data
- Dimentionality reduction
- Batch correction
- Cluster the cells
- Run non-linear dimensional reduction UMAP/tSNE
- Celltype identification
- Finding differentially expressed features
- Pathway analysis
- Scripts

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Single Cell RNA-seq data analysis

A step-by-step description of how to analyze Single Cell RNA-seq data using **Scrublet**, **Seurat**, **Harmony**, **SingleR**, and **FGSEA** tools.

- **STEP 1: (Doublets detction)**

  Using Scrublet we identify doublets in single-cell RNA-seq data. Scrublet calculates a doublet score for each cell for given a raw (unnormalized) UMI counts matrix counts_matrix with cells as rows and genes as columns.


- **STEP 2:(Setup Seurat Object)**

  Following doublets detection using Scrublet, the doublets are removed and only singlets considered for further analysis. We setup seurat object by loading and merging the data using "CreateSeuratObject" function.


- **STEP 3:(Quality control)**

  In this step, we explore QC metrics and filter cells based on the filtering criteria. Low-quality cells or empty droplets will often have very few genes and those cells are filtered out. Low-quality and dying cells often exhibit extensive mitochondrial contamination and those cells are also filtered out. We filter cells that have unique feature counts over 2,500 or less than 200. We filter cells that have >5% mitochondrial counts.


- **STEP 4:(Normalizing the data)**

  After removing unwanted cells from the dataset, the next step is to normalize the data. By default, Seurat apply a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.

- **STEP 5:(Identify variable features)**

  As a next step, Seurat calculates a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
  

- **STEP 6:(Scaling the data)**

  After identifying the variable features, we scale the data that is a standard pre-processing step prior to dimensional reduction techniques like PCA.


- **STEP 7:(Dimentionality reduction)**

  Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.


- **STEP 8:(Batch correction)**

  If the dataset come from different batches then we need to do a batch correction using Harmony. The most common way to run Harmony is on reduced dimensions such as PC embeddings from principal component analysis (PCA).

 
- **STEP 9:(Cluster the cells)**

  As a next step we cluster the cells. Seurat applies modularity optimization techniques such as the Louvain algorithm (default) or SLM to iteratively group cells together, with the goal of optimizing the standard modularity function. 

- **STEP 10: (Run non-linear dimensional reduction UMAP/tSNE)**

  Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space.

- **STEP 11: (Celltype identifcation)**

  We used SingleR a computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. SingleR’s annotations combined with Seurat, a processing and analysis package designed for scRNA-seq.

- **STEP 12: (Finding differentially expressed features)**

  Seurat find markers that define clusters via differential expression. By default, it identifies positive and negative markers of a single cluster compared to all other cells. FindAllMarkers() automates this process for all clusters.

- **STEP 13: (Pathway analysis)**

  We used FGSEA, an R-package for fast preranked gene set enrichment analysis (GSEA), for pathway analysis. This package allows calculate arbitrarily low GSEA P-values for a collection of gene sets.

## Scripts

- Doublets detection using Scrublet
  - Doublets_Detection_Scrublet.py

- Single Cell RNA-seq analysis (Human and Mouse)
  - Single-cell-RNAseq-analysis-Human.Rmd
  - Single-cell-RNAseq-analysis-Mouse.Rmd

- Single Cell RNA-seq analysis report

