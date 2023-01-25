#!/usr/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy as sc

sc.settings.verbosity = 2 
sc.settings.set_figure_params(dpi=150)

path = 'filtered_feature_bc_matrix/'
adata = sc.read(path + 'matrix.mtx', cache=True).T  # first time only
adata.var_names = pd.read_csv(path + 'features.tsv', header=None, sep='\t')[1]
adata.var_names_make_unique()

adata.obs['n_counts'] = adata.X.sum(1).A1
sc.pp.normalize_per_cell(adata, counts_per_cell_after=adata.obs['n_counts'].mean())
adata.raw = adata

sc.pp.filter_genes(adata, min_cells=3)
filter_result = sc.pp.filter_genes_dispersion(
    adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
print('{} genes pass filter'.format(filter_result['gene_subset'].sum()))
#sc.pl.filter_genes_dispersion(filter_result)

adata = adata[:, filter_result.gene_subset]
sc.pp.scale(adata)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=5, use_rep='X_pca')
#sc.tl.umap(adata)

import scrublet as scr
scrub = scr.Scrublet(adata.raw.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
scrub.plot_histogram()
#sc.pl.umap(adata, color='doublet_scores', color_map='Reds', size=40)

threshold = 0.25
doub_obs = adata.obs['doublet_scores']
adata.obs['predicted_doublets'] = pd.Categorical(['doublet' if x > threshold else 'singlet' for x in doub_obs])
#sc.pl.umap(adata, color='predicted_doublets', color_map='Reds', size=40)

adata.write('Sample.h5ad')
