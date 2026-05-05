#!/usr/bin/env python3
"""
Normalization, HVG selection, scaling, and PCA for scRNA-seq data.
Supports scran, seurat-style (CPM), and pearson residuals normalization.
"""

import argparse
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',         required=True)
    parser.add_argument('--output',        required=True)
    parser.add_argument('--normalization', default='scran',
                        choices=['scran', 'seurat', 'pearson'])
    parser.add_argument('--n-top-genes',   type=int, default=2000)
    parser.add_argument('--scale-max',     type=float, default=10)
    parser.add_argument('--n-pcs',         type=int, default=30)
    return parser.parse_args()


def normalize_scran(adata):
    """Size-factor normalization via scran pooling."""
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()

    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added='groups')

    # Compute size factors with scran
    r('suppressMessages(library(scran))')
    data_mat = adata_pp.X.T.toarray() if hasattr(adata_pp.X, 'toarray') else adata_pp.X.T
    r_data = anndata2ri.py2r(adata_pp)
    size_factors = r('sizeFactors(computeSumFactors(as(assay({d},"X"),"dgCMatrix"), clusters={c}))'
                     .format(d='r_data', c='adata_pp.obs["groups"].values'))
    adata.obs['size_factors'] = np.array(size_factors)
    adata.X = adata.X / adata.obs['size_factors'].values[:, None]
    sc.pp.log1p(adata)
    return adata


def normalize_seurat(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def normalize_pearson(adata):
    sc.experimental.pp.normalize_pearson_residuals(adata)
    return adata


def main():
    args = parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    adata.raw = adata.copy()  # keep raw counts

    print(f"Normalization: {args.normalization}")
    if args.normalization == 'scran':
        try:
            adata = normalize_scran(adata)
        except Exception as e:
            print(f"scran failed ({e}), falling back to seurat normalization")
            adata = sc.read_h5ad(args.input)
            adata.raw = adata.copy()
            adata = normalize_seurat(adata)
    elif args.normalization == 'seurat':
        adata = normalize_seurat(adata)
    else:
        adata = normalize_pearson(adata)

    print(f"Selecting top {args.n_top_genes} highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, subset=False)

    # Plot HVG
    fig, ax = plt.subplots(figsize=(8, 5))
    sc.pl.highly_variable_genes(adata, ax=ax, show=False)
    plt.tight_layout()
    plt.savefig('hvg_plot.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Scale only HVG
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=args.scale_max)

    # PCA on HVG
    print(f"Running PCA ({args.n_pcs} components)...")
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, svd_solver='arpack')
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca']    = adata_hvg.uns['pca']
    adata.varm['PCs']   = np.zeros((adata.n_vars, args.n_pcs))
    adata.varm['PCs'][adata.var['highly_variable']] = adata_hvg.varm['PCs']

    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    print("Preprocessing complete.")


if __name__ == '__main__':
    main()
