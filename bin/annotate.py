#!/usr/bin/env python3
"""
Automated cell type annotation using CellTypist.
Downloads model if not present, runs majority-vote prediction per cluster.
"""

import argparse
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',  required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--model',  default='Immune_All_Low.pkl',
                        help='CellTypist model name or local path')
    return parser.parse_args()


def main():
    args = parse_args()
    import celltypist
    from celltypist import models

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    # CellTypist expects log-normalized counts (10k per cell)
    adata_ct = adata.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    print(f"Downloading/loading CellTypist model: {args.model}...")
    try:
        model = models.Model.load(model=args.model)
    except Exception:
        models.download_models(model=args.model)
        model = models.Model.load(model=args.model)

    print("Running CellTypist prediction...")
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,
        over_clustering='cluster' if 'cluster' in adata.obs.columns else None
    )

    # Transfer annotations back to original adata
    adata.obs['celltypist_cell_type']        = predictions.predicted_labels['predicted_labels'].values
    adata.obs['celltypist_majority_vote']    = predictions.predicted_labels.get(
        'majority_voting', predictions.predicted_labels['predicted_labels']).values
    adata.obs['celltypist_conf_score']       = predictions.probability_matrix.max(axis=1).values

    # Save cell type table
    ct_table = adata.obs[['celltypist_cell_type',
                           'celltypist_majority_vote',
                           'celltypist_conf_score']].copy()
    if 'cluster' in adata.obs.columns:
        ct_table['cluster'] = adata.obs['cluster']
    ct_table.to_csv('celltypes.csv')

    # UMAP colored by cell type
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata, color='celltypist_majority_vote',
                   ax=axes[0], show=False, title='Cell type (majority vote)',
                   legend_loc='right margin')
        sc.pl.umap(adata, color='celltypist_conf_score',
                   ax=axes[1], show=False, title='Confidence score',
                   color_map='viridis')
        plt.tight_layout()
        plt.savefig('umap_celltypes.png', dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.figure()
        plt.savefig('umap_celltypes.png')
        plt.close()

    n_types = adata.obs['celltypist_majority_vote'].nunique()
    print(f"Annotated {adata.n_obs} cells into {n_types} cell types.")
    print(adata.obs['celltypist_majority_vote'].value_counts().to_string())

    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    print("Annotation complete.")


if __name__ == '__main__':
    main()
