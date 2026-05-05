#!/usr/bin/env python3
"""
Per-cluster marker gene detection using scanpy rank_genes_groups.
"""

import argparse
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',       required=True)
    parser.add_argument('--output',      required=True)
    parser.add_argument('--method',      default='wilcoxon',
                        choices=['wilcoxon', 't-test', 'logreg'])
    parser.add_argument('--min-logfc',   type=float, default=0.5)
    parser.add_argument('--pval-cutoff', type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    if 'cluster' not in adata.obs.columns:
        raise ValueError("No 'cluster' column found. Run clustering first.")

    print(f"Finding marker genes per cluster ({args.method})...")
    sc.tl.rank_genes_groups(
        adata,
        groupby='cluster',
        method=args.method,
        use_raw=True if adata.raw is not None else False,
        pts=True
    )

    # Extract and filter results
    result = sc.get.rank_genes_groups_df(
        adata,
        group=None,
        key='rank_genes_groups',
        log2fc_min=args.min_logfc,
        pval_cutoff=args.pval_cutoff
    )

    result.to_csv(args.output, index=False)
    print(f"Markers saved to {args.output} ({len(result)} entries).")

    # Dotplot of top 5 markers per cluster
    top_markers = (result.groupby('group')
                         .head(5)['names']
                         .drop_duplicates()
                         .tolist())
    if top_markers:
        fig, ax = plt.subplots(figsize=(12, max(4, len(top_markers) * 0.3)))
        sc.pl.dotplot(adata, top_markers, groupby='cluster',
                      ax=ax, show=False, standard_scale='var')
        plt.tight_layout()
        plt.savefig('dotplot_markers.png', dpi=150, bbox_inches='tight')
        plt.close()
    else:
        # Create empty placeholder so the output glob doesn't fail
        plt.figure()
        plt.savefig('dotplot_markers.png')
        plt.close()

    print("Marker detection complete.")


if __name__ == '__main__':
    main()
