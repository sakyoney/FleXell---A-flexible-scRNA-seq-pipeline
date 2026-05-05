#!/usr/bin/env python3
"""
Differential expression between biological conditions (e.g. control vs disease).
Groups cells by 'condition' column within each cluster, runs DE per cluster.
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',           required=True)
    parser.add_argument('--output',          required=True)
    parser.add_argument('--method',          default='wilcoxon',
                        choices=['wilcoxon', 't-test', 'logreg'])
    parser.add_argument('--logfc-threshold', type=float, default=0.5)
    parser.add_argument('--pval-cutoff',     type=float, default=0.05)
    parser.add_argument('--groupby',         default='condition')
    return parser.parse_args()


def volcano_plot(df, output_path):
    """Simple volcano plot: -log10(pval_adj) vs log2FC."""
    df = df.dropna(subset=['logfoldchanges', 'pvals_adj'])
    df['neg_log10_pval'] = -np.log10(df['pvals_adj'].clip(lower=1e-300))

    sig = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'].abs() > 0.5)]
    not_sig = df[~df.index.isin(sig.index)]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(not_sig['logfoldchanges'], not_sig['neg_log10_pval'],
               alpha=0.3, s=5, color='grey', label='NS')
    ax.scatter(sig['logfoldchanges'], sig['neg_log10_pval'],
               alpha=0.6, s=8, color='firebrick', label='Significant')
    ax.axvline(-0.5, ls='--', lw=0.8, color='steelblue')
    ax.axvline(0.5,  ls='--', lw=0.8, color='steelblue')
    ax.axhline(-np.log10(0.05), ls='--', lw=0.8, color='grey')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adj. p-value)')
    ax.set_title('Differential Expression')
    ax.legend(markerscale=2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs.columns:
        raise ValueError(f"Column '{args.groupby}' not found in adata.obs. "
                         f"Available: {list(adata.obs.columns)}")

    conditions = adata.obs[args.groupby].unique().tolist()
    if len(conditions) < 2:
        raise ValueError(f"Need at least 2 groups in '{args.groupby}', found: {conditions}")

    print(f"Running DE: {args.method}, groupby='{args.groupby}' ({conditions})...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=args.groupby,
        method=args.method,
        use_raw=True if adata.raw is not None else False,
        pts=True
    )

    result = sc.get.rank_genes_groups_df(
        adata,
        group=None,
        key='rank_genes_groups',
        log2fc_min=None,
        pval_cutoff=None
    )

    # Filter and sort
    result_filtered = result[
        (result['logfoldchanges'].abs() >= args.logfc_threshold) &
        (result['pvals_adj'] <= args.pval_cutoff)
    ].sort_values('pvals_adj')

    result.to_csv(args.output, index=False)
    sig_path = args.output.replace('.csv', '_significant.csv')
    result_filtered.to_csv(sig_path, index=False)

    print(f"Total DE genes: {len(result_filtered)} "
          f"(|log2FC|≥{args.logfc_threshold}, padj≤{args.pval_cutoff})")

    volcano_plot(result, 'volcano_plot.png')
    print("DE analysis complete.")


if __name__ == '__main__':
    main()
