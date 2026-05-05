#!/usr/bin/env python3
"""
Clustering (Leiden default, Louvain optional) and UMAP visualization.
Expects X_pca or X_integrated in obsm; builds neighbor graph if absent.
"""

import argparse
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',       required=True)
    parser.add_argument('--output',      required=True)
    parser.add_argument('--method',      default='leiden', choices=['leiden', 'louvain'])
    parser.add_argument('--resolution',  type=float, default=0.8)
    parser.add_argument('--n-neighbors', type=int,   default=15)
    parser.add_argument('--min-dist',    type=float, default=0.3)
    parser.add_argument('--n-pcs',       type=int,   default=30)
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Build neighbor graph if not already present (e.g. no integration step)
    if 'neighbors' not in adata.uns:
        use_rep = 'X_integrated' if 'X_integrated' in adata.obsm else 'X_pca'
        print(f"Building neighbor graph (use_rep={use_rep})...")
        sc.pp.neighbors(adata, use_rep=use_rep,
                        n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    if 'X_umap' not in adata.obsm:
        print("Computing UMAP...")
        sc.tl.umap(adata, min_dist=args.min_dist)

    print(f"Clustering: {args.method} (resolution={args.resolution})...")
    if args.method == 'leiden':
        sc.tl.leiden(adata, resolution=args.resolution, key_added='cluster')
    else:
        sc.tl.louvain(adata, resolution=args.resolution, key_added='cluster')

    n_clusters = adata.obs['cluster'].nunique()
    print(f"Found {n_clusters} clusters.")

    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color='cluster', ax=ax, show=False,
               title=f'{args.method.capitalize()} clusters (res={args.resolution})',
               legend_loc='on data')
    plt.tight_layout()
    plt.savefig('umap_clusters.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    print("Clustering complete.")


if __name__ == '__main__':
    main()
