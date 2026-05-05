#!/usr/bin/env python3
"""
Batch integration: Harmony (default) or scVI (optional).
Adds corrected embedding to obsm and runs UMAP for visualization.
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
    parser.add_argument('--method',      default='harmony',
                        choices=['harmony', 'scvi', 'bbknn', 'combat', 'scanorama'])
    parser.add_argument('--batch-key',   default='sample_id')
    parser.add_argument('--n-pcs',       type=int, default=30)
    parser.add_argument('--n-neighbors', type=int, default=15)
    parser.add_argument('--min-dist',    type=float, default=0.3)
    return parser.parse_args()


def run_harmony(adata, batch_key, n_pcs):
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm['X_pca'][:, :n_pcs],
                        adata.obs, batch_key)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    adata.obsm['X_integrated']  = adata.obsm['X_pca_harmony']
    return adata


def run_scvi(adata, batch_key):
    import scvi
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    model.train(max_epochs=400, early_stopping=True)
    adata.obsm['X_scvi']       = model.get_latent_representation()
    adata.obsm['X_integrated'] = adata.obsm['X_scvi']
    return adata


def run_bbknn(adata, batch_key, n_pcs, n_neighbors):
    import bbknn
    sc.pp.neighbors(adata, use_rep='X_pca', n_pcs=n_pcs, n_neighbors=n_neighbors)
    bbknn.bbknn(adata, batch_key=batch_key)
    adata.obsm['X_integrated'] = adata.obsm['X_pca']
    return adata


def run_combat(adata, batch_key):
    sc.pp.combat(adata, key=batch_key)
    sc.tl.pca(adata)
    adata.obsm['X_integrated'] = adata.obsm['X_pca']
    return adata


def run_scanorama(adata, batch_key, n_pcs):
    import scanorama
    batches = [adata[adata.obs[batch_key] == b].copy()
               for b in adata.obs[batch_key].unique()]
    corrected, _ = scanorama.correct_scanpy(batches, return_dimred=True)
    import numpy as np
    adata.obsm['X_scanorama']  = np.vstack([c.obsm['X_scanorama'] for c in corrected])
    adata.obsm['X_integrated'] = adata.obsm['X_scanorama']
    return adata


def main():
    args = parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    print(f"Batch integration: {args.method} (batch_key={args.batch_key})")
    if args.method == 'harmony':
        adata = run_harmony(adata, args.batch_key, args.n_pcs)
    elif args.method == 'scvi':
        adata = run_scvi(adata, args.batch_key)
    elif args.method == 'bbknn':
        adata = run_bbknn(adata, args.batch_key, args.n_pcs, args.n_neighbors)
    elif args.method == 'combat':
        adata = run_combat(adata, args.batch_key)
    elif args.method == 'scanorama':
        adata = run_scanorama(adata, args.batch_key, args.n_pcs)

    use_rep = 'X_integrated' if 'X_integrated' in adata.obsm else 'X_pca'
    sc.pp.neighbors(adata, use_rep=use_rep,
                    n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    sc.tl.umap(adata, min_dist=args.min_dist)

    # UMAP colored by batch
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=args.batch_key, ax=ax, show=False,
               title=f'UMAP after {args.method} integration')
    plt.tight_layout()
    plt.savefig('umap_batch.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    print("Integration complete.")


if __name__ == '__main__':
    main()
