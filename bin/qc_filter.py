#!/usr/bin/env python3
"""
Quality control and filtering for scRNA-seq data
"""

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='QC and filtering for scRNA-seq')
    parser.add_argument('--input', required=True, help='Input H5AD file')
    parser.add_argument('--output', required=True, help='Output filtered H5AD file')
    parser.add_argument('--outdir', required=True, help='Output directory for plots')
    parser.add_argument('--min-genes', type=int, default=200,
                       help='Minimum number of genes per cell')
    parser.add_argument('--max-genes', type=int, default=6000,
                       help='Maximum number of genes per cell')
    parser.add_argument('--min-cells', type=int, default=3,
                       help='Minimum number of cells per gene')
    parser.add_argument('--mt-percent-max', type=float, default=20,
                       help='Maximum mitochondrial percentage')
    parser.add_argument('--doublet-detection', action='store_true',
                       help='Enable doublet detection')
    return parser.parse_args()

def plot_qc_metrics(adata, outdir):
    """Generate QC plots"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Violin plots for QC metrics
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 
                         'pct_counts_mt', 'pct_counts_ribo'],
                jitter=0.4, multi_panel=True, ax=axes, show=False)
    
    plt.tight_layout()
    plt.savefig(outdir / 'qc_violin_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Scatter plots
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1], show=False)
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', ax=axes[2], show=False)
    
    plt.tight_layout()
    plt.savefig(outdir / 'qc_scatter_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Histogram of gene counts
    plt.figure(figsize=(10, 6))
    plt.hist(adata.obs['n_genes_by_counts'], bins=100, edgecolor='black')
    plt.xlabel('Number of genes')
    plt.ylabel('Number of cells')
    plt.title('Distribution of gene counts per cell')
    plt.savefig(outdir / 'gene_count_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()

def detect_doublets(adata):
    """Detect doublets using Scrublet"""
    try:
        import scrublet as scr
        
        print("Running doublet detection with Scrublet...")
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30
        )
        
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        
        print(f"Detected {predicted_doublets.sum()} doublets "
              f"({100*predicted_doublets.sum()/len(predicted_doublets):.2f}%)")
        
    except ImportError:
        print("Warning: Scrublet not installed. Skipping doublet detection.")
        adata.obs['doublet_score'] = 0
        adata.obs['predicted_doublet'] = False

def filter_cells_genes(adata, min_genes, max_genes, min_cells, mt_percent_max):
    """Apply filtering thresholds"""
    
    print("\nBefore filtering:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter by gene count upper limit
    adata = adata[adata.obs['n_genes_by_counts'] < max_genes, :].copy()
    
    # Filter by mitochondrial content
    adata = adata[adata.obs['pct_counts_mt'] < mt_percent_max, :].copy()
    
    print("\nAfter filtering:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median UMI/cell: {adata.obs['total_counts'].median():.0f}")
    print(f"  Median MT%: {adata.obs['pct_counts_mt'].median():.2f}%")
    
    return adata

def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    # Generate QC plots before filtering
    print("\nGenerating QC plots...")
    plot_qc_metrics(adata, args.outdir)
    
    # Doublet detection
    if args.doublet_detection:
        detect_doublets(adata)
        adata = adata[~adata.obs['predicted_doublet']].copy()
    
    # Apply filters
    print("\nApplying filters...")
    adata = filter_cells_genes(
        adata,
        args.min_genes,
        args.max_genes,
        args.min_cells,
        args.mt_percent_max
    )
    
    # Generate QC plots after filtering
    plot_qc_metrics(adata, Path(args.outdir) / 'filtered')
    
    # Save QC metrics
    qc_summary = adata.obs[[
        'n_genes_by_counts',
        'total_counts',
        'pct_counts_mt',
        'pct_counts_ribo'
    ]].describe()
    
    qc_summary.to_csv(Path(args.outdir) / 'qc_metrics.csv')
    
    # Save filtered data
    print(f"\nSaving filtered data to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    
    print("\nQC and filtering complete!")

if __name__ == '__main__':
    main()
