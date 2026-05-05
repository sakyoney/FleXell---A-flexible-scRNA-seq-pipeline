#!/usr/bin/env python3
"""
Generate a standalone HTML summary report for a FleXell pipeline run.
Embeds UMAP, QC, marker, and DE plots directly into the HTML.
"""

import argparse
import base64
import io
import os
from datetime import datetime

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad',             required=True)
    parser.add_argument('--markers',          required=True)
    parser.add_argument('--degs',             required=True)
    parser.add_argument('--multiqc',          required=True)
    parser.add_argument('--output',           required=True)
    parser.add_argument('--pipeline-version', default='1.0.0')
    return parser.parse_args()


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def make_umap_figure(adata):
    color_keys = [k for k in ['cluster', 'celltypist_majority_vote', 'condition', 'sample_id']
                  if k in adata.obs.columns]
    if not color_keys or 'X_umap' not in adata.obsm:
        return None
    ncols = min(len(color_keys), 3)
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    if ncols == 1:
        axes = [axes]
    for ax, key in zip(axes, color_keys[:ncols]):
        sc.pl.umap(adata, color=key, ax=ax, show=False, title=key)
    plt.tight_layout()
    return fig


def make_qc_figure(adata):
    metrics = [m for m in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
               if m in adata.obs.columns]
    if not metrics:
        return None
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics):
        ax.hist(adata.obs[m], bins=60, color='steelblue', edgecolor='none', alpha=0.8)
        ax.set_xlabel(m)
        ax.set_ylabel('Cells')
        ax.set_title(m.replace('_', ' '))
    plt.tight_layout()
    return fig


def summary_table(adata):
    rows = [
        ('Total cells', adata.n_obs),
        ('Total genes', adata.n_vars),
    ]
    if 'sample_id' in adata.obs:
        rows.append(('Samples', adata.obs['sample_id'].nunique()))
    if 'condition' in adata.obs:
        rows.append(('Conditions', adata.obs['condition'].nunique()))
    if 'cluster' in adata.obs:
        rows.append(('Clusters', adata.obs['cluster'].nunique()))
    if 'celltypist_majority_vote' in adata.obs:
        rows.append(('Cell types', adata.obs['celltypist_majority_vote'].nunique()))
    return rows


def html_table(rows):
    lines = ['<table class="summary"><tr><th>Metric</th><th>Value</th></tr>']
    for k, v in rows:
        lines.append(f'<tr><td>{k}</td><td>{v}</td></tr>')
    lines.append('</table>')
    return '\n'.join(lines)


def csv_preview(path, n=10):
    try:
        df = pd.read_csv(path, nrows=n)
        return df.to_html(index=False, classes='data-table', border=0)
    except Exception:
        return '<p>Could not load table.</p>'


CSS = """
body { font-family: Arial, sans-serif; margin: 40px; color: #222; }
h1 { color: #2c3e50; }
h2 { color: #34495e; border-bottom: 2px solid #bdc3c7; padding-bottom: 4px; }
.summary { border-collapse: collapse; margin-bottom: 20px; }
.summary td, .summary th { border: 1px solid #ccc; padding: 6px 14px; }
.summary th { background: #ecf0f1; }
.data-table { border-collapse: collapse; font-size: 12px; margin-bottom: 20px; }
.data-table td, .data-table th { border: 1px solid #ddd; padding: 4px 10px; }
.data-table th { background: #f5f5f5; }
img { max-width: 100%; margin: 10px 0; border: 1px solid #ddd; border-radius: 4px; }
.footer { color: #999; font-size: 12px; margin-top: 40px; }
"""


def main():
    args = parse_args()

    print(f"Loading {args.h5ad}...")
    adata = sc.read_h5ad(args.h5ad)

    sections = []

    # --- Summary ---
    rows = summary_table(adata)
    sections.append(f'<h2>Pipeline Summary</h2>{html_table(rows)}')

    # --- QC ---
    qc_fig = make_qc_figure(adata)
    if qc_fig:
        enc = fig_to_base64(qc_fig)
        plt.close(qc_fig)
        sections.append(f'<h2>QC Metrics</h2><img src="data:image/png;base64,{enc}"/>')

    # --- UMAP ---
    umap_fig = make_umap_figure(adata)
    if umap_fig:
        enc = fig_to_base64(umap_fig)
        plt.close(umap_fig)
        sections.append(f'<h2>UMAP</h2><img src="data:image/png;base64,{enc}"/>')

    # --- Top markers ---
    sections.append(f'<h2>Top Marker Genes (preview)</h2>{csv_preview(args.markers)}')

    # --- Top DEGs ---
    sections.append(f'<h2>Differentially Expressed Genes (preview)</h2>{csv_preview(args.degs)}')

    # --- MultiQC link ---
    multiqc_name = os.path.basename(args.multiqc)
    sections.append(f'<h2>Sequencing QC</h2>'
                    f'<p>MultiQC report: <a href="{multiqc_name}">{multiqc_name}</a></p>')

    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    html = f"""<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>FleXell Report</title>
<style>{CSS}</style></head>
<body>
<h1>FleXell scRNA-seq Pipeline Report</h1>
<p>Generated: {now} &nbsp;|&nbsp; Pipeline version: {args.pipeline_version}</p>
{''.join(sections)}
<div class="footer">FleXell &mdash; Flexible single-cell RNA-seq pipeline</div>
</body></html>"""

    with open(args.output, 'w') as f:
        f.write(html)
    print(f"Report written to {args.output}")


if __name__ == '__main__':
    main()
