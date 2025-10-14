#!/usr/bin/env python3
"""
Merge multiple scRNA-seq samples into a single AnnData object
Supports both 10x and Smart-seq2 protocols
"""

import argparse
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Merge scRNA-seq samples')
    parser.add_argument('--input', nargs='+', required=True,
                       help='Input count directories or files')
    parser.add_argument('--metadata', required=True,
                       help='Sample metadata CSV file')
    parser.add_argument('--protocol', choices=['10x', 'smartseq2', 'mixed'],
                       default='mixed', help='Sequencing protocol')
    parser.add_argument('--output', required=True,
                       help='Output H5AD file')
    return parser.parse_args()

def read_10x_data(sample_path, sample_id):
    """Read 10x Genomics data"""
    print(f"Reading 10x data from {sample_path}...")
    
    # Cell Ranger output structure
    if os.path.isdir(os.path.join(sample_path, 'filtered_feature_bc_matrix')):
        adata = sc.read_10x_mtx(
            os.path.join(sample_path, 'filtered_feature_bc_matrix'),
            var_names='gene_symbols',
            cache=True
        )
    else:
        adata = sc.read_10x_mtx(
            sample_path,
            var_names='gene_symbols',
            cache=True
        )
    
    # Add sample ID to cell barcodes
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
    adata.obs['sample_id'] = sample_id
    
    return adata

def read_smartseq2_data(count_file, sample_id):
    """Read Smart-seq2 count data"""
    print(f"Reading Smart-seq2 data from {count_file}...")
    
    # Read featureCounts output
    counts = pd.read_csv(count_file, sep='\t', comment='#', index_col=0)
    
    # Extract count column (usually last column)
    count_data = counts.iloc[:, -1:]
    count_data.columns = [sample_id]
    
    # Create AnnData
    adata = sc.AnnData(
        X=count_data.T,
        obs=pd.DataFrame(index=[sample_id]),
        var=pd.DataFrame(index=count_data.index)
    )
    
    adata.obs['sample_id'] = sample_id
    
    return adata

def merge_samples(input_paths, metadata_file, protocol):
    """Merge multiple samples into single AnnData object"""
    
    # Read metadata
    metadata = pd.read_csv(metadata_file)
    metadata = metadata.set_index('sample_id')
    
    adata_list = []
    
    for input_path in input_paths:
        # Extract sample ID from path
        sample_id = Path(input_path).stem
        if sample_id not in metadata.index:
            sample_id = Path(input_path).parent.name
        
        if sample_id not in metadata.index:
            print(f"Warning: {sample_id} not found in metadata, skipping...")
            continue
        
        sample_protocol = metadata.loc[sample_id, 'protocol'] if 'protocol' in metadata.columns else protocol
        
        # Read data based on protocol
        if sample_protocol == '10x':
            adata = read_10x_data(input_path, sample_id)
        elif sample_protocol == 'smartseq2':
            adata = read_smartseq2_data(input_path, sample_id)
        else:
            print(f"Unknown protocol {sample_protocol} for {sample_id}, skipping...")
            continue
        
        # Add metadata
        for col in metadata.columns:
            if col != 'protocol':
                adata.obs[col] = metadata.loc[sample_id, col]
        
        adata_list.append(adata)
    
    if not adata_list:
        raise ValueError("No valid samples found!")
    
    # Concatenate all samples
    print(f"\nMerging {len(adata_list)} samples...")
    adata = sc.concat(
        adata_list,
        join='outer',
        merge='same',
        label='batch',
        keys=[adata.obs['sample_id'].iloc[0] for adata in adata_list]
    )
    
    # Make variable names unique
    adata.var_names_make_unique()
    
    # Add basic metrics
    print("\nCalculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    print(f"\nMerged dataset:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Samples: {adata.obs['sample_id'].nunique()}")
    
    if 'condition' in adata.obs.columns:
        print(f"\nConditions:")
        print(adata.obs['condition'].value_counts())
    
    return adata

def main():
    args = parse_args()
    
    # Merge samples
    adata = merge_samples(args.input, args.metadata, args.protocol)
    
    # Save merged data
    print(f"\nSaving merged data to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    
    print("\nDone!")

if __name__ == '__main__':
    main()
