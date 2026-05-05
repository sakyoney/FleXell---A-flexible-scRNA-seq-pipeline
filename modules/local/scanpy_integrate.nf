process SCANPY_INTEGRATE {
    label 'process_high'

    publishDir "${params.outdir}/integrated", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "integrated.h5ad",   emit: h5ad
    path "umap_batch.png",    emit: umap
    path "versions.yml",      emit: versions

    script:
    """
    integrate.py \\
        --input ${h5ad} \\
        --output integrated.h5ad \\
        --method ${params.integration_method} \\
        --batch-key ${params.batch_key} \\
        --n-pcs ${params.n_pcs} \\
        --n-neighbors ${params.n_neighbors} \\
        --min-dist ${params.min_dist}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
        harmonypy: \$(python3 -c "import harmonypy; print(harmonypy.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """
}
