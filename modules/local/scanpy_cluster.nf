process SCANPY_CLUSTER {
    label 'process_medium'

    publishDir "${params.outdir}/clustering", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "clustered.h5ad",  emit: h5ad
    path "umap_clusters.png", emit: umap
    path "versions.yml",    emit: versions

    script:
    """
    cluster.py \\
        --input ${h5ad} \\
        --output clustered.h5ad \\
        --method ${params.clustering_method} \\
        --resolution ${params.resolution} \\
        --n-neighbors ${params.n_neighbors} \\
        --min-dist ${params.min_dist} \\
        --n-pcs ${params.n_pcs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
        leidenalg: \$(python3 -c "import leidenalg; print(leidenalg.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """
}
