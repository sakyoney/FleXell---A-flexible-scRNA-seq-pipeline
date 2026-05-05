process SCANPY_MARKERS {
    label 'process_medium'

    publishDir "${params.outdir}/markers", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "markers.csv",         emit: markers
    path "dotplot_markers.png", emit: dotplot
    path "versions.yml",        emit: versions

    script:
    """
    markers.py \\
        --input ${h5ad} \\
        --output markers.csv \\
        --method ${params.marker_method} \\
        --min-logfc ${params.min_logfc} \\
        --pval-cutoff ${params.marker_pval}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
