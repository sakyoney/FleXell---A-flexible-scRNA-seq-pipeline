process SCANPY_PREPROCESS {
    label 'process_medium'

    publishDir "${params.outdir}/preprocessed", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "preprocessed.h5ad", emit: h5ad
    path "hvg_plot.png",       emit: hvg_plot
    path "versions.yml",       emit: versions

    script:
    """
    preprocess.py \\
        --input ${h5ad} \\
        --output preprocessed.h5ad \\
        --normalization ${params.normalization} \\
        --n-top-genes ${params.n_top_genes} \\
        --scale-max ${params.scale_max_value} \\
        --n-pcs ${params.n_pcs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
