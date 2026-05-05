process SCANPY_DIFFEXP {
    label 'process_medium'

    publishDir "${params.outdir}/diffexp", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "degs.csv",            emit: degs
    path "volcano_plot.png",    emit: volcano
    path "versions.yml",        emit: versions

    script:
    """
    diffexp.py \\
        --input ${h5ad} \\
        --output degs.csv \\
        --method ${params.de_method} \\
        --logfc-threshold ${params.de_logfc_threshold} \\
        --pval-cutoff ${params.de_pval_cutoff} \\
        --groupby condition

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
