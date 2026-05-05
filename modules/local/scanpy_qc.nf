process SCANPY_QC {
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "filtered.h5ad",   emit: h5ad
    path "qc_plots/*",      emit: plots
    path "qc_metrics.csv",  emit: metrics
    path "versions.yml",    emit: versions

    script:
    def doublet_flag = params.doublet_detection ? '--doublet-detection' : ''
    """
    mkdir -p qc_plots

    qc_filter.py \\
        --input ${h5ad} \\
        --output filtered.h5ad \\
        --outdir qc_plots \\
        --min-genes ${params.min_genes} \\
        --max-genes ${params.max_genes} \\
        --min-cells ${params.min_cells} \\
        --mt-percent-max ${params.mt_percent_max} \\
        ${doublet_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
