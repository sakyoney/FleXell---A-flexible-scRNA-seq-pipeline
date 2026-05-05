process GENERATE_REPORT {
    label 'process_low'

    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    input:
    path h5ad
    path markers
    path degs
    path multiqc_report

    output:
    path "flexell_report.html", emit: report
    path "versions.yml",        emit: versions

    script:
    """
    generate_report.py \\
        --h5ad ${h5ad} \\
        --markers ${markers} \\
        --degs ${degs} \\
        --multiqc ${multiqc_report} \\
        --output flexell_report.html \\
        --pipeline-version ${workflow.manifest.version}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
