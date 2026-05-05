process SCANPY_MERGE {
    label 'process_medium'

    publishDir "${params.outdir}/merged", mode: params.publish_dir_mode

    input:
    path counts_list   // collected list of count dirs/files
    path samplesheet   // original CSV for metadata

    output:
    path "merged.h5ad", emit: h5ad
    path "versions.yml", emit: versions

    script:
    def protocol = params.protocol ?: 'mixed'
    """
    merge_samples.py \\
        --input ${counts_list} \\
        --metadata ${samplesheet} \\
        --protocol ${protocol} \\
        --output merged.h5ad

    python3 -c "import scanpy; print(scanpy.__version__)" > sc_version.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(cat sc_version.txt)
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
