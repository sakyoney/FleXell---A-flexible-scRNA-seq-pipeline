process FASTQC {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/fastqc/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path "versions.yml",             emit: versions

    script:
    def reads = reads_2 ? "${reads_1} ${reads_2}" : "${reads_1}"
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
