process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/quantification/smartseq2/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    tuple val(meta), path("${meta.id}_counts.txt"), emit: counts
    tuple val(meta), path("${meta.id}_counts.txt.summary"), emit: summary
    path "versions.yml", emit: versions

    script:
    def strandedness = params.strandedness == 'forward'  ? 1 :
                       params.strandedness == 'reverse'  ? 2 : 0
    def paired_flag  = params.smartseq_paired ? '-p' : ''
    """
    featureCounts \\
        -a ${gtf} \\
        -o ${meta.id}_counts.txt \\
        -T ${task.cpus} \\
        -s ${strandedness} \\
        ${paired_flag} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oP 'v\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}
