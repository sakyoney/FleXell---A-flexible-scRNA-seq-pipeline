// One FASTQ (or pair) per cell — plate-based Smart-seq2
process STAR_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/alignment/smartseq2/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path star_index

    output:
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}_Log.final.out"),                  emit: log
    path "versions.yml",                                                 emit: versions

    script:
    def reads      = reads_2 ? "${reads_1} ${reads_2}" : "${reads_1}"
    def strandness = params.strandedness ?: 'unstranded'
    """
    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM \\
        --outFilterMultimapNmax 1 \\
        --outFileNamePrefix ${meta.id}_ \\
        --runThreadN ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """
}
