process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}/quantification/10x/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path star_index

    output:
    tuple val(meta), path("${meta.id}/Solo.out/GeneFull/filtered"), emit: counts
    tuple val(meta), path("${meta.id}/*.bam"),                      emit: bam
    tuple val(meta), path("${meta.id}/Log.final.out"),              emit: log
    path "versions.yml",                                             emit: versions

    script:
    def whitelist = params.starsolo_whitelist ? "--soloCBwhitelist ${params.starsolo_whitelist}" : "--soloCBwhitelist None"
    def cb_len    = params.cb_len  ?: 16
    def umi_len   = params.umi_len ?: 12
    """
    STAR \\
        --soloType CB_UMI_Simple \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads_2} ${reads_1} \\
        ${whitelist} \\
        --soloCBstart 1 --soloCBlen ${cb_len} \\
        --soloUMIstart 17 --soloUMIlen ${umi_len} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --soloFeatures GeneFull \\
        --outFileNamePrefix ${meta.id}/ \\
        --runThreadN ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """
}
