process STAR_GENOME_BUILD {
    label 'process_high'

    // Cache the genome index so it's built once per run
    storeDir "${params.outdir}/genome_index/star"

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    def sa_bases = params.genome_sa_bases ?: 14
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --genomeSAindexNbases ${sa_bases} \\
        --runThreadN ${task.cpus}
    """
}
