process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    
    publishDir "${params.outdir}/quantification/10x/${meta.id}", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path fasta
    path gtf
    
    output:
    tuple val(meta), path("${meta.id}/outs/filtered_feature_bc_matrix"), emit: counts
    tuple val(meta), path("${meta.id}/outs/web_summary.html"), emit: web_summary
    tuple val(meta), path("${meta.id}/outs/possorted_genome_bam.bam"), emit: bam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def chemistry = params.chemistry ?: 'auto'
    def expect_cells = params.expect_cells ?: 5000
    def transcriptome = params.cellranger_transcriptome ?: 'transcriptome'
    
    """
    # Create Cell Ranger reference if not exists
    if [ ! -d "${transcriptome}" ]; then
        cellranger mkref \\
            --genome=reference \\
            --fasta=${fasta} \\
            --genes=${gtf} \\
            --nthreads=${task.cpus}
        transcriptome="reference"
    fi
    
    # Run Cell Ranger count
    cellranger count \\
        --id=${meta.id} \\
        --transcriptome=\${transcriptome} \\
        --fastqs=. \\
        --sample=${meta.id} \\
        --chemistry=${chemistry} \\
        --expect-cells=${expect_cells} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version | sed 's/cellranger cellranger-//')
    END_VERSIONS
    """
}
