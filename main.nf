#!/usr/bin/env nextflow

/*
========================================================================================
    scRNA-SEQ Analysis Pipeline
========================================================================================
    Github : https://github.com/sakyoney/scRNA-SEQ-Flow
    Author : scRNA-SEQ-Flow Team
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
if (!params.input) {
    exit 1, "Input samplesheet not specified! Please provide --input"
}

if (!params.protocol || !(params.protocol in ['10x', 'smartseq2'])) {
    exit 1, "Protocol must be specified as '10x' or 'smartseq2'. Current: ${params.protocol}"
}

if (!params.genome) {
    exit 1, "Genome must be specified! Use --genome GRCh38, GRCm39, or custom"
}

// Check genome configuration
if (params.genome == 'custom') {
    if (!params.fasta || !params.gtf) {
        exit 1, "Custom genome requires --fasta and --gtf parameters"
    }
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTQC                  } from './modules/local/fastqc'
include { MULTIQC                 } from './modules/local/multiqc'
include { CELLRANGER_COUNT        } from './modules/local/cellranger'
include { STARSOLO                } from './modules/local/starsolo'
include { STAR_ALIGN              } from './modules/local/star_smartseq'
include { FEATURECOUNTS           } from './modules/local/featurecounts'
include { SCANPY_MERGE            } from './modules/local/scanpy_merge'
include { SCANPY_QC               } from './modules/local/scanpy_qc'
include { SCANPY_PREPROCESS       } from './modules/local/scanpy_preprocess'
include { SCANPY_INTEGRATE        } from './modules/local/scanpy_integrate'
include { SCANPY_CLUSTER          } from './modules/local/scanpy_cluster'
include { SCANPY_MARKERS          } from './modules/local/scanpy_markers'
include { SCANPY_DIFFEXP          } from './modules/local/scanpy_diffexp'
include { SCANPY_ANNOTATE         } from './modules/local/scanpy_annotate'
include { GENERATE_REPORT         } from './modules/local/generate_report'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow SCRNA_SEQ_FLOW {

    // Parse input samplesheet
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id
            meta.condition = row.condition
            meta.protocol = row.protocol ?: params.protocol

            def fastq_1 = file(row.fastq_r1, checkIfExists: true)
            def fastq_2 = row.fastq_r2 ? file(row.fastq_r2, checkIfExists: true) : null

            return [meta, fastq_1, fastq_2]
        }

    // Get genome reference files
    ch_fasta = params.genome == 'custom' ?
        Channel.fromPath(params.fasta) :
        Channel.fromPath(params.genomes[params.genome].fasta)

    ch_gtf = params.genome == 'custom' ?
        Channel.fromPath(params.gtf) :
        Channel.fromPath(params.genomes[params.genome].gtf)

    // FASTQC - Quality control on raw reads
    FASTQC(ch_input)

    // Branch by protocol
    ch_input
        .branch {
            tenx: it[0].protocol == '10x'
            smartseq2: it[0].protocol == 'smartseq2'
        }
        .set { ch_branched }

    // 10x Genomics Processing
    if (params.use_cellranger) {
        CELLRANGER_COUNT(
            ch_branched.tenx,
            ch_fasta,
            ch_gtf
        )
        ch_10x_counts = CELLRANGER_COUNT.out.counts
    } else {
        STARSOLO(
            ch_branched.tenx,
            ch_fasta,
            ch_gtf
        )
        ch_10x_counts = STARSOLO.out.counts
    }

    // Smart-seq2 Processing
    STAR_ALIGN(
        ch_branched.smartseq2,
        ch_fasta,
        ch_gtf
    )

    FEATURECOUNTS(
        STAR_ALIGN.out.bam,
        ch_gtf
    )

    ch_smartseq2_counts = FEATURECOUNTS.out.counts

    // Merge all quantification results
    ch_all_counts = ch_10x_counts
        .mix(ch_smartseq2_counts)
        .collect()

    SCANPY_MERGE(ch_all_counts)

    // Quality Control and Filtering
    SCANPY_QC(SCANPY_MERGE.out.h5ad)

    // Normalization and Preprocessing
    SCANPY_PREPROCESS(SCANPY_QC.out.h5ad)

    // Optional: Batch Integration
    if (params.integration_method != 'none') {
        SCANPY_INTEGRATE(SCANPY_PREPROCESS.out.h5ad)
        ch_processed = SCANPY_INTEGRATE.out.h5ad
    } else {
        ch_processed = SCANPY_PREPROCESS.out.h5ad
    }

    // Clustering Analysis
    SCANPY_CLUSTER(ch_processed)

    // Find Marker Genes per Cluster
    SCANPY_MARKERS(SCANPY_CLUSTER.out.h5ad)

    // Differential Expression: Control vs Disease
    SCANPY_DIFFEXP(SCANPY_CLUSTER.out.h5ad)

    // Cell Type Annotation
    if (params.auto_annotate) {
        SCANPY_ANNOTATE(SCANPY_CLUSTER.out.h5ad)
        ch_final = SCANPY_ANNOTATE.out.h5ad
    } else {
        ch_final = SCANPY_CLUSTER.out.h5ad
    }

    // MultiQC - Aggregate all QC reports
    ch_multiqc_files = Channel.empty()
        .mix(FASTQC.out.zip.collect().ifEmpty([]))

    MULTIQC(ch_multiqc_files)

    // Generate Final Report
    GENERATE_REPORT(
        ch_final,
        SCANPY_MARKERS.out.markers,
        SCANPY_DIFFEXP.out.degs,
        MULTIQC.out.report
    )
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    SCRNA_SEQ_FLOW()
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    if (workflow.success) {
        log.info """
        ============================================
        Pipeline completed successfully!
        ============================================
        Results are in: ${params.outdir}
        """.stripIndent()
    } else {
        log.info """
        ============================================
        Pipeline completed with errors!
        ============================================
        """.stripIndent()
    }
}
