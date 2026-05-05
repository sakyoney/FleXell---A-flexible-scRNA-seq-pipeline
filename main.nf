#!/usr/bin/env nextflow

/*
========================================================================================
    FleXell — Flexible scRNA-seq Pipeline
========================================================================================
    Supports: 10x Genomics (STARsolo default / Cell Ranger optional)
              Smart-seq2 (plate-based, one FASTQ per cell)
    GitHub : https://github.com/sakyoney/FleXell---A-flexible-scRNA-seq-pipeline
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.input) {
    exit 1, "Input samplesheet not specified! Please provide --input"
}

if (!params.protocol || !(params.protocol in ['10x', 'smartseq2', 'mixed'])) {
    exit 1, "Protocol must be '10x', 'smartseq2', or 'mixed'. Current: ${params.protocol}"
}

if (!params.genome) {
    exit 1, "Genome must be specified! Use --genome GRCh38, GRCm39, or custom"
}

if (params.genome == 'custom' && (!params.fasta || !params.gtf)) {
    exit 1, "Custom genome requires --fasta and --gtf"
}

log.info """
============================================
  F L E X e l l  scRNA-seq Pipeline
============================================
  input       : ${params.input}
  protocol    : ${params.protocol}
  genome      : ${params.genome}
  outdir      : ${params.outdir}
  aligner     : ${params.use_cellranger ? 'cellranger' : 'starsolo'}
  integration : ${params.integration_method}
  annotation  : ${params.auto_annotate ? params.celltypist_model : 'disabled'}
============================================
""".stripIndent()

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTQC            } from './modules/local/fastqc'
include { MULTIQC           } from './modules/local/multiqc'
include { STAR_GENOME_BUILD } from './modules/local/star_genome'
include { CELLRANGER_COUNT  } from './modules/local/cellranger'
include { STARSOLO          } from './modules/local/starsolo'
include { STAR_ALIGN        } from './modules/local/star_smartseq'
include { FEATURECOUNTS     } from './modules/local/featurecounts'
include { SCANPY_MERGE      } from './modules/local/scanpy_merge'
include { SCANPY_QC         } from './modules/local/scanpy_qc'
include { SCANPY_PREPROCESS } from './modules/local/scanpy_preprocess'
include { SCANPY_INTEGRATE  } from './modules/local/scanpy_integrate'
include { SCANPY_CLUSTER    } from './modules/local/scanpy_cluster'
include { SCANPY_MARKERS    } from './modules/local/scanpy_markers'
include { SCANPY_DIFFEXP    } from './modules/local/scanpy_diffexp'
include { SCANPY_ANNOTATE   } from './modules/local/scanpy_annotate'
include { GENERATE_REPORT   } from './modules/local/generate_report'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow FLEXELL {

    // --- Parse samplesheet ---
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta      = [id: row.sample_id, condition: row.condition,
                             protocol: row.protocol ?: params.protocol]
            def fastq_1   = file(row.fastq_r1, checkIfExists: true)
            def fastq_2   = row.fastq_r2 ? file(row.fastq_r2, checkIfExists: true) : []
            return [meta, fastq_1, fastq_2]
        }

    // --- Reference genome ---
    ch_fasta = params.genome == 'custom'
        ? Channel.fromPath(params.fasta, checkIfExists: true)
        : Channel.fromPath(params.genomes[params.genome].fasta)

    ch_gtf = params.genome == 'custom'
        ? Channel.fromPath(params.gtf, checkIfExists: true)
        : Channel.fromPath(params.genomes[params.genome].gtf)

    // --- FASTQC ---
    FASTQC(ch_input)

    // --- Build STAR genome index once ---
    STAR_GENOME_BUILD(ch_fasta, ch_gtf)

    // --- Branch by protocol ---
    ch_input.branch {
        tenx:      it[0].protocol == '10x'
        smartseq2: it[0].protocol == 'smartseq2'
    }.set { ch_branched }

    // --- 10x quantification ---
    if (params.use_cellranger) {
        CELLRANGER_COUNT(ch_branched.tenx, ch_fasta, ch_gtf)
        ch_10x_counts = CELLRANGER_COUNT.out.counts
    } else {
        STARSOLO(ch_branched.tenx, STAR_GENOME_BUILD.out.index)
        ch_10x_counts = STARSOLO.out.counts
    }

    // --- Smart-seq2 quantification (one FASTQ per cell) ---
    STAR_ALIGN(ch_branched.smartseq2, STAR_GENOME_BUILD.out.index)
    FEATURECOUNTS(STAR_ALIGN.out.bam, ch_gtf)
    ch_smartseq2_counts = FEATURECOUNTS.out.counts

    // --- Merge all quantification results ---
    ch_all_counts = ch_10x_counts
        .mix(ch_smartseq2_counts)
        .map { meta, counts -> counts }
        .collect()

    ch_samplesheet = Channel.fromPath(params.input)

    SCANPY_MERGE(ch_all_counts, ch_samplesheet)

    // --- QC & Filtering ---
    SCANPY_QC(SCANPY_MERGE.out.h5ad)

    // --- Normalization & PCA ---
    SCANPY_PREPROCESS(SCANPY_QC.out.h5ad)

    // --- Batch Integration ---
    if (params.integration_method != 'none') {
        SCANPY_INTEGRATE(SCANPY_PREPROCESS.out.h5ad)
        ch_processed = SCANPY_INTEGRATE.out.h5ad
    } else {
        ch_processed = SCANPY_PREPROCESS.out.h5ad
    }

    // --- Clustering ---
    SCANPY_CLUSTER(ch_processed)

    // --- Marker genes & DE ---
    SCANPY_MARKERS(SCANPY_CLUSTER.out.h5ad)
    SCANPY_DIFFEXP(SCANPY_CLUSTER.out.h5ad)

    // --- Cell type annotation ---
    if (params.auto_annotate) {
        SCANPY_ANNOTATE(SCANPY_CLUSTER.out.h5ad)
        ch_final = SCANPY_ANNOTATE.out.h5ad
    } else {
        ch_final = SCANPY_CLUSTER.out.h5ad
    }

    // --- MultiQC ---
    ch_multiqc_files = FASTQC.out.zip.map { meta, zip -> zip }.collect()
    MULTIQC(ch_multiqc_files)

    // --- Final Report ---
    GENERATE_REPORT(
        ch_final,
        SCANPY_MARKERS.out.markers,
        SCANPY_DIFFEXP.out.degs,
        MULTIQC.out.report
    )
}

/*
========================================================================================
    ENTRY POINT
========================================================================================
*/

workflow {
    FLEXELL()
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
        ============================================
        FleXell completed successfully!
        Results : ${params.outdir}
        Report  : ${params.outdir}/report/flexell_report.html
        ============================================
        """.stripIndent()
    } else {
        log.error "FleXell completed with errors — check logs above."
    }
}
