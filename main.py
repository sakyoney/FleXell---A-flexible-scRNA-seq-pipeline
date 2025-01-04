import os
from snakemake.shell import shell

# Define paths to input and output directories
configfile: "config.yaml"

rule all:
    input:
        "results/differential_expression_results.csv"

rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "results/fastqc/{sample}_fastqc.html",
        "results/fastqc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} -o results/fastqc/"

rule trim_galore:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "results/trimmed/{sample}_trimmed.fq.gz"
    shell:
        "trim_galore {input} -o results/trimmed/"

rule align:
    input:
        "results/trimmed/{sample}_trimmed.fq.gz"
    output:
        "results/aligned/{sample}.bam"
    params:
        index = config["genome_index"]
    shell:
        "hisat2 -x {params.index} -U {input} | samtools sort -o {output}"

rule index_bam:
    input:
        "results/aligned/{sample}.bam"
    output:
        "results/aligned/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule count_genes:
    input:
        bam = "results/aligned/{sample}.bam",
        bai = "results/aligned/{sample}.bam.bai"
    output:
        "results/counts/{sample}_counts.txt"
    params:
        gtf = config["annotation_gtf"]
    shell:
        "featureCounts -a {params.gtf} -o {output} {input.bam}"

rule aggregate_counts:
    input:
        expand("results/counts/{{sample}}_counts.txt", sample=config["samples"])
    output:
        "results/counts/all_samples_counts.txt"
    shell:
        "python scripts/aggregate_counts.py {input} > {output}"

rule differential_expression:
    input:
        "results/counts/all_samples_counts.txt"
    output:
        "results/differential_expression_results.csv"
    script:
        "scripts/differential_expression_analysis.py"

