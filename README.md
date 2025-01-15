Snakemake Workflow for scRNA-seq Analysis
This workflow processes single-cell RNA-seq (scRNA-seq) data obtained from an Illumina sequencer and performs differential expression analysis. It is built using Snakemake and leverages popular bioinformatics tools like FastQC, Trim Galore!, HISAT2, featureCounts, and R-based analysis tools such as DESeq2 via Python's rpy2 library.
________________________________________
Table of Contents
1.	Workflow Overview
2.	Requirements
3.	Installation
4.	Configuration
5.	Running the Workflow
6.	Outputs
________________________________________
Workflow Overview
This pipeline performs the following steps:
1.	Quality Control (FastQC):
o	Generates quality control reports for raw FASTQ files.
2.	Trimming (Trim Galore!):
o	Removes adapters and low-quality bases from raw reads.
3.	Alignment (HISAT2):
o	Aligns trimmed reads to the reference genome to produce sorted BAM files.
4.	Gene Count Generation (featureCounts):
o	Counts reads mapped to each gene using a GTF annotation file.
5.	Aggregate Counts:
o	Combines gene count files from multiple samples into a single count matrix.
6.	Differential Expression Analysis (DESeq2):
o	Identifies differentially expressed genes between specified conditions using DESeq2 via Python's rpy2 interface.
________________________________________
Requirements
Software Dependencies
•	Snakemake (v7.0+)
•	Python (v3.8+) 
o	Required Python packages: 
pandas
rpy2
•	R (v4.0+) 
o	Required R packages: 
DESeq2
•	Bioinformatics Tools: 
o	FastQC
o	Trim Galore!
o	HISAT2
o	featureCounts (part of Subread package)
Installation of Bioinformatics Tools
For details on installing these tools, refer to their respective documentation or use the following commands:
Linux/MacOS
# Install Snakemake
conda install -c bioconda -c conda-forge snakemake

# Install FastQC
sudo apt-get install fastqc

# Install Trim Galore!
sudo apt-get install trim-galore

# Install HISAT2
sudo apt-get install hisat2

# Install Subread (featureCounts)
sudo apt-get install subread
Python Packages
pip install pandas rpy2
R Packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
________________________________________
Configuration
Edit the config.yaml file to define the input files, reference genome, and other parameters. Below is an example configuration:
# Path to genome index for HISAT2
genome_index: "/path/to/genome/index"

# Path to annotation file (GTF format)
annotation_gtf: "/path/to/annotation.gtf"

# List of sample names (without file extensions)
samples:
  - sample1
  - sample2
  - sample3
________________________________________
Running the Workflow
Step 1: Set Up Your Project Directory
Place your raw FASTQ files in a directory (e.g., data/raw) and ensure the file names match the sample names defined in config.yaml.
Step 2: Execute the Workflow
Run the following command in the directory containing the Snakefile:
snakemake --cores <number_of_cores>
For example, to run on 8 cores:
snakemake --cores 8
Step 3: Dry-Run (Optional)
To test the pipeline without running it:
snakemake --dry-run
________________________________________
Outputs
The pipeline generates the following key outputs:
1.	Quality Control Reports:
o	results/fastqc/: Contains FastQC reports (*_fastqc.html and *_fastqc.zip).
2.	Trimmed Reads:
o	results/trimmed/: Contains adapter-trimmed FASTQ files (*_trimmed.fq.gz).
3.	Aligned Reads:
o	results/aligned/: Contains sorted BAM files (*.bam) and their indices (*.bam.bai).
4.	Gene Counts:
o	results/counts/: Contains gene count files for each sample (*_counts.txt).
5.	Differential Expression Results:
o	results/differential_expression_results.csv: Contains the results of the differential expression analysis, including log2 fold changes and adjusted p-values.

________________________________________
Notes
•	Ensure all required tools and dependencies are installed and accessible in your PATH.
•	If using a cluster, Snakemake can be executed with cluster configuration to distribute tasks across nodes.
•	Adjust parameters in config.yaml and scripts as needed for your dataset and experiment design.
________________________________________
For questions or issues, please contact the workflow maintainer or refer to the official documentation of the tools used.

