#########################
# Define Environment
#########################

import os
import peppy
import pandas as pd

containerized: "container.sif"

#############################
# Define directories
#############################

DATA_DIR = "data"
RAWFASTQ_DIR = os.path.join(DATA_DIR, "fastq")
GENOMES_DIR = os.path.join(DATA_DIR, "genomes")
RESULTS_DIR = "results"
TRIMMED_DIR = os.path.join(RESULTS_DIR, "trimmed")
DEDUPPED_DIR = os.path.join(RESULTS_DIR, "dedupped")
ALIGNMENT_DIR = os.path.join(RESULTS_DIR, "alignments")
LOGS_DIR = os.path.join(RESULTS_DIR, "logs")

QC_DIR = os.path.join(RESULTS_DIR, "qc")
FASTQC_RAW_DIR = os.path.join(QC_DIR, "fastqc_raw")
FASTQC_PROCESSED_DIR = os.path.join(QC_DIR, "fastqc_processed")
QC_DIR_TRIMMING = os.path.join(QC_DIR, "trimming_reports")
QC_DIR_ALIGNMENT = os.path.join(QC_DIR, "alignment_reports")
QC_DIR_DEDUP = os.path.join(QC_DIR, "dedup")
QC_DIR_UMIEXTRACT = os.path.join(QC_DIR, "umi_extraction")
MULTIQC_DIR = os.path.join(QC_DIR, "multiqc")

#############################
# Load samples and metadata
#############################

# # Load project PEP
# project = peppy.Project("data/sample_metadata/PEP.yaml")

# # Get sample metadata
# sample_table = project.sample_table

# # Create a genomes key
# dict_list = []
# for i in range(len(sample_table)):
#     genome_names = sample_table['genome_names'].iloc[i].split(", ")
#     genome_accessions = sample_table['genome_accessions'].iloc[i].split(", ")
#     dict_list.append(dict(zip(genome_names, genome_accessions)))

# for i in range(len(dict_list)-1):
#     dict_list[i].update(dict_list[i+1])

# GENOMES = dict_list[0]

######################
# Define output files
#####################

rule all:
    input:
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam"), sample = sample_table['sample_name']),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam.bai"), sample = sample_table['sample_name']),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_{genome}.bam"), sample = sample_table['sample_name'], genome = GENOMES.keys()),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_{genome}.bam.bai"), sample = sample_table['sample_name'], genome = GENOMES.keys())

rule multiqc:
    output:
        os.path.join(QC_DIR, "multiqc/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc.out",
        err = "log/multiqc.err"
    params:
        QC_DIR = QC_DIR
    shell:
        "multiqc {params.QC_DIR} --outdir ($dirname {output}) --force"

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
include: "rules/get_genomic_data.smk"