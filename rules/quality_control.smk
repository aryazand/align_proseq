rule fastqc:
    input: 
        os.path.join(RAWFASTQ_DIR, "{sample}_{direction}.fastq.gz")
    output: 
        os.path.join(QC_DIR_FASTQC, "{sample}_{direction}_fastqc.html")
    wildcard_constraints:
        direction = "[1-2]"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}.out",
        err = "log/fastqc_{sample}_{direction}.err"
    shell: 
        "fastqc {input} -t {threads} -o $(dirname {output}) 2> {log.err} 1> {log.out}"

rule fastqc_on_processed_fastq:
    input: 
        os.path.join(TRIMMED_DIR, "{filename}.fastq.gz"),
    output: 
        os.path.join(QC_DIR_FASTQC, "{filename}_fastqc.html")
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{filename}.out",
        err = "log/fastqc_{filename}.err"
    shell: 
        "fastqc {input} -t {threads} -o $(dirname {output}) 2> {log.err} 1> {log.out}"

rule multiqc_fastqc:
    input:
        expand(os.path.join(QC_DIR_FASTQC, "{sample}_{direction}_fastqc.html"), sample = sample_table['sample_name'], direction = ["1", "2"]),
        expand(os.path.join(QC_DIR_FASTQC, "{sample}_{direction}_trimmed_umi_fastqc.html"), sample = sample_table['sample_name'], direction = ["1", "2"]),
    output:
        os.path.join(MULTIQC_DIR, "fastqc/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc_fastqc.out",
        err = "log/multiqc_fastqc.err"
    params:
        QC_DIR = QC_DIR_FASTQC
    shell:
        "multiqc {params.QC_DIR} --outdir $(dirname {output}) --force 2> {log.err} 1> {log.out}"

rule multiqc_trimming:
    input:
        expand(os.path.join(QC_DIR_TRIMMING, "{sample}_{direction}_trimming_report.txt"), sample = sample_table['sample_name'], direction = ["1", "2"]),
    output:
        os.path.join(MULTIQC_DIR,  "trimming/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc_trimming.out",
        err = "log/multiqc_trimming.err"
    params:
        QC_DIR = QC_DIR_TRIMMING
    shell:
        "multiqc {params.QC_DIR} --outdir $(dirname {output}) --force 2> {log.err} 1> {log.out}"

rule multiqc_umi_extraction:
    input:
        expand(os.path.join(QC_DIR_UMIEXTRACT, "{sample}_umi_report.txt"), sample = sample_table['sample_name']),
    output:
        os.path.join(MULTIQC_DIR, "umi_extraction/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc_umi_extraction.out",
        err = "log/multiqc_umi_extraction.err"
    params:
        QC_DIR = QC_DIR_UMIEXTRACT
    shell:
        "multiqc {params.QC_DIR} --outdir $(dirname {output}) --force 2> {log.err} 1> {log.out}"

rule multiqc_alignment:     
    input:
        expand(os.path.join(QC_DIR_ALIGNMENT, "{sample}.stats"), sample = sample_table['sample_name']),
    output:
        os.path.join(MULTIQC_DIR, "alignment/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc_alignment.out",
        err = "log/multiqc_alignment.err"
    params:
        QC_DIR = QC_DIR_ALIGNMENT
    shell:
        "multiqc {params.QC_DIR} --outdir $(dirname {output}) --force 2> {log.err} 1> {log.out}"

rule multiqc_dedup:
    input:
        expand(os.path.join(QC_DIR_DEDUP, "{sample}.out"), sample = sample_table['sample_name']),
    output:
        os.path.join(MULTIQC_DIR, "dedup/multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc_dedup.out",
        err = "log/multiqc_dedup.err"
    params:
        QC_DIR = QC_DIR_DEDUP
    shell:
        "multiqc {params.QC_DIR} --outdir $(dirname {output}) --force 2> {log.err} 1> {log.out}"

rule multiqc:
    input:
        os.path.join(MULTIQC_DIR, "fastqc/multiqc_report.html"),
        os.path.join(MULTIQC_DIR, "trimming/multiqc_report.html"),
        os.path.join(MULTIQC_DIR, "umi_extraction/multiqc_report.html"),
        os.path.join(MULTIQC_DIR, "alignment/multiqc_report.html"),
        os.path.join(MULTIQC_DIR, "dedup/multiqc_report.html")