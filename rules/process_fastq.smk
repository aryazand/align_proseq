rule trim_reads_pe:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz"
    output:
        r1=temp(os.path.join(TRIMMED_DIR, "{sample}_1_trimmed.fastq.gz")),
        r2=temp(os.path.join(TRIMMED_DIR, "{sample}_2_trimmed.fastq.gz")),
        report_r1=os.path.join(QC_DIR_TRIMMING, "{sample}_1_trimming_report.txt"),
        report_r2=os.path.join(QC_DIR_TRIMMING, "{sample}_2_trimming_report.txt")
    conda:
        "../envs/trim-galore.yml"
    log:
        out = "log/trim_reads_{sample}.out",
        err = "log/trim_reads_{sample}.err"
    params:
        cores = config["trim_galore"]["threads"],
        adaptor = config["trim_galore"]["adaptor"],
        stringency = config["trim_galore"]["stringency"],
        additional = config["trim_galore"]["additional"]
    shell:
	    """
        # make directory if it does not exist
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.report_r1})

        # run trim galore
        trim_galore --paired --cores {params.cores} --gzip \
            {params.adaptor} --stringency {params.stringency} \
            {params.additional} \
            --output_dir $(dirname {output.r1}) \
            {input.r1} {input.r2} \
            2> {log.err} 1> {log.out}
        
        # rename files
        mv $(dirname {output.r1})/{wildcards.sample}_1_val_1.fq.gz {output.r1}
        mv $(dirname {output.r2})/{wildcards.sample}_2_val_2.fq.gz {output.r2}

        # move trimming reports
        mv $(dirname {output.r1})/{wildcards.sample}_1.fastq.gz_trimming_report.txt {output.report_r1}
        mv $(dirname {output.r2})/{wildcards.sample}_2.fastq.gz_trimming_report.txt {output.report_r2}
        """

rule extract_umi: 
    input: 
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_trimmed.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_trimmed.fastq.gz")
    output:
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_trimmed_umi.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_trimmed_umi.fastq.gz"),
        report = os.path.join(QC_DIR_UMIEXTRACT, "{sample}_umi_report.txt")
    log:
        err = "log/extract_umi_{sample}.err"
    conda: 
        "../envs/umitools.yml"
    params: 
        bc_pattern = config["umi_tools_extract"]["bc_pattern"],
        additional = config["umi_tools_extract"]["additional"]
    shell:
        """
        umi_tools extract -I {input.r1} --read2-in={input.r2} \
        --bc-pattern={params.bc_pattern} --bc-pattern2={params.bc_pattern} \
        {params.additional} \
        --stdout={output.r1} --read2-out={output.r2} \
        2> {log.err} 1> {output.report}
        """