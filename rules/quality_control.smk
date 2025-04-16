rule fastqc:
    input: 
        os.path.join(RAWFASTQ_DIR, "{sample}_{direction}.fastq.gz")
    output: 
        os.path.join(FASTQC_RAW_DIR, "{sample}_{direction}_fastqc.html")
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
        os.path.join(FASTQC_PROCESSED_DIR, "{filename}_fastqc.html")
    wildcard_constraints:
        direction = "[1-2]"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{filename}.out",
        err = "log/fastqc_{filename}.err"
    shell: 
        "fastqc {input} -t {threads} -o $(dirname {output}) 2> {log.err} 1> {log.out}"