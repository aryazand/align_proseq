rule download_genome:
    # download genome from NCBI
    output: 
        fna = "data/genome/{species}/{accession}.fna",
        gff = "data/genome/{species}/{accession}.gff"
    conda:
        "../envs/get-genome.yml"
    log:
        "log/download_genome_{species}_{accession}.log",
    shell:
        """
        datasets download genome accession {wildcards.accession} --filename {wildcards.species}.zip --include gff3,genome
        unzip {wildcards.species}.zip -d {wildcards.species}
        mv {wildcards.species}/ncbi_dataset/data/{wildcards.accession}/*.fna {output.fna}
        mv {wildcards.species}/ncbi_dataset/data/{wildcards.accession}/*.gff {output.gff}
        sed -i -re 's/(>\\S*)\\s.*/\\1/' {output.fna}
        rm {wildcards.species}.zip
        rm -r {wildcards.species}
        """

rule get_chrom_sizes:
    # get chromosome sizes
    input:
        "data/genome/{species}/{genome}.fna"
    output:
        "data/genome/{species}/{genome}.chrom.sizes"
    conda:
        "../envs/get-genome.yml"
    log:
        "log/get_chrom_sizes_{species}_{genome}.log"
    shell:
        """
        bioawk -cfastx '{{ print $name, length($seq) }}' {input} > {output}
        """

rule concatenate_genomes:
    input:
        expand("data/genome/{species}/{genome}.fna", zip, species = GENOMES.keys(), genome = GENOMES.values())
    output:
        "data/genome/{combined_species_names}.fna"
    log:
        out = "log/concatenate_genomes_{combined_species_names}.out",
        err = "log/concatenate_genomes_{combined_species_names}.err"
    shell:
        """
        cat {input} | sed '/^>/ s/[[:space:]]/\\_/g' > {output}
        """        