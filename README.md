# Pipeline to process PRO-seq FASTQ to mapped reads
A custom Snakemake pipeline to process PRO-seq results. Currently limited to paired-end sequencing. 

Steps involved include:
1. Download FASTQ from NCBI SRA database if necessary
2. Trim adapters
3. UMI extraction 
4. Download genomes(s) and build bowtie2 indices
5. Map reads using Bowtie 2
6. Extract individual BAM files for each genome 

## Setting up project

### Sample Metadata
This pipeline uses PEP files to define sample metadata. PEP involves two files: PEP.yaml and PEP.csv. See [pep.databio.org](https://pep.databio.org/) for more information. The PEP should contain 4 pieces of info for each sample 

1. Sample name
2. SRR (if downloading from SRA database)
3. Genomes used for each sample (can be defined in PEP.yaml)
4. Genome GenBank accessions (can be defined in PEP.yaml)
 
## Dependencies
You can make apptainer image that contains all dependencies. If you're not going to use the apptainer image, then recommend running snakemake with `--use-conda` flag and that you need to have `pd` and `peppy` python packages installed. 

### How to make apptainer image from container.def file
The container.def file was made as follows, which requires snakemake and spython to generate. It takes all of conda environment files to build a container definition file

```
snakemake --containerize > Dockerfile
spython recipe Dockerfile container.def 
```

To then create an apptainer image, we use the following: 

```
apptainer build container.sif container.def 
```

## TO DO
1. Make UMI extraction optional
2. Allow for both PE and SE reads
3. Make it so that files don't have to be manually named
4. Create a config file for settings