#!/bin/sh
# Bioinformatics for Biologists: Analysing and Interpreting Genomics Datasets

# Install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda list
export PATH=~/miniconda3/bin:$PATH

# Import the conda environment from provided yml file
unzip MOOC.yml.zip
conda env create -n MOOC --file MOOC.yml
source activate MOOC

# Download the sequence data from the European Nucleotide Archive using the accession number
fastq-dump --split-files ERR5743893

# File organisation
mkdir data
mv ERR5743893_* data/

# Quality control
mkdir qc_reports
fastqc data/ERR5743893_1.fastq data/ERR5743893_2.fastq --outdir qc_reports/
cd qc_reports/
multiqc .

# Index the reference genome
bwa index ref_gen/MN908947.fasta

# Create a fasta.fai file
## Required for visualisation of the .bam file we will later generate in the IGV
samtools faidx ref_gen/MN908947.fasta

# Map the sequences from the target sample to the reference genome
bwa mem ref_gen/MN908947.fasta data/ERR5743893_1.fastq data/ERR5743893_2.fastq > mapping/ERR5743893.sam
cd mapping/
ls -lhrt

# Convert .sam file to a .bam file
samtools view -@ 20 -S -b mapping/ERR5743893.sam > mapping/ERR5743893.bam

# Sorting the file to have the reads ordered based on where they align in the reference genome, instead of keeping the default order of reads from the fastq files
samtools sort -@ 10 -o mapping/ERR5743893.sorted.bam mapping/ERR5743893.bam

# Index the sorted .bam file
samtools index mapping/ERR5743893.sorted.bam
## The .bam file can now be visualised on IGV.

# Variant Calling using FreeBayes
freebayes -f ref_gen/MN908947.fasta mapping/ERR5743893.sorted.bam > mapping/ERR5743893.vcf

## Zip and index the .vcf file
bgzip mapping/ERR5743893.vcf
tabix mapping/ERR5743893.vcf.gz

## Create and view the stats .vcf file
bcftools stats ERR5743893.vcf.gz > ERR5743893_stat.vcf.txt

# Install Nextflow
java -version
curl -fsSL get.nextflow.io | bash
## Conda channels have already been set, so we will create an environment called nextflow and install Nextflow in it right away
conda create --name nextflow nextflow

## Activate nextflow
conda activate nextflow
nextflow help

# Install Singularity
sudo apt install -y runc cryptsetup-bin
wget -O singularity.deb https://github.com/sylabs/singularity/releases/download/v3.11.4/singularity-ce_3.11.4-jammy_amd64.deb
sudo dpkg -i singularity.deb
rm singularity.deb

## ENA accessions to analyse
echo "ERR5556343

SRR13500958

ERR5743893

ERR5181310

ERR5405022" > samples.txt

# Run fastq-dump on each of the ENA accessions
for i in $(cat samples.txt);do fastq-dump --split-files $i;done

## Compress and organise the fastq files
gzip *.fastq
mkdir data/viralrecon
mv *.fastq.gz data/viralrecon

# Create the sample sheet with the sample names and location of the fastq files
## nf-core provide a python script that generates the sample sheet automatically
wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py
python3 fastq_dir_to_samplesheet.py data/viralrecon samplesheet.csv -r1 _1.fastq.gz -r2 _2.fastq.gz

# Run viralrecon
conda activate nextflow
nextflow run nf-core/viralrecon -profile singularity \
--max_memory '12.GB' --max_cpus 4 \
--input samplesheet.csv \
--outdir results/viralrecon \
--protocol amplicon \
--genome 'MN908947.3' \
--primer_set artic \
--primer_set_version 3 \
--skip_kraken2 \
--skip_assembly \
--skip_pangolin \
--skip_nextclade \
--platform illumina

# Comparing FreeBayes and viralrecon variant calling results
## FreeBayes identified 76 variants in the ERR5743893 accession, including 75 SNPs and 1 indel, whereas viralrecon identified 49 variants, including 27 SNPs, 1 indel and 21 missense variants.
