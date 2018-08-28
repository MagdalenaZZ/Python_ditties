#!/bin/bash

cd /scratch/DCT/BIOINFCO/jamesc/pipelines/end-seq-pipline

module load anaconda/3/latest

#python end-seq-pipeline.py \
#--sample_file /scratch/DCT/BIOINFCO/jamesc/pipelines/end-seq-pipline/samples.txt \
#--output_dir /scratch/DCT/BIOINFCO/jamesc/pipelines/end-seq-pipline/testing \
#--genome_file /scratch/readonly/bwa/07x/Mus_musculus.GRCm38.dna.chromosome.combined.fa.gz

python end-seq-pipeline.py \
--sample_file /scratch/DCT/BIOINFCO/jamesc/ibc0002/canela_2016_mol_cell_fastq/samples.txt \
--output_dir /scratch/DCT/BIOINFCO/jamesc/pipelines/end-seq-pipline/testing \
--genome_file /scratch/readonly/bwa/07x/Mus_musculus.GRCm38.dna.chromosome.combined.fa.gz



