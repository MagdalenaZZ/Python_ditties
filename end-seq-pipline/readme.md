END-seq Pipeline
================

A pipeline for alignment and reporting of END-seq data

Writes a job submission script for LSF that

1. load all modules (need to set versions
2. run bwa mem to produce a sam file
3. run samtools view to generate a bam
4. run samtools sort to sort and index (Picard?)
5. run bedtools merge to find all overlapping regions
6. run coverageBed to find depth of reads in each overlap
7. may need some other post-processing to tidy up tables
8. remove intermediate files and generally clean up




