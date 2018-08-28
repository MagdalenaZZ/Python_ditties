#!/usr/bin/python 

import argparse
import os
import ngs

"""
run END-seq jobs
================
Read a file with fastq paths and sample line
names. Start LSF jobs for each pair of
files.
"""



# get file names from the command line args
parser = argparse.ArgumentParser()
parser.add_argument('--sample_file', help='location of a file containing sample names and fastq files to align')
parser.add_argument('--output_dir', help='location of the output directory root')
parser.add_argument('--genome_file', help='location of a file containing the reference genome, indexed for BWA')

args = parser.parse_args()
sample_file = args.sample_file
output_dir = args.output_dir
genome_file = args.genome_file

samples = ngs.read_samples(sample_file)

for sample in samples.keys():
	fastq1 = samples[sample].pop()
	#fastq2 = samples[sample].pop()
	ngs.ensure_dir("%s/%s/" % (output_dir, sample))
	
	'''
	commands needs to contain the job submission script for LSF
	Steps are:
		1. load all modules
		2. run bwa mem to produce a sam file
		3. run samtools view to generate a bam
		4. run samtools sort to sort and index (Picard?)
		5. run bedtools merge to find all overlapping regions
		6. run coverageBed to find depth of reads in each overlap
		7. may need some other post-processing to tidy up tables
		8. remove intermediate files and generally clean up
	'''
	
	commands = '''\
#!/bin/bash
#BSUB -J %(sample)s
#BSUB -o %(output_dir)s/%(sample)s/endseq.out
#BSUB -e %(output_dir)s/%(sample)s/endseq.err
#BSUB -n 8
#BSUB -q normal
#BSUB -P DROFXHABM
#BSUB -W 168:00
#BSUB -R "span[ptile=8]"

module load bwa/0.7.12 bedtools/2.25.0 samtools/1.3

cd %(output_dir)s/%(sample)s

bwa mem \\
%(genome_file)s \\
%(fastq1)s \\
| gzip -3 \\
> %(sample)s.sam.gz \\
&& samtools view \\
-b \\
%(sample)s.sam.gz \\
> %(sample)s.bam \\
&& samtools sort \\
%(sample)s.bam \\
> %(sample)s.sorted.bam \\
&& samtools index \\
%(sample)s.sorted.bam \\
&& bamToBed \\
-i %(sample)s.sorted.bam \\
> %(sample)s.sorted.bed \\
&& mergeBed \\
-i %(sample)s.sorted.bed \\
> %(sample)s.sorted.merged.bed \\
&& coverageBed \\
-counts \\
-a %(sample)s.sorted.merged.bed \\
-b %(sample)s.sorted.bam \\
> %(sample)s.sorted.merged.coverage.txt 

''' % locals()
	output = open("%(output_dir)s/%(sample)s/run_endseq.sh" % locals(), 'w')
	output.write(commands)
	output.close
	
	bsubcmd = "bsub < %(output_dir)s/%(sample)s/run_endseq.sh" % locals()
	#return_code = os.system( bsubcmd )
	return_code = os.popen(bsubcmd).readlines()
	print("return code from bsub: %s" % (return_code))

