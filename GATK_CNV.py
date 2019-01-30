#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import re
from subprocess import call

"""

Script for parsing a VCF file

"""


epi = ('\
    \n\
	Run GATK CNV pipeline from tumour and normal sample\n\
    Your genome.intervals_list is output from GATK pipeline, and a panel-of-normals (from PON.py) \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running GATK CNV pipeline from tumour and normal data', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-g', '--genome', default=None, dest='genome', action='store', required=True, help="genome.interval_list")
parser.add_argument('-d', '--genome.dict', default=None, dest='gdic', action='store', required=True, help="genome.dict file")
parser.add_argument('-n', '--normal_bam', default=None, dest='norm', action='store', required=True, help="Normal.bam file")
parser.add_argument('-t', '--tumour_bam', default=None, dest='tum', action='store', required=True, help="Tumour.bam file")
parser.add_argument('-p', '--pon', default=None, dest='pon', action='store', required=True, help="Panel of Normal CNV hdf5 file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.genome)==True:
    print("Cannot find input file ",args.genome)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.norm)==True:
    print("Cannot find input file ",args.norm)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.tum)==True:
    print("Cannot find input file ",args.tum)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.gdic)==True:
    print("Cannot find input file ",args.gdic)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.pon)==True:
    print("Cannot find input file ",args.pon)
    sys.exit(1)

# Create prefixes and files
norm=args.norm.rstrip()
normpx = norm.split("/")[-1]
normpx = ''.join(normpx.split(".bam")[0:-1])
tum=args.tum.rstrip()
tumpx = tum.split("/")[-1]
tumpx = ''.join(tumpx.split(".bam")[0:-1])


# Create versions of the genome, and link to here
genome=args.genome.split("/")[-1]
genomepx = '.'.join(genome.split(".")[0:-1])
px=tumpx


# Create hdf5 files for tumour and normal
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (normpx,normpx,normpx,norm,args.genome,normpx))
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (tumpx,tumpx,tumpx,tum,args.genome,tumpx))


# Denoise the read counts for all tumour and normal samples, using the panel of normals
print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--count-panel-of-normals %s.pon.hdf5  --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv' % (normpx,normpx,normpx,normpx,px,normpx,normpx))

print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--count-panel-of-normals %s.pon.hdf5  --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv' % (tumpx,tumpx,tumpx,tumpx,px,tumpx,tumpx))


# Use customised R
print ('source activate r-3.5.0\nmkdir Plots')


# Plot the two files showed above to highlight differences and sanity check
print ('gatk PlotDenoisedCopyRatio  --standardized-copy-ratios %s.standardizedCR.tsv  --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (normpx,normpx,genomepx,normpx))
print ('gatk PlotDenoisedCopyRatio  --standardized-copy-ratios %s.standardizedCR.tsv  --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (tumpx,tumpx,genomepx,tumpx))


# ModelSegments create segments across the genome for paired tumour and normal sample
print ("mkdir Segments")
print ('gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios %s.denoisedCR.tsv \
--allelic-counts %s.allelicCounts.tsv \
--normal-allelic-counts %s.allelicCounts.tsv \
--output Segments \
--output-prefix %s ' % (tumpx,tumpx,normpx,tumpx) )



# Call copy-neutral, amplified and deleted segments with CallCopyRatioSegments

print ('gatk CallCopyRatioSegments \
--input Segments/%s.cr.seg \
--output Segments/%s.called.seg ' % (tumpx,tumpx) )



# PlotModeledSegments visualizes copy and allelic ratio segmentation results.

print ('gatk PlotModeledSegments \
--denoised-copy-ratios %s.denoisedCR.tsv \
--allelic-counts %s.hets.tsv \
--segments %s.modelFinal.seg \
--sequence-dictionary %s \
--minimum-contig-length 46709983 \
--output Plots \
--output-prefix %s ' % (tumpx,tumpx,tumpx,args.gdic,tumpx) )


# Check outputs



quit()