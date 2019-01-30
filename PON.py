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
	Create PON from genome and list of normals\n\
    Your genome needs to have a .fai and .dict file available too \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script creates a GATK PON from genome and normal data', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-g', '--genome', default=None, dest='genome', action='store', required=True, help="genome.fa")
parser.add_argument('-l', '--list', default=None, dest='list', action='store', required=True, help="List of normal.bam files")
parser.add_argument('-p', '--prefix', default="test", dest='prefix', action='store', required=False, help="Output prefix")


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
if not os.path.isfile(args.list)==True:
    print("Cannot find input file ",args.list)
    sys.exit(1)

list = open(args.list, "r")
list2 = open(args.list, "r")
list3 = open(args.list, "r")

# Create versions of the genome, and link to here
genome=args.genome.split("/")[-1]
px=args.prefix
#print(genome,args.genome,"REF")

if ((str(genome)==str(args.genome))==False):
    #call(["ln", "-s",args.genome])
    print ("ln -s %s\nln -s %s.fai" % (args.genome,args.genome))

genomepx = '.'.join(genome.split(".")[0:-1])


# Create genome dict and fai if it doesnt exist
print('bsub -q bio -P Analysis -J CSD%s -o CSD%s.o -e CSD%s.e gatk -n 1 --java-options "-Xmx2g" CreateSequenceDictionary -R %s -O %s.dict ' % (px,px,px,genome,genomepx))


"""
# Create intervals from genome # FAST, easy, less than 1 minute
print ('gatk PreprocessIntervals  -R %s  --bin-length 1000  --padding 0  -O %s.interval_list' % (genome,genome ))


# Preprocess intervals; exomes
#print ("  %s " % (args.genome))
print ('gatk PreprocessIntervals \
-L %s.interval_list \
-R %s \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
-O %s.preprocessed.interval_list \
' % (genome, genome, genome) )

"""

# To generate consecutive bins of 1000 bases from the reference (e.g., for whole genome sequencing analyses):
print ('gatk PreprocessIntervals \
 -R %s \
 --bin-length 1000 \
 --padding 0 \
 -O %s.interval_list' % (genome,genomepx) )

# Filter the result to only include major chromosomes
print ('cat  %s.interval_list | grep chr | grep -v chrEBV | grep -v chrUn | grep -v random | grep -v ChrM | grep -v ChrX | grep -v ChrY  > %s.preprocessed_intervals.interval_list' %(genomepx,genomepx))


hdf5f=[]
# Make a small fake panel-of-normals side-by-side file # Takes long 7519 sec, is big 10G, preferably paralellise
for norm in list:
    norm=norm.rstrip()
    normpx = norm.split("/")[-1]
    normpx = ''.join(normpx.split(".bam")[0:-1])
    print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=5500]" -q bio  gatk --java-options "-Xmx5g" CollectReadCounts -I %s -L %s.preprocessed_intervals.interval_list --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (normpx,normpx,normpx,norm,genomepx,normpx))
    hdf5f.append("%s.hdf5" %(normpx))



# Create a PON intervals file from multiple samples input
print ('bsub -P Analysis -J PON%s -o PON%s.o -e PON%s.e  -n 1  -R "rusage[mem=10000]" -q bio gatk --java-options "-Xmx9g" \
CreateReadCountPanelOfNormals --minimum-interval-median-percentile 5.0 -O %s.pon.hdf5 ' %(px,px,px,px), end=" " )
for file in hdf5f:
    file = file.rstrip()
    print('\ \n -I  %s   ' % (file), end=" ")


# Denoise the read counts for all normal samples, using the panel of normals

for norm in list2:
    norm=norm.rstrip()
    normpx = norm.split("/")[-1]
    normpx = ''.join(normpx.split(".bam")[0:-1])
    print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=8500]" -q bio \
gatk --java-options "-Xmx8g" DenoiseReadCounts -I %s.hdf5 \
--count-panel-of-normals %s.pon.hdf5  --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv' % (normpx,normpx,normpx,normpx,px,normpx,normpx))



# Use customised R
print ('\nsource activate r-3.5.0\n')
print("\nmkdir Plots%s\n" %(px))

# Plot the two files showed above to highlight differences
for norm in list3:
    norm=norm.rstrip()
    normpx = norm.split("/")[-1]
    normpx = ''.join(normpx.split(".bam")[0:-1])
    print ('gatk PlotDenoisedCopyRatios --standardized-copy-ratios %s.standardizedCR.tsv  --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots%s --output-prefix %s' % (normpx,normpx,genomepx,px,normpx))



quit()
