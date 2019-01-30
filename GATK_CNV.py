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
parser.add_argument('-f', '--genome.fa', default=None, dest='gfa', action='store', required=True, help="genome.fa file")
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

print ('gatk AnnotateIntervals -L 2_Segments/genome.preprocessed_intervals.interval_list -R genome.fa --sequence-dictionary genome.dict --interval-merging-rule OVERLAPPING_ONLY -O genome.GC.AnnotateIntervals')

# Create hdf5 files for tumour and normal
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (normpx,normpx,normpx,norm,args.genome,normpx))
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (tumpx,tumpx,tumpx,tum,args.genome,tumpx))


# Denoise the read counts for all tumour and normal samples, using the panel of normals
print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--annotated-intervals genome.GC.AnnotateIntervals --count-panel-of-normals %s  --standardized-copy-ratios %s.standardizedCR.tsv \
--denoised-copy-ratios %s.denoisedCR.tsv' % (normpx,normpx,normpx,normpx,args.pon,normpx,normpx))

print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--annotated-intervals genome.GC.AnnotateIntervals --count-panel-of-normals %s  --standardized-copy-ratios %s.standardizedCR.tsv \
--denoised-copy-ratios %s.denoisedCR.tsv' % (tumpx,tumpx,tumpx,tumpx,args.pon,tumpx,tumpx))


# Use customised R
print ('source activate r-3.5.0\nmkdir Plots')


# Plot the two files showed above to highlight differences and sanity check
print ('bsub -P Analysis -J PCR%s -o PCR%s.o -e PCR%s.e  -n 1 -q bio \
gatk PlotDenoisedCopyRatio --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (normpx,normpx,normpx,normpx,normpx,genomepx,normpx))

print ('bsub -P Analysis -J PCR%s -o PCR%s.o -e PCR%s.e  -n 1 -q bio \
gatk PlotDenoisedCopyRatio --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (tumpx,tumpx,tumpx,tumpx,tumpx,genomepx,tumpx))



#Collect counts at germline variant sites for the matched-control

print ('bsub -P Analysis -J CAC%s -o CAC%s.o -e CAC%s.e -n 1 -R "rusage[mem=30000]" -q bio \
gatk --java-options "-Xmx29g" CollectAllelicCounts  -L /home/mzarowiecki/scratch/REF/af-only-gnomad.hg38.01.chr.vcf.gz \
 -I %s -R %s -O %s.allelicCounts.tsv' % (normpx,normpx,normpx,args.norm,args.gfa,normpx))

#Collect counts at the same sites for the case sample

print ('bsub -P Analysis -J CAC%s -o CAC%s.o -e CAC%s.e -n 1 -R "rusage[mem=30000]" -q bio \
gatk --java-options "-Xmx29g" CollectAllelicCounts -L /home/mzarowiecki/scratch/REF/af-only-gnomad.hg38.01.chr.vcf.gz \
-I %s -R %s -O %s.allelicCounts.tsv' % (tumpx,tumpx,tumpx,args.tum,args.gfa,tumpx) )


# ModelSegments create segments across the genome for paired tumour and normal sample
print ("mkdir Segments")
print ('bsub -P Analysis -J MS%s -o MS%s.o -e MS%s.e  -n 1  -R "rusage[mem=30000]" -q bio gatk --java-options "-Xmx29g" ModelSegments \
--denoised-copy-ratios %s.denoisedCR.tsv \
--allelic-counts %s.allelicCounts.tsv \
--normal-allelic-counts %s.allelicCounts.tsv \
--output Segments \
--output-prefix %s ' % (tumpx,tumpx,tumpx,tumpx,tumpx,normpx,tumpx) )



# Call copy-neutral, amplified and deleted segments with CallCopyRatioSegments

print ('bsub -P Analysis -J CRS%s -o CRS%s.o -e CRS%s.e  -n 1 -q bio gatk CallCopyRatioSegments \
--input Segments/%s.cr.seg \
--output Segments/%s.called.seg ' % (tumpx,tumpx,tumpx,tumpx,tumpx) )


print ('source activate r-3.5.0\nmkdir SegPlots')

# PlotModeledSegments visualizes copy and allelic ratio segmentation results.

print ('bsub -P Analysis -J PMS%s -o PMS%s.o -e PMS%s.e  -n 1 -q bio gatk PlotModeledSegments \
--denoised-copy-ratios %s.denoisedCR.tsv \
--allelic-counts Segments/%s.hets.tsv \
--segments Segments/%s.modelFinal.seg \
--sequence-dictionary %s \
--minimum-contig-length 46709983 \
--output SegPlots \
--output-prefix %s ' % (tumpx,tumpx,tumpx,tumpx,tumpx,tumpx,args.gdic,tumpx) )


# Tag germline events

print '(gatk TagGermlineEvents)'


# Check outputs



##################
#
# Continue with the GATK ACNV workflow
# https://gatkforums.broadinstitute.org/gatk/discussion/7387/description-and-examples-of-the-steps-in-the-acnv-case-workflow
#
##################


print ('gatk GetBayesianHetCoverage --reference %s --snpIntervals /home/mzarowiecki/scratch/REF/af-only-gnomad.hg38.01.chr.vcf.gz --tumor %s --tumorHets %s.tum.hets.tsv --normal %s --normalHets %s.norm.hets.tsv --hetCallingStringency 30' % (args.gfa,args.tum,tumpx,args.norm,normpx))


print ('gatk AllelicCNV  --tumorHets %s.tum.hets.tsv --tangentNormalized <coverage_profile> --segments <called_segments> --outputPrefix <output_prefix>' % (tumpx))


print ('gatk CallCNLoHAndSplits  --tumorHets <tumor_het_pulldown>
    --segments <acnv_segments> --tangentNormalized <coverage_profile> --outputDir <output_dir>
    --rhoThreshold 0.2 --numIterations 10  ')



quit()
