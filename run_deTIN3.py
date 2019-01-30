#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import re
from subprocess import call

"""

Script for preparing and running deTIN

"""


epi = ('\
    \n\
	Make the preparation for deTIN, test tumour and normal sample\n\
     \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running deTIN tumour-in-normal subtraction', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-v', '--vcf', default=None, dest='vcf', action='store', required=True, help="Somatic VCF file")
parser.add_argument('-f', '--gfasta', default=None, dest='gfa', action='store', required=True, help="genome.fa")
parser.add_argument('-g', '--genome_intervals', default=None, dest='genome', action='store', required=True, help="genome.interval_list")
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
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)
"""
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

"""
# Create prefixes and files
norm=args.norm.rstrip()
normpx = norm.split("/")[-1]
normpx = ''.join(normpx.split(".bam")[0:-1])
tum=args.tum.rstrip()
tumpx = tum.split("/")[-1]
tumpx = ''.join(tumpx.split(".bam")[0:-1])


# Create versions of the genome, and link to here
#gfa=args.genome
genome=args.gfa.split("/")[-1]
genomepx = '.'.join(genome.split(".")[0:-1])


# Make a fake Mutect1 file from a VCF
print('source activate py27.12; python ~/bin/Python_ditties/convert_VCF_to_Mutect1.py -i %s' % (args.vcf))


#Collect counts at germline variant sites for the matched-control
######- there can only be one sample in the VCF file produced
###### so take the someatic VCF, split out so only one sample is in there
###### index the sample with 
###### gatk IndexFeatureFile -F /genomes/scratch/magz/TIN_on_Tx/4_doTIN/LN3999999-DNA_E01.tum.vcf.gz
###### then run that

print ('bsub -P Analysis -J CACn%s -o CACn%s.o -e CACn%s.e -n 1 -R "rusage[mem=30000]" -q bio \
gatk --java-options "-Xmx29g" CollectAllelicCounts  -L %s \
 -I %s -R %s -O %s.norm.allelicCounts.tsv' % (normpx,normpx,normpx,args.vcf,args.norm,args.gfa,normpx))

#Collect counts at the same sites for the case sample
print ('bsub -P Analysis -J CACt%s -o CACt%s.o -e CACt%s.e -n 1 -R "rusage[mem=30000]" -q bio \
gatk --java-options "-Xmx29g" CollectAllelicCounts -L %s \
-I %s -R %s -O %s.tum.allelicCounts.tsv' % (tumpx,tumpx,tumpx,args.vcf,args.tum,args.gfa,tumpx) )


## Or now do for normal
#Collect counts at germline variant sites for the matched-control
# print ('bsub -P Analysis -J CACng%s -o CACng%s.o -e CACng%s.e -n 1 -R "rusage[mem=30000]" -q bio \
# gatk --java-options "-Xmx29g" CollectAllelicCounts  -L %s \
# -I %s.vcf.gz -R %s -O %s.normg.allelicCounts.tsv' % (normpx,normpx,normpx,normpx,args.norm,args.gfa,normpx))

#Collect counts at the same sites for the case sample
#print ('bsub -P Analysis -J CACtg%s -o CACtg%s.o -e CACtg%s.e -n 1 -R "rusage[mem=30000]" -q bio \
#gatk --java-options "-Xmx29g" CollectAllelicCounts -L %s \
#-I %s.vcf.gz -R %s -O %s.tumg.allelicCounts.tsv' % (tumpx,tumpx,tumpx,normpx,args.tum,args.gfa,tumpx) )


# Create hdf5 files for tumour and normal
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (normpx,normpx,normpx,norm,args.genome,normpx))
print ('bsub -P Analysis -J CRC%s -o CRC%s.o -e CRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio  gatk --java-options "-Xmx9g" \
CollectReadCounts -I %s -L %s --interval-merging-rule OVERLAPPING_ONLY -O %s.hdf5' % (tumpx,tumpx,tumpx,tum,args.genome,tumpx))


# Denoise the read counts for all tumour and normal samples, using the panel of normals
print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--count-panel-of-normals %s  --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv' % (normpx,normpx,normpx,normpx,args.pon,normpx,normpx))

print('bsub -P Analysis -J DRC%s -o DRC%s.o -e DRC%s.e  -n 1  -R "rusage[mem=9500]" -q bio gatk --java-options "-Xmx9g" DenoiseReadCounts -I %s.hdf5 \
--count-panel-of-normals %s  --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv' % (tumpx,tumpx,tumpx,tumpx,args.pon,tumpx,tumpx))

# Use customised R
print ('source activate r-3.5.0\nmkdir Plots')


# Plot the two files showed above to highlight differences and sanity check
print ('bsub -P Analysis -J PCR%s -o PCR%s.o -e PCR%s.e  -n 1 -q bio \
gatk PlotDenoisedCopyRatio --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (normpx,normpx,normpx,normpx,normpx,genomepx,normpx))

print ('bsub -P Analysis -J PCR%s -o PCR%s.o -e PCR%s.e  -n 1 -q bio \
gatk PlotDenoisedCopyRatio --standardized-copy-ratios %s.standardizedCR.tsv --denoised-copy-ratios %s.denoisedCR.tsv \
--sequence-dictionary %s.dict --minimum-contig-length 46709983 --output Plots --output-prefix %s' % (tumpx,tumpx,tumpx,tumpx,tumpx,genomepx,tumpx))


# ModelSegments create segments across the genome for paired tumour and normal sample
print ("mkdir Segments")
print ('bsub -P Analysis -J MS%s -o MS%s.o -e MS%s.e  -n 1  -R "rusage[mem=6500]" -q bio gatk --java-options "-Xmx6g" ModelSegments \
--denoised-copy-ratios %s.denoisedCR.tsv \
--allelic-counts %s.allelicCounts.tsv \
--normal-allelic-counts %s.allelicCounts.tsv \
--output Segments \
--output-prefix %s ' % (tumpx,tumpx,tumpx,tumpx,tumpx,normpx,tumpx) )



# Remove the headers
cat LP3000396-DNA_E01.tum.allelicCounts.tsv.44.tsv | grep -v '@' > 

# Run deTIN after all files have been created and fixed
print ('PYTHONPATH=/home/mzarowiecki/bin/Python_package/lib/python2.7/site-packages:/home/mzarowiecki/bin/deTiN/deTiN')

print ('python /home/mzarowiecki/bin/deTiN/deTiN/deTiN.py  \
--mutation_data_path LP3001098-DNA_A04_LP3001095-DNA_B02.somatic.vcf.gz.mutect1.stats  \
--cn_data_path LP3001095-DNA_B02.modelFinal.seg \
--tumor_het_data LP3001095-DNA_B02.hets.tsv \
--normal_het_data LP3001095-DNA_B02.hets.normal.tsv \
--exac_data_path ExaC.pickle \
--output_name test \
--output_dir example_data' % () )


quit()



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
