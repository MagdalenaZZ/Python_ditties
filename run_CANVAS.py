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
	Make the preparation for CANVAS, test tumour and normal sample\n\
     \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running deTIN tumour-in-normal subtraction', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-t', '--tumour_bam', default=None, dest='tum', action='store', required=True, help="Tumour.bam file full path")
parser.add_argument('-v', '--vcf', default=None, dest='vcf', action='store', required=True, help="Somatic VCF file")
parser.add_argument('-n', '--nvcf', default=None, dest='nvcf', action='store', required=True, help="Normal variants VCF file")
parser.add_argument('-o', '--out_folder', default=None, dest='out', action='store', required=True, help="Output folder full path")


gc_prof = '/home/mzarowiecki/bin/hartwigmedicalfoundation/GC_profile.hg38.1000bp.cnp'
jar_loc = '/home/mzarowiecki/bin/hartwigmedicalfoundation/hmftools_pipeline_v4_3'
circos_loc = '/home/mzarowiecki/bin/circos-0.69-6/bin/circos'
# Needs a BED file with human variants
#human_bed = '/home/mzarowiecki/scratch/REF/af-only-gnomad.hg38.01.chr.bed'
human_bed = '/home/mzarowiecki/scratch/REF/GermlineHetPon.hg38.bed'

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()


# Check if input files exist
if not os.path.isfile(args.tum)==True:
    print("Cannot find input file ",args.tum)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)

# Check if NVCF file exits 
if not os.path.isfile(args.nvcf)==True:
    print("Cannot find input file ",args.nvcf)
    sys.exit(1)

# Output will be created if it doesnt exist

# Create prefixes and files
tum=args.tum.rstrip()
tumpx = tum.split("/")[-1]
tumpx = ''.join(tumpx.split(".bam")[0:-1])

sh=args.out + '/' + tumpx + '.sh'
f = open(sh, 'w')

# Print bsub header
print ('#!/bin/bash \n\
\n\
# a template to run canvas-spw on trio samples\n\
# when multisample SNV vcfs already created for B allele frequencies\n\
\n\
# the filename for STDOUT\n\
#BSUB -o %s/%s.o\n\
#BSUB -e %s/%s.e\n\
\n\
# The queue to which the job is to be submitted\n\
#BSUB -q bio\n\
# project code for analysis, research and development activities\n\
#BSUB -P Analysis\n\
\n\
# Memory usage configuration\n\
#BSUB -R "span[hosts=1]" \n\
#BSUB -R "select[largedsk]" \n\
#BSUB -R "hname!=\'hpc-prod-grid-lsfexec-001\' && hname!=\'hpc-prod-grid-lsfexec-002\'" \n\
#BSUB -n 4\n\
\n\
# add canvas application profile so that jobs are not preemtied too many times (see https://jira.extge.co.uk/browse/INFRA-6931)\n\
#BSUB -app canvas\n\
' % (args.out,tumpx,args.out,tumpx), file=f )

# Print file locations

print ('\n\
source /etc/profile.d/modules.sh\n\
\n\
module load canvas/1.39.0.1598\n\
module load bcftools/1.5\n\
\n\
TODAY=`date +%%Y-%%m-%%d`\n\
\n\
# refence files\n\
\n\
GENOME=/home/mzarowiecki/scratch/REF/reference_no_alt\n\
KMER=/home/mzarowiecki/scratch/REF/reference_no_alt/kmer.GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa\n\
#KMER=/genomes/scratch/dkasperaviciute/sv_cnv/canvas-spw/edico/reference_no_alt/kmer.GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa\n\
#GENOME=/genomes/scratch/dkasperaviciute/sv_cnv/canvas-spw/edico/reference_no_alt\n\
#KMER=/genomes/scratch/dkasperaviciute/sv_cnv/canvas-spw/edico/reference/kmer.GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fa\n\
#GENOME=/genomes/scratch/dkasperaviciute/sv_cnv/canvas-spw/edico/reference\n\
FILTER=/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/hg38-NSv6/Annotation/Canvas/filter13.bed\n\
PLO=/home/mzarowiecki/scratch/Benchmark_CNV_callers/Canvas/Canvas139/male_ploidy.vcf.gz\n\
\n\
#input files\n\
TBAM=%s\n\
VCF=%s\n\
NVCF=%s\n\
NAME=%s\n\
\n\
# output folder\n\
CANVAS_OUTPUT_DIR=%s\n\
mkdir -p $CANVAS_OUTPUT_DIR\n\
cd $CANVAS_OUTPUT_DIR\n\
\n\
' % (args.tum, args.vcf, args.nvcf, tumpx, args.out), file=f)



# Germline free

print ('canvas Somatic-WGS -b $TBAM  --somatic-vcf=$VCF --sample-b-allele-vcf=$NVCF -n $NAME -o $CANVAS_OUTPUT_DIR/%s -r $KMER -g $GENOME -f $FILTER --ploidy-vcf=$PLO' % (tumpx), file=f)

f.close()

print ('bsub -J %s < %s' % (tumpx, sh) ) 
            

quit()








