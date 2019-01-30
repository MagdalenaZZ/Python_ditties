#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os.path
import argparse
import csv
import pysam
import re
from subprocess import call
from pathlib import Path


"""

Script for converting a strelka2 VCF to a formatted somatic VCF

"""


epi = ('\
    \n\
	File parser,  VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a Strelka2 VCF file and formats it', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, nargs = '+' , help="VCF.gz file(s)")
parser.add_argument('-b', '--bed', default=None, dest='bed', action='store', required=False, help="BED boundary file, put DEF for default Tx")
#parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF.gz file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)


args  = parser.parse_args()

# Check if input files exist and create an index, if the index does not exist

#print (args.vcf)


"""
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)
if not (os.path.isfile(args.vcf+".tbi")==True or os.path.isfile(args.vcf+".csi")==True ):
    call(["bcftools","index",args.vcf])
"""

bed='/genomes/scratch/magz/REF/SureSelectV5_TRACERx_Edition.padded.reduced.hg38.bed'


if args.bed is None:
    print ("BED is undefined")
elif args.bed is not None:
    if re.match( r'DEF', args.bed):
        print ("BED is default : %s" % (bed))
    else:
        print ("BED is custom file :  %s " % (args.bed))
        if not os.path.isfile(args.bed)==True:
            print("Cannot find input file ",args.bed)
            sys.exit(1)



prefix= args.vcf[0].split('/')[-1]
px=prefix.split('.')[0]


# create an object of new bed file and open in to write data.
output=px+".doTx.sh"
out = open(output, 'w')

#print (prefix,px)

vcfs= ' '.join(args.vcf)

# Collate input files, if there are several
print ("bcftools concat -a -D -O z -o %s.somatic.vcf.gz  %s" % (px, vcfs),file=out)
print ("bcftools index  %s.somatic.vcf.gz" % (px),file=out)



# Filter by BED, if BED exists
if args.bed is None:
    print ("ln -s %s.somatic.vcf.gz %s.filter.vcf.gz" % (px,px),file=out)
else:
    #bcftools view -R /genomes/scratch/magz/REF/SureSelectV5_TRACERx_Edition.padded.reduced.hg38.bed -O z -o ../3_VCFs_2/LN3999999-DNA_A01_LP3000396-DNA_A01.TxPASS.vcf.gz ../3_VCFs_2/LN3999999-DNA_A01_LP3000396-DNA_A01.somatic.vcf.gz
    print ("bcftools view -R %s -O z -o %s.filter.vcf.gz %s.somatic.vcf.gz" % (bed,px,px),file=out)



# Keep only the tumour
# bcftools view -s LP3000396-DNA_E01 -O z -o LN3999999-DNA_E01_LP3000396-DNA_E01.Txfil.vcf.gz LN3999999-DNA_E01_LP3000396-DNA_E01.TxPASS.vcf.gz
print ("bcftools view -s TUMOR -O z -o %s.filt.to.vcf.gz %s.filter.vcf.gz" % (px,px ),file=out)



# Reformat to general VCF format - to be compatile with Tx
# python ~/git/TRACERx_validation/VCF_our_reform_to_general.py -i LN3999999-DNA_A01_LP3000396-DNA_A01.Txfil.vcf.gz
print ("python ~/git/TRACERx_validation/VCF_our_reform_to_general.py -i %s.filt.to.vcf.gz" % (px),file=out)



# Remove redundant fields
# bcftools  annotate -x "INFO/ALTMAP,INFO/ALTPOS,INFO/DP,INFO/IC,INFO/IHP,INFO/MQ,INFO/MQ0,INFO/NT,INFO/OVERLAP,INFO/PNOISE,INFO/PNOISE2,INFO/QSI,INFO/QSI_NT,INFO/QSS_NT,INFO/RC,INFO/ReadPosRankSum,INFO/RU,INFO/SGT,INFO/SNVSB,INFO/TQSI,INFO/TQSI_NT,INFO/TQSS,INFO/TQSS_NT,INFO/AF1000G,INFO/AA,INFO/GMAF,INFO/cosmic,INFO/clinvar,INFO/EVS,INFO/RefMinor,INFO/phyloP,INFO/CSQT,INFO/CSQR,INFO/OLD_MULTIALLELIC,INFO/OLD_VARIANT,FORMAT/AU,FORMAT/CU,FORMAT/GU,FORMAT/TU,FORMAT/FDP,FORMAT/SDP,FORMAT/SUBDP,FORMAT/DP2,FORMAT/TAR,FORMAT/TIR,FORMAT/TOR,FORMAT/DP50,FORMAT/FDP50,FORMAT/SUBDP50" LN3999999-DNA_A01_LP3000396-DNA_A01.Txfil.vcf.gz.pysam.vcf.gz | bgzip > LN3999999-DNA_A01_LP3000396-DNA_A01.TIN.vcf.gz
print ("bcftools  annotate -x \"INFO/DP,INFO/IC,INFO/IHP,INFO/MQ,INFO/MQ0,INFO/NT,INFO/OVERLAP,INFO/PNOISE,INFO/PNOISE2,INFO/QSI,INFO/QSI_NT,INFO/QSS_NT,INFO/RC,INFO/ReadPosRankSum,INFO/RU,INFO/SGT,INFO/SNVSB,INFO/TQSI,INFO/TQSI_NT,INFO/TQSS,INFO/TQSS_NT,FORMAT/AU,FORMAT/CU,FORMAT/GU,FORMAT/TU,FORMAT/FDP,FORMAT/SDP,FORMAT/SUBDP,FORMAT/DP2,FORMAT/TAR,FORMAT/TIR,FORMAT/TOR,FORMAT/DP50,FORMAT/FDP50,FORMAT/SUBDP50\" %s.filt.to.vcf.gz.pysam.vcf.gz | bgzip > %s.pruned.vcf.gz"  % (px, px),file=out)
print ("vt sort -o %s.compTx.vcf.gz %s.pruned.vcf.gz" % (px,px),file=out)
print ("bcftools index  %s.compTx.vcf.gz" % (px),file=out)
print ("rm %s.filter.vcf.gz %s.filt.to.vcf.gz %s.filt.to.vcf.gz.pysam.vcf.gz %s.pruned.vcf.gz" % (px,px,px,px),file=out)


out.close()


print ("bash %s" % (output))


exit(0)


