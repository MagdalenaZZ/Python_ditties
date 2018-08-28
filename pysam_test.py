#!/usr/bin/env python
from __future__ import division
import pysam
import sys
import os.path
import argparse

"""

Script for parsing and editing a VCF file

"""


epi = ('\
    \n\
	VCF file parser, allowing for custom advanced filtering of VCF files\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and applies custom filtering', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF file")
#parser.add_argument('-f', '--filter-string', default=None, dest='fi', action='store', required=True, help="String of filters to apply")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)

# Read input
output=args.vcf+".pysam.vcf"
print("Input: ",args.vcf, "Output: " , output)


# read the input file
myvcf = pysam.VariantFile(args.vcf, "r")
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)

for r in myvcf:
    vcf_out.write(r)






