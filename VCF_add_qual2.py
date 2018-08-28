#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import pysam
import sys
import os.path
import argparse
import subprocess


"""

Script for parsing and editing a VCF file

"""

epi = ('\
    \n\
	VCF file converter, moves QUAL to FORMAT\n\
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
output=args.vcf+".qadd.vcf"
output=output.replace('.vcf.gz','')
#outc= args.vcf+".qadd.vcf.gz"
print("Input: ",args.vcf, "Output: ")


# read the input file
myvcf = pysam.VariantFile(args.vcf, "r")


# Add the QU field to header. Say its a string and can take any values.
myvcf.header.formats.add("QU", ".", "String", "Quality")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)

# Add the new HP format to all samples
for r in myvcf:

    #print (r.samples)

    if 'DP' in r.samples[0].keys():
        r.samples[0]['QU'] = r.filter.keys()
        r.samples[1]['QU'] = r.filter.keys()
        vcf_out.write(r)



#print("bcftools", "view","-O","z","-o",outc,output)
subprocess.call(["bgzip",output])

exit(0)
