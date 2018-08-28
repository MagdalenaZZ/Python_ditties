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


# Add the HP field to header. Say its a string and can take any values. It
# depends what format you want to give.
myvcf.header.formats.add("FQ", ".", "String", "Frequency of alternative allele")
myvcf.header.formats.add("GT", ".", "String", "Genotype")
myvcf.header.formats.add("QU", ".", "String", "Quality")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)
#with open(output, "a") as out:


# Add the new HP format to all samples
for r in myvcf:

    # Only put things that have PASS
    #if (ref.qual !='PASS'):
    #    continue

    refb = r.ref
    altb = r.alts[0]
    refd = 0
    altd = 0

    #print(refb, altb)

    # If is SNP
    if 'CU' in r.samples[0].keys():

        if (r.ref == 'A'):
            refd = r.samples[0]['AU'][0]
        elif (r.ref == 'C'):
            refd = r.samples[0]['CU'][0]
        elif (r.ref == 'G'):
            refd = r.samples[0]['GU'][0]
        elif (r.ref == 'T'):
            refd = r.samples[0]['TU'][0]
        else:
            print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

        if (r.alts[0] == 'A'):
            altd = r.samples[0]['AU'][0]
        elif (r.alts[0] == 'C'):
            altd = r.samples[0]['CU'][0]
        elif (r.alts[0] == 'G'):
            altd = r.samples[0]['GU'][0]
        elif (r.alts[0] == 'T'):
            altd = r.samples[0]['TU'][0]
        else:
            print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

        if ( (altd + refd) < 1):
            print ("WARN: This SNP call is not supported by reads", str(r) )
            continue

        freq = altd / (altd + refd)

        # Insert genotype
        gtref = 0
        gtalt = 0

        if ( 0.1 < freq  ):
            gtalt=1
        gt = (gtref, gtalt)
        #print (gt)

        r.samples[0]['FQ'] = str(freq)
        r.samples[0]['GT'] = gt
        r.samples[0]['QU'] = str(r.filter.items()[0][0])
        r.samples[0]['DP'] = (altd + refd)
        #print ("Rqual", str(r.filter.items()[0][0]))


        #if str(r).startswith( 'chr' ):
        vcf_out.write(r)
        #else:
        #    print("WARN: is not valid input", str(r))

    # If is indel
    elif 'TIR' in r.samples[0].keys():

        # print(refb, refd, altb, altd)
        # Insert frequency
        # vaf = (TIR / float(TAR + TIR))
        #print ( r.samples[0]['TAR'][0],r.samples[0]['TIR'][0] )
        if ((r.samples[0]['TAR'][0]+r.samples[0]['TIR'][0]) < 1):
            print ("WARN: This INDEL call is not supported by reads", str(r) )
            #r.samples[0]['TAR'][0]=[1, 1]
            continue
        if ((r.samples[0]['TAR'][0]+r.samples[0]['TIR'][0]) > r.samples[0]['DP']):
            print ("WARN: This call DP is smaller than exp ", str(r) )
            #continue

        freq = r.samples[0]['TIR'][0] / (r.samples[0]['TAR'][0]+r.samples[0]['TIR'][0])

        # Insert genotype
        gtref = 0
        gtalt = 0

        if (0.1 < freq):
            gtalt=1
        gt = (gtref, gtalt)

        r.samples[0]['FQ'] = str(freq)
        r.samples[0]['GT'] = gt
        r.samples[0]['QU'] = str(r.filter.items()[0][0])
        r.samples[0]['DP'] = (r.samples[0]['TAR'][0]+r.samples[0]['TIR'][0])
        #print ("Rqual", str(r.filter.items()[0][0]))

        vcf_out.write(r)

    else:
        print ("Warn: line format not as expected")




exit(0)



"""

        hp_value = 'Hola'
        for sample in r.samples:
            #for val in variant.samples[sample]['DP']:
            hp_value += str(r.pos) + '-'
            if r.samples[sample]['DP'] >= 100:
                hp_value = "larger"
                r.samples[sample]['HP'] = hp_value
                # print hp_value
            if r.samples[sample]['DP'] < 100:
                hp_value = "smaller"
                r.samples[sample]['HP'] = hp_value
        out.write(str(r))

"""