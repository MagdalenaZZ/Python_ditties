#!/usr/bin/env python
from __future__ import division
import pysam
import sys
import os.path
#import os.system
import argparse
import subprocess

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
parser.add_argument('-f', '--filter', default=None, dest='fi', action='store', required=False, help="Variant filter to apply")

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
myvcf.header.formats.add('FQ', ".", "Float", "Frequency of alternative allele")
myvcf.header.formats.add('GT', ".", "String", "Genotype")
myvcf.header.formats.add('QU', ".", "String", "Quality")
myvcf.header.formats.add('AD', "2", "String", "Allele depth")
#myvcf.header.formats.add('AD', "2" , "Integer", "Unfiltered allele depth")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)

# Add the new HP format to all samples
for r in myvcf:

    refb = r.ref
    try:
        altb = r.alts[0]
    # If the ALT allele is .
    except TypeError:
        altb = "N"
        r.alts=('N',)

    refd = 0
    altd = 0

    # If is SNP
    if 'CU' in r.samples[0].keys():

        for i in r.samples:
            #print ("SAM", i)
            if (r.ref == 'A'):
                refd = r.samples[i]['AU'][0]
            elif (r.ref == 'C'):
                refd = r.samples[i]['CU'][0]
            elif (r.ref == 'G'):
                refd = r.samples[i]['GU'][0]
            elif (r.ref == 'T'):
                refd = r.samples[i]['TU'][0]
            else:
                print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

            if (altb == 'A'):
                altd = r.samples[i]['AU'][0]
            elif (altb == 'C'):
                altd = r.samples[i]['CU'][0]
            elif (altb == 'G'):
                altd = r.samples[i]['GU'][0]
            elif (altb == 'T'):
                altd = r.samples[i]['TU'][0]
            # If ALT is unknown/N, add all tier1 counts together and subtract REF counts
            elif (altb == 'N'):
                altd = r.samples[i]['AU'][0]+ r.samples[i]['CU'][0]+ r.samples[i]['GU'][0]+ r.samples[i]['TU'][0]-refd
            else:
                print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

            freq=0
            if ( (altd + refd) < 1):
                print ("WARN: This SNP call is not supported by reads", altd, refd, freq, str(r) )
                #altd=altd+0.001
                #refd=refd+0.001
                #freq = altd / (altd + refd)
                #continue
            else:
                freq = altd / (altd + refd)

            # Insert genotype
            gtref = 0
            gtalt = 0

            if ( 0.1 < freq  ):
                gtalt=1
            gt = (gtref, gtalt)

            #print (gt, r.chrom, r.pos, r.ref, altb, freq, gtref, gtalt, str(refd), str(altd), i)

            # Very strangely, I'm not allowed to put strings of different lengths in, so I have to pad them out
            # so they are all exactly the same character length
            refc = "$$$$" + str(refd)
            refco = refc[-4:]
            altc = "$$$$" + str(altd)
            altco = altc[-4:]
            #ad = str(''.join(refco + ',' + altco))

            # Add the new values
            r.samples[i]['FQ'] = float(freq)
            r.samples[i]['QU'] = str(r.filter.items()[0][0])
            r.samples[i]['DP'] = (altd + refd)
            r.samples[i]['GT'] = gt
            r.samples[i]['AD'] = (refco,altco)

        vcf_out.write(r)

    # If is indel
    elif 'TIR' in r.samples[i].keys():
        for i in r.samples:
            freq=0

            # See support for indel
            if ((r.samples[i]['TAR'][0]+r.samples[i]['TIR'][0]) < 1):
                print ("WARN: This INDEL call is not supported by reads", str(r) )
            else:
                freq = r.samples[i]['TIR'][0] / (r.samples[i]['TAR'][0]+r.samples[i]['TIR'][0])

            # Insert genotype
            gtref = 0
            gtalt = 0

            if (0.1 < freq):
                gtalt=1
            gt = (gtref, gtalt)

            # Very strangely, I'm not allowed to put strings of different lengths in, so I have to pad them out
            # so they are all exactly the same character length
            refc = "$$$$" + str(r.samples[i]['TAR'][0])
            refco = refc[-4:]
            altc = "$$$$" + str(r.samples[i]['TIR'][0])
            altco = altc[-4:]

            r.samples[i]['FQ'] = float(freq)
            r.samples[i]['GT'] = gt
            r.samples[i]['QU'] = str(r.filter.items()[0][0])
            r.samples[i]['DP'] = (r.samples[i]['TAR'][0]+r.samples[i]['TIR'][0])
            r.samples[i]['AD'] = (refco, altco)

            vcf_out.write(r)

    else:
        print ("Warn: line format not as expected")



vcf_out.close()

#subprocess.call(["zcat",output, ])




exit(0)


