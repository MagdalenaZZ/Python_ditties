#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

"""

Script for parsing and editing a VCF file

"""

import sys
import argparse
import os.path
import cyvcf2
from cyvcf2 import VCF
import pandas as pd


epi = ('\
    \n\
	Compares variant calls in a merged VCF-file\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script compares variant calls in a merged VCF-file', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default='921.somatic.vcf.gz', dest='inf', action='store', required=True, help="VCF file")
#parser.add_argument('-f', '--filter-string', default=None, dest='fi', action='store', required=True, help="String of filters to apply")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.inf)==True:
    print("Cannot find input file ",args.inf)
    sys.exit(1)





vcf = cyvcf2.VCF(args.inf)


# create a new vcf Writer using the input vcf as a template.
w = Writer(f, vcf)

# Create other output
output=args.inf+".stats"

df = pd.DataFrame(columns=vcf.samples) #creates a new dataframe that's empty


v=-1

for variant in VCF(args.inf): # or VCF('some.bcf')
    v=v+1
    alt = [item.encode('utf-8') for item in variant.ALT]
    #print(variant.REF, alt) # e.g. REF='A', ALT=['C', 'T']


    # Somehow assessing the number of alternative alleles
    if (len(alt)==1):
        pass
    # Multiple alternative alleles
    elif (len(alt) >1):
        #print("Long",str(variant))
        pass
    # No alernative variant - possibly a SN deletion
    elif (len(alt) < 1):
        #print("Null",str(variant))
        pass
    else:
        print(str(variant))








    # Get a numpy array of the depth per sample:
    dp = variant.format('DP')
    #print(str(variant))
    #print(dp[0], dp[1], dp[2], dp[3], dp.size)


    # parse through the format for each sample, splitting normal and tumour

    for i in range(0,dp.size):

        if i % 2 == 0:
            #print("N", dp[i])
            if (dp[i]>0):
                #s = int(dp[i][0])
                df.loc[v, vcf.samples[i]] = int(dp[i][0])
            pass  # Even
        else:
            #print("T", dp[i])
            if (dp[i]>0):
                #s= dp[i]
                df.loc[v, vcf.samples[i]] = int(dp[i][0])
            pass  # Odd



    # Split indels and SNPs

    if (':'.join(variant.FORMAT)=="DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50" ):
        print ("Is indel ",':'.join(variant.FORMAT), variant.REF, alt )
    elif ( ':'.join(variant.FORMAT)=="DP:FDP:SDP:SUBDP:AU:CU:GU:TU"):
        #print("Is SNV ", ':'.join(variant.FORMAT), variant.REF, alt)
        pass
    else:
        #print ("Is other ", variant.format('DP'), variant.REF, alt)
        pass
    if variant.is_indel:
        print("Indel is", ':'.join(variant.FORMAT), variant.REF, alt)



# Print output dataframe
df.to_csv(output,sep='\t',mode='w')




exit(0)
