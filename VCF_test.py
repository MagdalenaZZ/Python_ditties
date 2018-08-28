#!/usr/bin/env python:w
from __future__ import division


"""

Script for parsing and editing a VCF file

"""

import vcf
import sys
import argparse
import os.path
import vcf



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
o1=args.vcf+".pyVCF.vcf"
output=args.vcf+".filtered.vcf"
#print(args.vcf, output)

vcf_trans = vcf.Reader(filename=args.vcf)
vcf_reader = vcf.Reader(filename=args.vcf)

#vcf_writer = vcf.Writer(open('/dev/null', 'w'), vcf_reader)
vcf_wt = vcf.Writer(open(o1, 'w'), vcf_trans)
vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

# Write the file with nothing done, just the changes introduced by pyVCF
for r in vcf_trans:
    vcf_wt.write_record(r)


for r in (vcf_reader):
    try:
        #print(r.FORMAT)

        if (r.FORMAT == 'DP:FDP:SDP:SUBDP:AU:CU:GU:TU'):

            #print(r.REF,r.ALT)
            refb=r.REF
            altb=r.ALT
            refd=0
            altd=0

            if (r.REF=='A'):
                refd=r.samples[0]['AU'][0]
            elif (r.REF=='C'):
                refd = r.samples[0]['CU'][0]
            elif (r.REF=='G'):
                refd = r.samples[0]['GU'][0]
            elif (r.REF=='T'):
                refd = r.samples[0]['TU'][0]
            else:
                print("WARN: ", r.REF , " is not A,C,G or T :", r.ID)

            if (r.ALT[0]=='A'):
                altd = (r.samples[0]['AU'][0])
            elif (r.ALT[0]=='C'):
                altd = r.samples[0]['CU'][0]
            elif (r.ALT[0]=='G'):
                altd = r.samples[0]['GU'][0]
            elif (r.ALT[0]=='T'):
                altd = r.samples[0]['TU'][0]
            else:
                print("WARN: ", r.REF , " is not A,C,G or T :", r.ID)

            r.FORMAT='DP:FDP:SDP:SUBDP:AU:CU:GU:TU:GT:FQ'
            freq=altd/(altd+refd)
            #r.samples[0]['FQ'] = freq
            #r.FORMAT['FQ']=freq
            #print("ref", refd, r.REF ,"samp", r.samples[0])
            #print("alt", altd, r.ALT, "samp", r.samples[0])
            #print(freq, refd, altd)
            #print ("INFO",r.FORMAT)
            print ("BF",r.samples[0])
            r.samples[0].add_field('FOO', 'BAR')
            print ("AF",r.samples[0])
            #for s in r.samples:
            #    s.add_field('FOO', 'BAR')
            #    #s['DP']="Hi"
            #    print("New ",s)

            #print ("SNP",r.FORMAT)
        elif (r.FORMAT=='DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50'):
            pass
            #print ("INDEL", r.FORMAT)
        else:
            print ("Warn: not sure which type of variant this is ", r.FORMAT)

        vcf_writer.write_record(r)

    except:
        pass
        #print ("EXCEPT: ", r.CHROM, r.POS)

