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

Script for parsing a VCF with ALT in INFO fields and outputting hom/het/mixed and confidence

"""


epi = ('\
    \n\
	File parser,  VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and scores ALT-presence', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF.gz file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist and create an index, if the index does not exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)
if not (os.path.isfile(args.vcf+".tbi")==True or os.path.isfile(args.vcf+".csi")==True ):
    call(["bcftools","index",args.vcf])

# Merge the file with ALT variants
alts='/home/mzarowiecki/scratch/REF/allASDPs.SNV.50_10.valid.vcf.gz'

call(["bcftools", "merge","--force-samples","-O","z","-o",args.vcf+".asdp.vcf.gz",args.vcf,alts])




# read the input file
myvcf = pysam.VariantFile(args.vcf+".asdp.vcf.gz", "r")

# create an object of new bed file and open in to write data.
output=args.vcf+".asdp.res.vcf.gz"
out = open(output +'.review', 'w')
vaf = open(output +'.vaf', 'w')



myvcf.header.info.add("ALT", "1", "String", "Is variant on ALT or primary")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)



# First parse through VCF file and pick out all SNVs in ALTs.

res={}
# res[alt][hom/het/amb/none]=number


for r in myvcf:

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    chr = r.chrom
    pos = r.pos
    id = str(r.id)
    #varID=':'.join([id.split(":")[0],id.split(":")[1]])
    #altb = r.ref
    #altb = r.alts
    score = r.qual
    filter = r.filter
    info = r.info
    format = r.format
    samples = r.samples
    end = r.stop # r.info["END"]
    strand='.'
    r.info['ALT']="NAN"
    
    # Test if asdp exists
    #if 'AL' in r.info.keys():
        #print (r.info['AL'])
    # Test if we have that variant
#    if 'DP' in r.format.keys():
#        print (r.format.keys())


    # Test if both info exists
    if 'AL' in r.info.keys() and 'AD' in r.format.keys():


        # Filter complex variants
        if len( r.samples[0]['AD']) > 2:
            continue
            #print(r.samples[0]['AD'])

        # Classify variant
        fq = r.samples[0]['AD'][1]/(r.samples[0]['AD'][1]+r.samples[0]['AD'][0])

        # If it is certainly a REF
        if fq<0.1:
            r.info['ALT']="REF"
        # If it is certainly a HOM ALT
        elif fq>=0.9:
            r.info['ALT']="HOM"
        # If it is certainly a HET PA/ALT
        elif fq>0.30 and fq<0.7:
            r.info['ALT']="HET"
        else:
            r.info['ALT']="AMB"

        print(r.info['AL'], r.info['ALT'],fq,r.samples[0]['AD'][0],r.samples[0]['AD'][1], file=vaf, sep='\t')
        vcf_out.write(r)
        

    elif 'AL' in r.info.keys():
        r.info['ALT']="REF"

    else:
        continue

# Populate summary
    if r.info['AL'] in res:
            
        if r.info['ALT'] in res[r.info['AL']]:
            res[r.info['AL']][r.info['ALT']] += 1
        else:
            res[r.info['AL']][r.info['ALT']] = 1

    else:
        res[r.info['AL']]={}
        #res[r.info['AL']]['NON']=0
        res[r.info['AL']]['HOM']=0
        res[r.info['AL']]['HET']=0
        res[r.info['AL']]['AMB']=0
        res[r.info['AL']]['REF']=0
        res[r.info['AL']][r.info['ALT']] = 1



# Then summarise the score

print ("Alt","Max","AMB","HET","HOM","REF","Sum","Verdict",sep='\t',end='',file=out)

for key in sorted(res):
    #print(key,end='\t',file=out) 
    max2 = max(res[key], key=res[key].get)
    print('\n',key,'\t',max2,sep='',end='',file=out)
    sum2=0
    for key2 in sorted(res[key]):
                print('\t',res[key][key2],sep='',end='',file=out) 
                sum2+=res[key][key2]

    print('\t',sum2,sep='',end='',file=out)


    # Try to classify
    fq=res[key][max2]/sum2

    # The max value is larger than all the others together
    if fq>0.8:
        print('\t',max2,sep='',end='',file=out)
    elif fq>0.6 and max2=='HET':
        print('\t',max2,sep='',end='',file=out)
    else:
        print('\t','AMBI',sep='',end='',file=out)


    

print ('\n',sep='\t',end='',file=out)







 
out.close()
vaf.close()

vcf_out.close()

exit(0)
