#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from intervals import IntInterval
import sys
import os.path
import argparse
#import csv
import pysam
import re
#from subprocess import call
import subprocess
from pprint import pprint

"""

Script for parsing a small CNV vcf and outputting stats of regions covered by BED

"""


epi = ('\
    \n\
	File parser, allowing counting of CNV variant from VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and extracts regions from a BED file', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', nargs='+',action='append', type=str, help="VCF files")
parser.add_argument('-b', '--bed', default=None, dest='bed', action='store', required=True, help="BED file")
parser.add_argument('-o', '--outprefix', default=None, dest='px', action='store', required=True, help="Out prefix")
# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.bed)==True:
    print("Cannot find input file ",args.bed)
    sys.exit(1)


#print ("BED file: ",args.bed)
#print ("VCF files: ", args.vcf)


# Create summary output

#print (args.vcf[0][0])
vcf=args.vcf[0][0]
output=''.join([args.px,'.R'])
R=open(output, 'w')


# Define the visualisation area

# Parse the input BED
bed= open(args.bed, 'r')

# Get the total length of reference
bedtot=0
b={}
bedfrags=[]
bednames=[]
bededge=[]
bedfrags.append(str(0))
bededge.append(str(0))
bednames.append('')

for entry in bed:
    entry=entry.rstrip()
    ent = entry.split("\t")
    #print (ent[2],ent[1],int(ent[2])-int(ent[1])+1,sep='\t')
    ints=(int(ent[1]), int(ent[2]))
    if ent[0] in b:
        b[ent[0]][ints] = bedtot
    else:
        b[ent[0]]={}
        b[ent[0]][ints] = bedtot
    flen=(int(ent[2]) - int(ent[1]) + 1)
    bedtot += flen
    #print(ent[0],bedtot)
    bednames.append(ent[11])
    bednames.append('')
    bedfrags.append(str(bedtot-(flen/2)))
    bedfrags.append(str(bedtot))
    bededge.append(str(bedtot))
bed.close()


bedfrags=','.join(bedfrags)
bededge=','.join(bededge)
bednames='","'.join(bednames)


print('pdf(file="%s.pdf", useDingbats=FALSE)' % args.px, file=R)
print('par(mar=c(4, 4, 4, 8.1), xpd=TRUE)', file=R)
print ('plot.new()','plot.window(xlim = c(0,%s), ylim = c(0,20))' %(bedtot), sep='\n',file=R)
print('title(main="Copy numbers across region", xlab="Genomic position", ylab="Coverage")', sep='\n',file=R)
print('axis(2, pos=0, tck=1,lty=3,at = seq(0, 30, by = 2))',file=R)
print('axis(2, pos=0, at = seq(0, 30, by = 2))',file=R)
print('axis(1, pos=0, tck=1, labels=FALSE,lty=3,at=c(%s))' %(bededge),file=R)
print('axis(1, pos=0, at=c(%s),labels=c("%s"),las=2)' % (bedfrags,bednames) , file=R)
print('cols = c( \
\'#800000\',  #1 marron \n\
\'#cd0000\',  #2 red \n\
\'#e6194b\',  #3 pink red \n\
\'#ed5e81\',  #4 dark pink \n\
\'#fabebe\',  #5 pink \n\
\'#cf5c0a\',  #6 dark orange \n\
\'#f58231\',  #7 orange \n\
\'#ffcb98\',  #8 peach \n\
\'#ffe119\',  #9 yellow \n\
\'#fffac8\',  #10 pale yellow \n\
\'#bcf60c\',  #11 lime \n\
\'#87c362\', #12 light green \n\
\'#3cb44b\',  #13 green \n\
\'#808000\',  #14 bottle green \n\
\'#aaffc3\',  #15 mint green \n\
\'#46f0f0\',  #16 turqoise \n\
\'#008080\',  #17 emerald \n\
\'#8ea1e7\',  #18 light blue \n\
\'#4363d8\',  #19 blue \n\
\'#000075\',  #20 dark blue \n\
\'#480f5a\',  #21 dark purple \n\
\'#461eb4\',  #22 purple \n\
\'#e6beff\',  #23 lilac \n\
\'#f032e6\',  #24 fuschia \n\
\'#9a6324\',  #25 brown \n\
\'#5c3b15\',  #26 dark brown \n\
\'#ffffff\',  #27 white \n\
\'#cecece\',  #28 light grey \n\
\'#a6a6a6\',  #29 grey \n\
\'#6D6D6D\',  #30 dark grey \n\
\'#0C0C0C\'   #31 blackish \n\
      ) \
', file=R )

#print(b)
#print(bedtot)

## For each VCF file, do bedtools intersect

#quit()


# Make output dict
o={}


# Get comparison values for everything

args.vcf=args.vcf[0]
colo=0
legend=[]
colour=[]

for bed in args.vcf:

    #print(bed,colo)
    colo = colo + 1
    bedf = open(bed, 'r')
    legend.append(str(bed))
    cols=''.join(['cols[',str(colo),']'])
    colour.append(cols)

    # For each line in the output, pick out those that overlap
    # Save the BED range for each exon
    for el in bedf:
        el=el.rstrip()
        ele = el.split("\t")
        inss = ([int(ele[1]), int(ele[2])])
        #print(el)

        if ele[0] in b:

            # Is there an overap with any ref
            for region, value3 in b[ele[0]].items():

                maxes=max(region[0],inss[0])+b[ele[0]][region]
                mines=min(region[1],inss[1])+b[ele[0]][region]

                #try:
                    #ovls =(region & inss)
                    #print(ele)
                    #xend=b[ele[0]][region]+ ovls.length
                depth=float(ele[9])+(float(colo)/100)
                #    sstart=b[ele[0]][region]
                print('segments( %s, %s,x1=%s, y1=%s, col=cols[%s],lend=1)' % (maxes,depth ,mines, depth ,colo),file=R)
                #print(maxes,depth,mines,depth,colo,inss,region,ovls)
                #except:
                #    print("FAIL",ele[0], maxes, mines, inss, region)
                # If start and




# Make legend

leg='","'.join(legend)
col=','.join(colour)

#print("legendleg,col)
print('legend("topright",legend = c("%s"),col=c(%s), lwd = c(1,1), cex=0.3,inset=c(-0.2,0))' % (leg,col), file=R)

print('dev.off()', file=R)



R.close()
