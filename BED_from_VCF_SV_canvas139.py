#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os.path
import argparse
import csv
import pysam
import re

"""

Script for parsing a small SV vcf and outputting BED

"""


epi = ('\
    \n\
	File parser, allowing counting of SV variaint from VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and converts it to a BED file', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)

# read the input file
myvcf = pysam.VariantFile(args.vcf, "r")

# create an object of new bed file and open in to write data.
output=args.vcf+".bed"
out = open(output, 'w')


for r in myvcf:

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    chr = r.chrom
    pos = r.pos
    id = str(r.id)
    varID=':'.join([id.split(":")[0],id.split(":")[1]])
    #altb = r.ref
    #altb = r.alts
    score = r.qual
    filter = r.filter
    info = r.info
    format = r.format
    samples = r.samples
    end = r.stop # r.info["END"]
    strand='.'
    svtype='NA'

    if 'SVTYPE' in r.info.keys():
        svtype = r.info.get('SVTYPE', "")

    #for key in r.info.keys():
    #    data = r.info.get(key, "")
    #    print (key,data)
    #    svtype=r.info['SVTYPE']

    # FORMAT
    #['PR', 'SR', 'RC', 'BC', 'CN', 'MCC']


    # Split out Manta calls

    if re.match(r'Manta', varID):
        pass
        #print("Manta",chr, pos, end, varID, score, strand, sep='\t')




    #print (list((r.header.filters)))
    #print(list((r.header.formats)))

    elif re.match(r'Canvas', varID):

        # Extract relevant information
        cn='NA'
        mcc='NA'
        filter='NA'
        rc='NA'

        if 'CN' in r.samples[0].keys():
            cn=r.samples[0]['CN']
        if 'MCC' in r.samples[0].keys():
            cn=r.samples[0]['MCC']
        if 'RC' in r.samples[0].keys():
            rc=r.samples[0]['RC']
        for key in r.filter.keys():
            filter=key


        print (chr,pos,end,varID,score,strand,filter,svtype,rc,cn,mcc,sep='\t', file=out)

    else:

        print("Unknown",chr, pos, end, varID, score, strand, sep='\t')



"""

# open up file
with open(args.vcf, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter ='\t')
    for row in reader:
        tsvout =""
        # pass header
        if row[0].startswith("#"):
            pass
        # CANVAS
        # only count canvas if Gain or Loss
        elif "Canvas" in row[2] and "REF" not in row[2]:
            # add to tsv file
            varID= row[2].split(":")
            end = row[7].split('END=')[-1]
            # if fist line do not add line break
            if len(tsvout) < 1:
            # sample chr start end QUAL Filter tool change
                tsvout += sample + "\t"+ row[0] + "\t"+ row[1]+ "\t"+ end +"\t"+ row[5]+"\t"+ row[6] +\
                      "\t"+ varID[0]+ "\t"+ varID[1]
            else:
                tsvout += "\n"+ sample + "\t"+ row[0] + "\t"+ row[1]+ "\t"+ end +"\t"+ row[5]+"\t"+ row[6] + \
                          "\t" + varID[0]+ "\t"+ varID[1]
            # add to tallies
            if "PASS" in row:
                canvasT += 1
                if "LOSS" in row[2]:
                    canvasLoss += 1
                elif "GAIN" in row[2]:
                    canvasGain += 1
            print(tsvout)
        # MANTA
        # tally up Manta outputs
        elif "Manta" in row[2]:
            # find SV type
            varID = row[2].split(":")
            QUAL = "."
            vcffilter = "."
            end = "."
            # identify translocations
            if "MantaBND" in row[2]:
            # No Ends in BND cases as they are translocations - add N/A to end field for the TSV
                end = "."
                QUAL = "."
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaBND += 1

            # identify Dels - can be problematic
            elif "MantaDEL" in row[2]:
                # some Dels columns out of sync with no ref or QUAL field so need to find end value from appropated field
                if len(row) ==11:
                    end = row[7].split(';')[0]
                    QUAL = row[5]
                    vcffilter = row[6]

                elif len(row) < 11:
                    #check which column end and Filter criteria are in
                    if "END=" in row[5]:
                        end = row[5].split(';')[0]
                        QUAL = "."
                        vcffilter = row[4]
                    elif "END=" in row[6]:
                        end = row[6].split(';')[0]
                        QUAL = "."
                        vcffilter = row[4]
                    else:
                        pass
                else:
                    pass
                end = "".join(i for i in end if i.isdigit())
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaDEL += 1

            # Finds invs
            elif "MantaINV" in row[2]:
                end = row[7].split(';')[0]
                end = "".join(i for i in end if i.isdigit())
                QUAL = row[5]
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaINV += 1
            # find DUPs
            elif "MantaDUP" in row[2]:
                end = row[7].split(';')[0]
                end = "".join(i for i in end if i.isdigit())
                QUAL = row[5]
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaDUP += 1

            # write to tsv. If first record do not add the \n
            if len(tsvout) < 1:
                # sample chr start end QUAL Filter tool change
                tsvout += sample + "\t" + row[0] + "\t" + row[1] + "\t" + end + "\t" + QUAL + "\t" + vcffilter + \
                          "\tManta\t" + varID[0]
            else:
                tsvout += "\n" + sample + "\t" + row[0] + "\t" + row[1] + "\t" + end + "\t" + QUAL + "\t" + vcffilter + \
                          "\tManta\t" + varID[0]
            print(tsvout)
        else:
            pass


tally_out = open(output, 'w')
tally_out.write("Sample\tcanvas_Total\tcanvas_Gain\tcanvas_Loss\tManta_Total\tManta_BND\tManta_DEL\tManta_INV\tMantaDUP\n")
tally_out.write(str(sample) + "\t" + str(canvasT)+ "\t" + str(canvasGain) + "\t" + str(canvasLoss)+ "\t" +
                str(MantaT) + "\t" + str(MantaBND)+ "\t" + str(MantaDEL)+ "\t" +str(MantaINV) +"\t" + str(MantaDUP) + "\n")
tally_out.close()

"""


out.close()

exit(0)
