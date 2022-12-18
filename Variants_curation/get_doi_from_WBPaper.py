#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import subprocess
import requests
import csv
import re
import json
from wbpreader.get_fulltext import fulltext_wbp
from wbpreader.get_fulltext import fulltext_pmid
from wbpreader.id_conversions import pmid2doi


#print ("Finished importing subs")

'''

Script to parse the variants file and highlight sentences


'''

epi = ('\
    \n\
	Give an input file with variants IDs\n\
    Your genome.intervals_list is output from GATK pipeline, and a panel-of-normals (from PON.py) \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running GATK CNV pipeline from tumour and normal data', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inp', action='store', required=True, help="http://svm.textpresso.org/celegans/svm_results/20201206/120620_091120_structcorr")
#parser.add_argument('-t', '--tempdir', default='Tmpdir', dest='tmp', action='store', required=True, help="Directory for temporary results")
parser.add_argument('-o', '--output', default=None, dest='out', action='store', required=True, help="output file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

'''
if not os.path.exists(args.tmp):
    os.makedirs(args.tmp)
    print ("Created folder", args.tmp)

# Retrieve input file from URL
url = args.inp
r = requests.get(url, allow_redirects=True)
#fo = open("foo.txt", "wb")
open(args.out, 'wb').write(r.content)
#print (r.content, url)
#res.close()

# Check if output files exist
if not os.path.isfile(args.out)==True:
    print("Cannot find output file ",args.out)
    sys.exit(1)
'''

# Open output file
outf = open(args.out,'w')
        #print (com1, file=curlf)
        #outf.close()
#print ("Trying to open", args.out)


# Read in variants file
with open(args.inp, 'rt') as csvfile:
    reader = csv.reader(csvfile, delimiter ='\t')
    for row in reader:
        #print ("Row", row)
        wbp = row[0].split('\t')
        #print (row[0], wbp[0])
        seen=dict()

        # Is a WBPapaer
        if re.match( 'WBPaper', wbp[0]):
            print ("Match", wbp[0])
            ids =  pmid2doi(wbp[0])

            # Try to retrieve the full-text document from WB, and look for valid sentences
            print(ids)
            print (row,ids, sep='\t', file=outf)


        else:
            print ("NO ",wbp[0] , args.tmp)


outf.close()


quit()




