#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
#import subprocess
#import requests
import csv
import re
#import json
#from wbpreader.readerdef import fulltext_wbp
#from wbpreader.readerdef import fulltext_pmid
#print ("Finished importing subs")

'''

Script to parse the variants file and highlight sentences


'''

epi = ('\
    \n\
	Give an input file the data saved from a VFP email\n\
     \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes a tab-delimited file from VFP data', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

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
    '''
'''

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
        print ("Row", row)
        
        if len(row)>0:
            str1 = " "
            str1 = ' '.join([str(elem) for elem in row])
            lis1=str1.split('?')
            #wbp = row[0]
            #wb= wbp[0]
            #w= wb[0]
            print ('00', lis1[0], lis1[1],lis1[2], lis1[3],lis1[4], sep='\t' )
        else:
            print ("Cannot split", len(row))
        # Is a WBPapaer
        if re.match( 'WBPaper', '1'):
           #print ("Match",wbp[0] )
           pass
        else:
            pass
            #print ("NO " )


outf.close()


quit()



