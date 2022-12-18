#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
#import subprocess
#import requests
#import csv
import re
#import json
#from wbpreader.readerdef import fulltext_wbp
#from wbpreader.readerdef import fulltext_pmid
#print ("Finished importing subs")

'''

Script to parse an Alliance schema and show it as a json form schema


'''

epi = ('\
    \n\
	Give an input file Alliance json\n\
     \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script takes a json schema and turns it into a jsonform compliant schema', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inp', action='store', required=True, help="VFP email")
#parser.add_argument('-t', '--tempdir', default='Tmpdir', dest='tmp', action='store', required=True, help="Directory for temporary results")
parser.add_argument('-o', '--output', default=None, dest='out', action='store', required=True, help="output file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()
