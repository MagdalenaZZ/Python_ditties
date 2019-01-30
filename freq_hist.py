#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os.path
import argparse
#import csv
import pysam
import re
#from subprocess import call
import subprocess

import numpy as np
import seaborn as sns
"""

Script for parsing a small CNV vcf and outputting stats of regions covered by BED

"""


epi = ('\
    \n\
	File parser, allowing counting of CNV variaint from VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a tab-delim file and outputs a freq hist', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inf' ,action='store', type=str, required=True, help="Tab-delim file with column(s) with numbers")
parser.add_argument('-c', '--column', default=None, dest='col', action='store', type=str,required=True, help="Which column to draw freq-hist from")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.inf)==True:
    print("Cannot find input file ",args.inf)
    sys.exit(1)





#out=open(args.inf + '.res', 'w')

for elem in args.col.split(','):

    d = np.empty(0)
    elem=int(elem)
    #print("ELEM:",elem)
    inf= open(args.inf, 'r')
    
    for line in inf:
        #print("LINE",line)
        a=line.split('\t')
        try:
            val=float(a[elem-1])
            #if isinstance(line[elem+1], (int, long, float, complex)):
            d = np.append(d,val)
            #print("NUM",val)
        except:
            continue
            #print("NOM",a[elem-1])
    inf.close()
    sns_plot=sns.distplot(d)
    fig = sns_plot.get_figure()
    fig.savefig(args.inf+'.'+str(elem)+".png")


exit(0)
