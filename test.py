#!/usr/bin/env python

"""

Nice descriptive header

"""




import sys
import argparse

# Describe what the script does
parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('integers', metavar='N', type=int, nargs='+',
                   help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')

# Save the input integer in args
args = parser.parse_args()

# Import some modules
import numpy as np
import pandas as pd
import sys

tbx = sys.argv[1]

ts = "out.txt"
#fh = ts.nopen(1)

fields = ["AC_Adj", "AN_Adj", "AC_AFR", "AN_AFR", "AC_AMR", "AN_AMR", "AN_EAS",
        "AC_EAS", "AN_FIN", "AN_NFE", "AN_OTH", "AN_SAS", "AC_FIN", "AC_NFE",
        "AC_OTH", "AC_SAS"]

for toks in fields:
	print(toks)






