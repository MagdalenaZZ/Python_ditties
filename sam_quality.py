#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import pysam
import sys
import os
import os.path
import argparse
import re
import numpy as np

"""

Script for parsing and editing a BAM file

"""


epi = ('\
    \n\
	BAM file parser, allowing for custom advanced filtering of BAM files\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a BAM file and applies custom filtering', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='bam', action='store', required=True, help="BAM file")
#parser.add_argument('-f', '--filter-string', default=None, dest='fi', action='store', required=True, help="String of filters to apply")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.bam)==True:
    print("Cannot find input file ",args.bam)
    sys.exit(1)

# Read input
output=args.bam+".res.bam"
output2=args.bam+".splits"
o2 = open(output2, 'w')
print("Input: ",args.bam, "Output: " , output)


# read the input file
#mysam = pysam.VariantFile(args.bam, "r")
samfile = pysam.AlignmentFile(args.bam, "rb")



# create an object of new vcf file and open in to write data.
#bam_out = pysam.VariantFile(output, 'w', header=myvcf.header)

pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)

#acc = {}
M = np.zeros(300)
N = np.zeros(300)
D = np.zeros(300)
I = np.zeros(300)
S = np.zeros(300)

p = re.compile('^\d*M$')



for read in samfile.fetch():
    #if read.is_paired:
    pairedreads.write(read)
    # Filter all that are perfect mapping
    p = re.compile('^\d*M$')

    # If the read is unmapped
    if not read.cigarstring:
        pass

    # If the read is all match
    elif p.match(read.cigarstring):
        pass
        #print ("Is",read.cigarstring)
    # Anything else
    else:
        spl = re.split('(\d+)', read.cigarstring)
        spl=spl[1:]
        start = int(0)
        #print(spl)
        for i in range(1, len(spl), 2):
            if (re.match('^N$',spl[i])):
                pass
            else:
                start = start + int(spl[i - 1])
            o2.write('\t'.join([str(spl[i]), str(spl[i - 1]), str(start)] ) + os.linesep)
            #print(spl[i], spl[i - 1], file=o2)

'''
        start=int(0)
        
        for i in range(0, len(spl), 2):
            # Add to
            #print(spl[i], spl[i-1])

            # M
            if (re.match('^M$',spl[i])):
                M[start:int(spl[i - 1])] = M[start:int(spl[i-1])] +1
                start=start+int(spl[i-1])
                print (M)
            # N
            elif (re.match('^N$',spl[i])):
                N[start = N[start] +1
                print (N)
            # D
            elif (re.match('^D$',spl[i])):
                D[start:int(spl[i - 1])] = D[start:int(spl[i-1])] +1
                start=start+int(spl[i-1])
                print (D)
            # I
            elif (re.match('^I$',spl[i])):
                I[start:int(spl[i - 1])] = I[start:int(spl[i-1])] +1
                start=start+int(spl[i-1])
                print (I)
            print ("Not" , read.cigarstring, start, spl)
'''





pairedreads.close()
samfile.close()

o2.close()

quit()






