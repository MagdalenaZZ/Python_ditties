#!/Users/mz3/miniconda/envs/py376/bin python
from __future__ import division
from __future__ import print_function
#import sys
import os.path
import argparse
from argparse import ArgumentParser
#import re
#import gffutils
from wbpreader.translation import translate_dict
from wbpreader.translation import vars_from_GFF
#from Bio import SeqIO
#from pyfaidx import Fasta

'''
Script for generating variants Ace from mutation
'''




epi = ('\
    \n\
	Takes input variant and transcripts, and generates ace-format data\n\
    The inupt file should contain 3 columns: \n\
    1. Gene name eg WBGene00007063\n\
    2. Transcript name eg 2L52.1b\n\
    3. Variant code eg V600E K342* or W256amber\n\
    Example: \n\
    WBGene00007065	3R5.1b	    Q214L \n\
    WBGene00007064	2RSSE.1a	L10K \n\
    WBGene00007064	2RSSE.1b	A4ochre \n\
    WBGene00007064	2RSSE.1b	S7amber \n\
    \n\
    Note: The genome fasta, protein fasta and GFF needs to be from the same source, and pertain  \n\
    to the genome and annotation the variant originates from \n\
    Note: If the GFF contains variants with VEP results, a smaller ACE file will be created for  \n\
    those new ones overlapping with existing ones\n\
    Note: Wormbase GFFs have to be pre-processed using the script \n\
    \n\
    \n\
    \n\
')



# Describe what the script does
parser = argparse.ArgumentParser(description='This script generates ace files from mutations', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

#print ("HERE2")



# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inp', action='store', required=True, help="File with genes and variants eg K342* or W256amber")
parser.add_argument('-p', '--protein', default=None, dest='pro', action='store', required=True, help="Reference protein fasta file")
parser.add_argument('-gff', '--gff', default=None, dest='gff', action='store', required=True, help="GFF file with full path, containing genes and possibly existing variants")
parser.add_argument('-g', '--genome', default=None, dest='genome', action='store', required=True, help="Reference genome fasta file")



args  = parser.parse_args()

#print("ARGS",args)

def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        print("The file does not exist!", arg)
        exit(1)
    else:
        print ("File exists", arg)
        #return open(arg, 'r')  # return an open file handle

#print ("Do", args, args.inp)

is_valid_file(args, args.inp)
is_valid_file(args, args.pro)
is_valid_file(args, args.gff)
is_valid_file(args, args.genome)

# Check if input files exist
if not os.path.isfile(args.inp)==True:
    print("Cannot find input file ",args.inp)
    sys.exit(1)

# read the input file
inp = open(args.inp, 'r')

# Get a translation table
#prots=translate_dict('all')
#print (prots)

# Get all known variants from GFF
known_vars=vars_from_GFF(args.gff)



exit(0)


