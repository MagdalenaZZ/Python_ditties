#!/usr/bin/env python


import sys
import argparse

# Describe what the script does
parser = argparse.ArgumentParser(description='This script does blah', epilog= 'Happy use!')



parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max, help='sum the integers (default: find the max)')
parser.add_argument('-in', '--input', default=None, dest='fastq', action='store', help="fastq.gz files")
parser.add_argument('-o', "--output", default=None, dest='name', action='store', help="Sample name (output prefix)")
parser.add_argument('-b', "--bed", default=None, dest='bed', action='store', help="Bed file")
parser.add_argument('-r', default=None, dest='reference', action='store', help="Reference genome file")
parser.add_argument('bams', metavar='BAM-files', type=str , nargs='+', help='The BAM-files to summarise')

args  = parser.parse_args()


print(args)


