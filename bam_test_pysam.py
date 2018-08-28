#!/usr/bin/env python3

import argparse
import pysam

parser = argparse.ArgumentParser(
    description = 'Counts records in a bam file',
    usage = '%(prog)s [options] <in.bam> <outfile>')
parser.add_argument('bam_in', help='Name of input bam file')
options = parser.parse_args()

samfile = pysam.Samfile(options.bam_in, "rb")
count = 0

for x in samfile.fetch(until_eof=True):
    a = x.qname
    count += 1

print(count)
