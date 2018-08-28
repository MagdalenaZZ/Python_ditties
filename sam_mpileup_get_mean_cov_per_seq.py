#!/usr/bin/env python3.3

import argparse
import mpileup
import utils

parser = argparse.ArgumentParser(
    description = 'Gets the mean coverage per sequence from a samtools mpileup file',
    usage = '%(prog)s <mpileup file> <reference.fa.fai> <outfile>')
parser.add_argument('mpileup', help='Name of mpileup file', metavar='mpileup file')
parser.add_argument('fai', help='Name of fai file of reference fasta', metavar='reference.fa.fai')
parser.add_argument('outfile', help='Name of output file', metavar='outfile')
options = parser.parse_args()
means = mpileup.get_mean_coverage_per_seq(options.mpileup, options.fai)
f = utils.open_file_write(options.outfile)
for id, cov in sorted(means.items()):
    print(id, cov, sep='\t', file=f)
utils.close(f)
