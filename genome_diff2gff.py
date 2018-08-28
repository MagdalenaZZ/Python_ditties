#!/usr/bin/env python3.3

import argparse
import utils
import genome_diff

parser = argparse.ArgumentParser(
    description = 'Converts a genome diff file to gff file',
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('infile', help='Name of input genome diff file')
parser.add_argument('outfile', help='Name of output gff file')
options = parser.parse_args()

gd = genome_diff.GenomeDiff(options.infile)
gd.write_gff(options.outfile)
