#!/usr/bin/env python3

import argparse
import iva

parser = argparse.ArgumentParser(
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('infile', help='Name of input file')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()


a = iva.assembly.Assembly(contigs_file=options.infile)
a._merge_overlapping_contigs(list(a.contigs.keys()))
a.write_contigs_to_file(options.outfile)
