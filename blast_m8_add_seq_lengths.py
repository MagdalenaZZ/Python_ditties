#!/usr/bin/env python3.3

import argparse
import sys
import blast

parser = argparse.ArgumentParser(
    description = 'Given a blast m8 or m9 file, makes new file with two extra columns of the lengths of the qry and ref sequence',
    usage = '%(prog)s <blast m8/9 in> <query fai file> <reference fai file> <output file>')
parser.add_argument('blast_in', help='Name of blast m8/9 file', metavar='blast m8/9 in')
parser.add_argument('qry_fai', help='Name of query fai file', metavar='query fai file')
parser.add_argument('ref_fai', help='Name of reference fai file', metavar='reference fai file')
parser.add_argument('outfile', help='Name of output file', metavar='Output file')
options = parser.parse_args()
blast.add_sequence_lengths(options.blast_in, options.ref_fai, options.qry_fai, options.outfile)

