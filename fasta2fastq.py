#!/usr/bin/env python3.3

import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Given a fasta and qual file, makes a fastq file',
    usage = '%(prog)s <fasta in> <qual in> <fastq out>')
parser.add_argument('fasta', help='Name of fasta file to be read', metavar='fasta in')
parser.add_argument('qual', help='Name of qual file to be read', metavar='qual in')
parser.add_argument('outfile', help='Name of output fastq file', metavar='fastq out')
options = parser.parse_args()
fastn.fasta_to_fastq(options.fasta, options.qual, options.outfile)
