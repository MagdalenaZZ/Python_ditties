#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Makes a qual file from a fasta file, all quality scores the same',
    usage = '%(prog)s [options] <fasta_in> <qual_out>')
parser.add_argument('--qual_score', type=int, help='Score to use for all quality scores [%(default)s]', default=40)
parser.add_argument('fasta_in', help='Name of input fasta file')
parser.add_argument('qual_out', help='Name of output qual file')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.fasta_in)
qual_out = utils.open_file_write(options.qual_out)

for seq in seq_reader:
    print ('>' + seq.id, ' '.join([str(options.qual_score)] * len(seq)), sep='\n', file=qual_out)

utils.close(qual_out)

