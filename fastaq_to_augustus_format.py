#!/usr/bin/env python3

import argparse
import sys
from fastaq import *

parser = argparse.ArgumentParser(
    description = 'Converts fastq file to augstus fastq format',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

seq_reader = sequences.file_reader(options.infile)
f_out = utils.open_file_write(options.outfile)

for seq in seq_reader:
    if seq.id.endswith('/1'):
        seq.id = seq.id[:-2] + '-1'
    elif seq.id.endswith('/2'):
        seq.id = seq.id[:-2] + '-2'
    else:
        print("Didn't get seq id ending with /1 or /2. Cannot continue", seq.id, file=sys.stderr)

    print(seq, file=f_out)

utils.close(f_out)
