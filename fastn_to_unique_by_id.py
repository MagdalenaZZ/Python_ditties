#!/usr/bin/env python3.3

import copy
import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Removes duplicate sequences from a fasta/q file, based on their names. If the same name is found more than once, then the longest sequence is kept',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()

# hash the sequences
seq_reader = fastn.file_reader(options.infile)
seqs = {}
ids_in_order = []

for seq in seq_reader:
    if len(seq) == 0:
       continue
    if seq.id not in seqs:
        seqs[seq.id] = copy.copy(seq)
        ids_in_order.append(seq.id)
    elif len(seqs[seq.id]) < len(seq):
        seqs[seq.id] = copy.copy(seq)

# write the output
f_out = utils.open_file_write(options.outfile)
for id in ids_in_order:
    print(seqs[id], file=f_out)
utils.close(f_out)
