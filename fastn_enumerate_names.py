#!/usr/bin/env python3.3

import argparse
import fastn
import utils
import sys

parser = argparse.ArgumentParser(
    description = 'Renames sequences in a file by using the same prefix then adding 1,2,3... etc',
    usage = '%(prog)s [options] <prefix> <fasta/q in> <fasta/q out>')
parser.add_argument('--start_index', type=int, help='Starting number [%(default)s]', default=1)
parser.add_argument('--rename_file', help='If used, will write a file of old name to new name')
parser.add_argument('--keep_suffix', action='store_true', help='Use this to keep a /1 or /2 suffix at the end of each name')
parser.add_argument('prefix', help='Prefix to use for all sequence names')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
fout_seqs = utils.open_file_write(options.outfile)
counter = options.start_index

sequence_suffixes = ['/1', '/2']

if options.rename_file:
    fout_rename = utils.open_file_write(options.rename_file)
    print('#old\tnew', file=fout_rename)

for seq in seq_reader:
    old_id = seq.id
    seq.id = options.prefix + str(counter)

    for suff in sequence_suffixes:
        if old_id.endswith(suff):
            seq.id += suff
            break

    if options.rename_file:
        print(old_id, seq.id, sep='\t', file=fout_rename)

    print(seq, file=fout_seqs)
    counter += 1

utils.close(fout_seqs)

if options.rename_file:
    utils.close(fout_rename)

