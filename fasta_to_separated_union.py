#!/usr/bin/env python3

import argparse
from fastaq import *


parser = argparse.ArgumentParser(
    description = 'Makes union flie from FASTA, with separator between sequences to stop reads mapping',
    usage = '%(prog)s [options] <fasta in> <outfiles prefix>')
parser.add_argument('infile', help='Name of input fasta file')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('--seqname', help='Name of output sequence [%(default)s]', default='union')
options = parser.parse_args()
seq_reader = sequences.file_reader(options.infile)
seq_out = sequences.Fasta('union', '')
gff_lines = []
total_bases = 0
filler_seq = 'n' * 100 + 'g' * 100 + 'n' * 100

for seq in seq_reader:
    seq_out.seq += seq.seq
    gff_lines.append('\t'.join([
        options.seqname,
        'Contig',
        seq.id,
        str(total_bases + 1),
        str(total_bases + len(seq)),
        '.',
        '.',
        '.',
        'color=3',
    ]))

    seq_out.seq += filler_seq
    total_bases += len(seq)

    gff_lines.append('\t'.join([
        options.seqname,
        'Fake',
        'Fake',
        str(total_bases + 1),
        str(total_bases + len(filler_seq)),
        '.',
        '.',
        '.',
        'color=2',
    ]))

    total_bases += len(filler_seq)

f = utils.open_file_write(options.outprefix + '.fa')
print(seq_out, file=f)
utils.close(f)

f = utils.open_file_write(options.outprefix + '.gff')
for line in gff_lines:
    print(line, file=f)
utils.close(f)

