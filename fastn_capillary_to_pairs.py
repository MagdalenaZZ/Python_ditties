#!/usr/bin/env python3.3

import copy
import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Given a fasta/q file of capillary reads, makes an interleaved file of read pairs (where more than read from same ligation, takes the longest read) and a file of unpaired reads. Replaces the .p1k/.q1k part of read names to denote fwd/rev reads with /1 and /2',
    usage = '%(prog)s <infile> <outfiles prefix>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outprefix', help='Prefix of output files', metavar='outfiles prefix')
options = parser.parse_args()

# hash the sequences, only taking longest where an end has been sequenced more than once
seq_reader = fastn.file_reader(options.infile)
fwd_seqs = {}
rev_seqs = {}
unpaired_seqs = {}

for seq in seq_reader:
    id_info = seq.split_capillary_id()
    if id_info['dir'] == 'fwd':
        seq.id = id_info['prefix'] + '/1'
        h = fwd_seqs
    elif id_info['dir'] == 'rev':
        seq.id = id_info['prefix'] + '/2'
        h = rev_seqs
    else:
        seq.id = id_info['prefix']
        h = unpaired_seqs

    key = id_info['prefix']

    if key not in h or len(h[key]) < len(seq):
        h[key] = copy.copy(seq)

# write the output files
f_pe = utils.open_file_write(options.outprefix + '.paired.fq.gz')
f_up = utils.open_file_write(options.outprefix + '.unpaired.fq.gz')

for id in fwd_seqs:
    if id in rev_seqs:
        print(fwd_seqs[id], file=f_pe)
        print(rev_seqs[id], file=f_pe)
        del rev_seqs[id]
    else:
        print(fwd_seqs[id], file=f_up)

for seq in rev_seqs.values():
    print(seq, file=f_up)

for seq in unpaired_seqs.values():
    print(seq, file=f_up)

utils.close(f_pe)
utils.close(f_up)
