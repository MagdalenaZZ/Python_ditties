#!/usr/bin/env python3.3

import argparse
import fastn
import utils
import genome_intervals

parser = argparse.ArgumentParser(
    description = 'Makes a random genome with sequence lengths and names determined by an fai file. IMPORTANT: not really random, at the moment every base will be an A (or an N if --gaps_file used)',
    usage = '%(prog)s [options] <fai file> <outfile>')
parser.add_argument('--gaps_file', help='File of gaps, each line in the form: "chr start end" (tab separated)')
parser.add_argument('fai_file', help='Name of fai file')
parser.add_argument('outfile', help='Name of output fasta file')
options = parser.parse_args()

gaps = {}
if options.gaps_file:
    f = utils.open_file_read(options.gaps_file)
    for line in f:
        (id, start, end) = line.rstrip().split('\t')
        gap = genome_intervals.Interval(int(start) - 1, int(end) - 1)
        if id not in gaps:
            gaps[id] = []
        gaps[id].append(gap)
    utils.close(f)

f_in = utils.open_file_read(options.fai_file)
f_out = utils.open_file_write(options.outfile)

for line in f_in:
    a = line.rstrip().split()
    fa = fastn.Fasta(a[0], 'A' * int(a[1]))

    if fa.id in gaps:
        fa.seq = list(fa.seq)
        for gap in gaps[fa.id]:
            fa.seq[gap.start:gap.end+1] = ['N'] * len(gap)
        fa.seq = ''.join(fa.seq)

    print(fa, file=f_out)


utils.close(f_in)
utils.close(f_out)

