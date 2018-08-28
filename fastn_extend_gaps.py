#!/usr/bin/env python3.3

import argparse
import fastn
import utils
import sys

parser = argparse.ArgumentParser(
    description = 'Extends the length of all gaps (and trims the start/end of chromosomes) in a fasta/q file. Does this by replacing a set number of bases either side of each gap with Ns',
    usage = '%(prog)s [options] <fasta/q in> <fasta/q out>')
parser.add_argument('--trim_number', type=int, help='Number of bases to trim around each gap, and off ends of each sequence [%(default)s]', default=100)
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
fout = utils.open_file_write(options.outfile)

for seq in seq_reader:
    if len(seq) < 2 * options.trim_number:
        print('Warning! Ignoring sequence', seq.id, 'because it\'s too short', file=sys.stderr)
        continue

    gaps = seq.gaps()
    bases = list(seq.seq)

    # extend the length of each gap
    for gap in gaps:
        left_start = max(gap.start - options.trim_number, 0)
        right_end = min(gap.end + options.trim_number + 1, len(seq))

        for i in range(left_start, gap.start):
            bases[i] = 'N'

        for i in range(gap.end, right_end):
            bases[i] = 'N'

    seq.seq = ''.join(bases)

    # trim start/end bases and tidy up any resulting Ns at either end of the trimmed seq
    seq.trim(options.trim_number, options.trim_number)
    seq.trim_Ns()

    print(seq, file=fout)

utils.close(fout)

