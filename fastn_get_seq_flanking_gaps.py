#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Gets the sequence either side of gaps in a fasta/q file',
    usage = '%(prog)s [options] <fasta/q in> <fasta/q out>')
parser.add_argument('--left', type=int, help='Number of bases to get to left of gap [%(default)s]', default=25)
parser.add_argument('--right', type=int, help='Number of bases to get to right of gap [%(default)s]', default=25)
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.infile)
fout = utils.open_file_write(options.outfile)

print('#id', 'gap_start', 'gap_end', 'left_bases', 'right_bases', sep='\t', file=fout)

for seq in seq_reader:
    gaps = seq.gaps()

    for gap in gaps:
        left_start = max(gap.start - options.left, 0)
        right_end = min(gap.end + options.right + 1, len(seq))
        print(seq.id, gap.start + 1, gap.end + 1, seq.seq[left_start:gap.start], seq.seq[gap.end + 1:right_end], sep='\t', file=fout)


utils.close(fout)

