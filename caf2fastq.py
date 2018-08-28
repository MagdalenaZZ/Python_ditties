#!/usr/bin/env python3.3

import sys
import fastn
import utils
import caf
import argparse

parser = argparse.ArgumentParser(
    description = 'Gets all sequences as fastq from a caf file. Writes one fastq file per ligation',
    usage = '%(prog)s <infile> <outfiles prefix>')
parser.add_argument('--min_length', type=int, help='Minimum length of sequence to output [%(default)s]', default = 50)
parser.add_argument('--trim', action='store_true', help='If trimming info is in file, then trim the sequences. Default is to not trim')
parser.add_argument('infile', help='Name of caf file to be read')
parser.add_argument('outprefix', help='Name of prefix of output files', metavar='outfiles prefix')
options = parser.parse_args()

caf_reader = caf.file_reader(options.infile)
outfiles = {} # name -> filehandle


for c in caf_reader:
    if options.trim:
        if c.clip_start is not None and c.clip_end is not None:
            c.seq.seq = c.seq.seq[c.clip_start - 1: c.clip_end]
            c.seq.qual = c.seq.qual[c.clip_start - 1: c.clip_end]
        else:
            print('Warning: no clipping info for sequence', c.id, file=sys.stderr)

    filename = '.'.join([str(x) for x in [options.outprefix, c.ligation, c.insert_min, c.insert_max, 'fq.gz']])

    if filename not in outfiles:
        outfiles[filename] = utils.open_file_write(filename)


    if len(c.seq) >= options.min_length:
        print(c.seq, file=outfiles[filename])

for fh in outfiles.values():
    utils.close(fh)

