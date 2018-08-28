#!/usr/bin/env python3.3

import argparse
import fastn
import sam
import os
import sys
import utils
import external_progs

parser = argparse.ArgumentParser(
    description = 'Makes a single tabix indexed artemis plot file of proper read pair coverage from a list of BAM files',
    usage = '%(prog)s [options] <reference.fasta.fai> <outfile.gz> <in.1.bam> [in.2.bam [in.3.bam ...] ]')
parser.add_argument('fai_in', help='fai file of reference fasta', metavar='reference.fasta.fai')
parser.add_argument('outfile', help='Name of output file. End it with .gz since it will be bgzipped and has to have this for Artemis to work', metavar='outfile.bgz')
parser.add_argument('bam_files', nargs='+', help='List of BAM files', metavar='in.1.bam ...')
options = parser.parse_args()

ref_lengths = []

f = utils.open_file_read(options.fai_in)
for line in f:
    a = line.rstrip().split()
    ref_lengths.append((a[0], int(a[1])))

utils.close(f)


if not options.outfile.endswith('.gz'):
    options.outfile += '.gz'
    print('Appending .gz to output file. Now called', options.outfile, file=sys.stderr)
tmp_out = options.outfile + '.tmp.bgz'
fout = utils.open_file_write(tmp_out)
#print('#test', file=fout)


def get_cov_for_one_seq(seq_name, seq_length, bam_files, f_out):
    previous_pos = 0
    pos = 0

    try:
        fin = os.popen(external_progs.samtools + ' mpileup -r ' + seq_name + ' ' + ' '.join(bam_files))
    except IOError:
        raise Error("Error opening for reading BAM file '" + bam_file + "'")

    for line in fin:
        a = line.rstrip().split()
        pos = int(a[1])

        while previous_pos < pos - 1:
            print(a[0], previous_pos, '\t'.join(['0'] * len(bam_files)), sep='\t', file=f_out)
            previous_pos += 1

        print(a[0], pos, '\t'.join([a[i] for i in range(3,len(a),3)]), sep='\t', file=f_out)
        previous_pos = pos

    fin.close()

    pos += 1

    while pos <= seq_length:
        print(seq_name, pos, '\t'.join(['0'] * len(bam_files)), sep='\t', file=f_out)
        pos += 1

for seq, length in ref_lengths:
    get_cov_for_one_seq(seq, length, options.bam_files, fout)


utils.close(fout)
os.rename(tmp_out, options.outfile)
utils.syscall('tabix -b 2 -e 2 ' + options.outfile)
