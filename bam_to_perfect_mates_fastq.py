#!/usr/bin/env python3.3

import argparse
import fastn
import sam
import utils
import sys


parser = argparse.ArgumentParser(
    description = 'Filters bam for reads where both of the pair are mapped. Writes pair of fastq files, where the sequence of each read is taken from the reference, so it is perfect',
    usage = '%(prog)s [options] <reference.fasta> <in.bam> <outprefix>')
parser.add_argument('--length', type=int, help='length of output reads [%(default)s]', default=100)
parser.add_argument('ref_in', help='Name of input reference fasta file')
parser.add_argument('bam_in', help='Name of input bam file')
parser.add_argument('outprefix', help='Prefix of output files. Will write outprefix_1.fq.gz and outprefix_2.fq.gz')
options = parser.parse_args()

ref_seqs = {}
fastn.file_to_dict(options.ref_in, ref_seqs)

sam_reader = sam.file_reader(options.bam_in)
reads = {}


f1 = utils.open_file_write(options.outprefix + '_1.fq')
f2 = utils.open_file_write(options.outprefix + '_2.fq')

for s in sam_reader:
    if s.is_mapped() and s.is_mate_mapped():
        end = min(s.pos + options.length - 1, len(ref_seqs[s.rname]) - 1)
        start = max(end - options.length + 1, 0)
        read = fastn.Fastq(s.id, ref_seqs[s.rname][start:end+1] , 'I' * (end - start + 1))

        if not s.is_forward_strand():
            read.revcomp()

        if read.id in reads:
            if s.is_first_of_pair():
                read.id += '/1'
                mate = reads[s.id]
                mate.id += '/2'
                print(read, file=f1)
                print(mate, file=f2)
            else:
                read.id += '/2'
                mate = reads[s.id]
                mate.id += '/1'
                print(read, file=f2)
                print(mate, file=f1)

            del reads[s.id]
        else:
            reads[s.id] = read


print('remaining reads:', len(reads), file=sys.stderr)

utils.close(f1)
utils.close(f2)
