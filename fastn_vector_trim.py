#!/usr/bin/env python3.3

import argparse
import os
import fastn
import utils
import sam
import genome_intervals
import external_progs

parser = argparse.ArgumentParser(
    description = 'Given a fasta/q file of reads, and a second fasta of vector sequences, trims the vectors off the reads. Made specifically for assembled capillary read pairs - uses BWA for mapping. Untested on short reads or unassembled read pairs',
    usage = '%(prog)s [options] <reads fasta/q> <vectors fasta> <outprefix>')
parser.add_argument('--join_distance', type=int, help='Join hits at most this many bases apart [%(default)s]', metavar='INT', default=100)
parser.add_argument('reads_in', help='Name of input fasta/q file of reads', metavar='reads fasta/q')
parser.add_argument('vectors_in', help='Name of input fasta file of vectors', metavar='vectors fasta')
parser.add_argument('outprefix', help='Prefix of names of ouput files')
options = parser.parse_args()

bwa_index = options.outprefix + '.bwa_index'
bwa_sam = options.outprefix + '.map_reads.sam'
utils.syscall(' '.join([external_progs.bwa, 'index -p', bwa_index, options.vectors_in]))
utils.syscall(' '.join([external_progs.bwa, 'bwasw -f', bwa_sam, bwa_index, options.reads_in]))

read_hit_coords = {} # id -> [(start, end), (start, end), ...]

sam_reader = sam.file_reader(bwa_sam)

for sam_record in sam_reader:
    if not sam_record.is_mapped():
        continue

    if not sam_record.is_forward_strand():
        sam_record.cigar.reverse()

    hit_start = 1
    hit_end = len(sam_record.seq)

    if sam_record.cigar.operations[0].operator == 'S':
        hit_start = sam_record.cigar.operations[0].number

    if sam_record.cigar.operations[-1].operator == 'S':
        hit_end = len(sam_record.seq) - sam_record.cigar.operations[-1].number

    if sam_record.id not in read_hit_coords:
        read_hit_coords[sam_record.id] = []

    read_hit_coords[sam_record.id].append(genome_intervals.Interval(hit_start - 1, hit_end - 1))

external_progs.bwa_index_clean(bwa_index)
os.unlink(bwa_sam)


seq_reader = fastn.file_reader(options.reads_in)
f_fa = utils.open_file_write(options.outprefix + '.fq')
f_log = utils.open_file_write(options.outprefix + '.log')

for seq in seq_reader:
    if seq.id not in read_hit_coords:
        print(seq, file=f_fa)
        print(seq.id, 'no hit', sep='\t', file=f_log)
    else:
        hits = read_hit_coords[seq.id]
        genome_intervals.merge_overlapping_in_list(hits)
        i = 0

        while i < len(hits) - 1:
            if hits[i+1].start - hits[i].end <= options.join_distance:
                hits[i] = hits[i].union_fill_gap(hits[i+1])
                hits.pop(i+1)
            else:
                i += 1

        possible_coords = []
        for i in range(len(hits)-1):
            possible_coords.append(genome_intervals.Interval(hits[i].end + 1, hits[i+1].start - 1))

        if len(possible_coords) == 1:
            seq.trim(possible_coords[0].start, len(seq) - possible_coords[0].end)
            print(seq, file=f_fa)
            print(seq.id, possible_coords[0].start + 1, possible_coords[0].end + 1, sep='\t', file=f_log)
        else:
            print(seq.id, 'multiple hits:', hits, sep='\t', file=f_log)


utils.close(f_fa)
utils.close(f_log)

