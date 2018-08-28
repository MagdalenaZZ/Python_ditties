#!/usr/bin/env python3.3

import argparse
import os
import copy

import utils
import fastn
import nucmer
import genome_intervals


parser = argparse.ArgumentParser(
    description = 'Given a nucmer coords file, reports the regions of the query that have no nucmer hit. Doesn''t report gaps - i.e. assumes that all gaps had a hit.',
    usage = '%(prog)s [options] <nucmer.coords> <query.fasta/q> <outfile>')
parser.add_argument('nucmer_coords', help='Name of nucmer coords file', metavar='nucmer.coords')
parser.add_argument('query_file', help='Name of query fasta or fastq file', metavar='query.fasta/q')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.query_file)
seq_lengths = {}  # id -> sequence length
covered_regions = {} # id -> list of covered regions


# get query sequence lengths and gap positions - add each gap coord to the
# list of covered positions for each sequence
for seq in seq_reader:
    assert seq.id not in seq_lengths
    seq_lengths[seq.id] = len(seq)
    covered_regions[seq.id] = seq.gaps()



nucmer_reader = nucmer.file_reader(options.nucmer_coords)

for hit in nucmer_reader:
    assert hit.qry_name in seq_lengths

    # gaps are stored with coords starting from zero. Nucmer starts at 1, so need to decrement the coords
    start, end = sorted([hit.qry_start - 1, hit.qry_end - 1])
    covered_regions[hit.qry_name].append(genome_intervals.Interval(start, end))


# merge the covered regions
for l in covered_regions.values():
    genome_intervals.merge_overlapping_in_list(l)


f = utils.open_file_write(options.outfile)

# get the regions that are not covered
for id, covered in covered_regions.items():
    not_covered = []

    if len(covered) == 0:
        not_covered = [[1, seq_lengths[id]]]
    else:
        if covered[0].start != 0:
            not_covered.append([1, covered[0].start])

        for i in range(len(covered) - 1):
            not_covered.append([covered[i].end + 2, covered[i+1].start])

        if covered[-1].end + 1 != seq_lengths[id]:
            not_covered.append([covered[-1].end + 2, seq_lengths[id]])

    for t in not_covered:
        print(id, seq_lengths[id], t[0], t[1], t[1] - t[0] + 1, sep='\t', file=f)


utils.close(f)

