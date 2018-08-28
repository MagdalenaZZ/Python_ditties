#!/usr/bin/env python3

import argparse
import utils
import genome_intervals

parser = argparse.ArgumentParser(
    description = 'Takes a file of genome regions, merges any overlapping regions',
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('infile', help='Name of input file of regions')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()


def file2regions(fname):
    regions = {}

    f = utils.open_file_read(fname)

    for line in f:
        if line.startswith('#'):
            continue

        (chr, start, end) = line.rstrip().split()
        if chr not in regions:
            regions[chr] = []

        regions[chr].append(genome_intervals.Interval(start, end))

    utils.close(f)
    return regions


regions = file2regions(options.infile)
f = utils.open_file_write(options.outfile)


for chr, l in sorted(regions.items()):
    genome_intervals.merge_overlapping_in_list(l)
    for region in l:
        print(chr, region.start, region.end, sep='\t', file=f)

utils.close(f)

