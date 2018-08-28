#!/usr/bin/env python3.3

import argparse
import utils
import genome_intervals

parser = argparse.ArgumentParser(
    description = 'Takes 2 files of genome regions, returns the intersection',
    usage = '%(prog)s [options] <infile1> <infile2> <outfile>')
parser.add_argument('infile1', help='Name of first file of regions')
parser.add_argument('infile2', help='Name of second file of regions')
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


regions1 = file2regions(options.infile1)
regions2 = file2regions(options.infile2)

f = utils.open_file_write(options.outfile)

for chr in regions1:
    if chr in regions2:
        intersection = genome_intervals.intersection(regions1[chr], regions2[chr])
        for region in intersection:
            print(chr, region.start, region.end, sep='\t', file=f)

utils.close(f)

