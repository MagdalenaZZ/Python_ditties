#!/usr/bin/env python3

import argparse
import sam

parser = argparse.ArgumentParser(
    description = 'Counts records in a bam file',
    usage = '%(prog)s [options] <in.bam> <outfile>')
parser.add_argument('bam_in', help='Name of input bam file')
options = parser.parse_args()

sam_reader = sam.file_reader(options.bam_in)
count = 0

for sam_record in sam_reader:
    a = sam_record.id
    count += 1

print(count)
