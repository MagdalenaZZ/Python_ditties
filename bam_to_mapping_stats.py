#!/usr/bin/env python3.3

import argparse
import fastn
import sam
import utils

parser = argparse.ArgumentParser(
    description = 'Gathers stats from a BAM file',
    usage = '%(prog)s [options] <in.bam> <outfile>')
parser.add_argument('bam_in', help='Name of input bam file')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

sam_reader = sam.file_reader(options.bam_in)

max_edit = 20
edit_distances = {x:0 for x in range(max_edit+1)}
max_soft_clip = 20
soft_clipped = {x:0 for x in range(max_soft_clip+1)}


stats = ['Reads',
         'Paired',
         'First of pair',
         'Second of pair',
         'Mapped',
         'Unmapped',
         'Duplicate',
         'Proper pair',
         'Perfectly aligned',
         'Perfectly aligned and proper pair']

counts = {x:0 for x in stats}



for sam_record in sam_reader:
    counts['Reads'] += 1

    if sam_record.is_mapped():
        counts['Mapped'] += 1
        edit_dist = 0
        soft_clip = sam_record.cigar.soft_clipped_bases()

        if 'NM' in sam_record.tags:
            edit_dist = sam_record.tags['NM'][1]
            edit_distances[min(max_edit, edit_dist)] += 1

        soft_clipped[min(max_soft_clip, soft_clip)] += 1

        if soft_clip == edit_dist == 0:
            counts['Perfectly aligned'] += 1

        if sam_record.is_proper_pair():
            counts['Proper pair'] += 1

            if soft_clip == edit_dist == 0:
                counts['Perfectly aligned and proper pair'] += 1


    else:
        counts['Unmapped'] += 1

    if sam_record.is_paired():
        counts['Paired'] += 1
        if sam_record.is_first_of_pair():
            counts['First of pair'] += 1
        else:
            counts['Second of pair'] += 1

    if sam_record.is_duplicate():
        counts['Duplicate'] += 1

f = utils.open_file_write(options.outfile)

for s in stats:
    print(counts[s], round(100.0 * counts[s] / counts['Reads'], 2), s, sep='\t', file=f)

for x in sorted(edit_distances):
    print('edit_distance', x, edit_distances[x], sep='\t', file=f)

for x in sorted(soft_clipped):
    print('soft_clipped', x, soft_clipped[x], sep='\t', file=f)

utils.close(f)

