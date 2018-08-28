#!/usr/bin/env python3.3

import argparse
import copy
import fastn
import sam
import utils
import sys


def get_distance_to_end(sam_record, ref_lengths):
    if not sam_record.is_mapped():
        return -1
    elif sam_record.query_strand == '+':
        return ref_lengths[sam_record.rname] - (sam_record.pos + 1)
    else:
        return sam_record.ref_hit_end_position()


def get_ordered_sam_pair(in1, in2):
    if in1.is_first_of_pair():
        sam1 = in1
        sam2 = in2
    else:
        sam2 = in1
        sam1 = in2

    assert sam1.is_first_of_pair()
    assert sam2.is_second_of_pair()

    return(sam1, sam2)


parser = argparse.ArgumentParser(
    description = 'Reports the possible links between sequences in a BAM file and other handy things, by taking read pairs where each read in the pair is mapped to a different sequence. Assumes FR orientation of reads',
    usage = '%(prog)s [options] <in.bam> <outprefix>',
    epilog = r'''
Output files:
1) .links:
   showing links between contig ends. e.g. a line
   contig1  contig2  +  -  10
   means that there are 10 read pairs joining the right of contig1
   to the left of contig2 (with reads in the correct orientation)

2) .matrix.tsv:
   same info as in 1), but as a matrix (a tsv file that can be viewed in your
   faviourite spreadsheet program).

3) .orphaned_ends:
   contigs where there are reads mapped at the end with their mate unmapped
   (and the mapped read in the correct orientation). e.g.
   contig1  +  42
   means that there were 42 such reads on the right end of contig1

4) .orphaned_ends.sam:
   a file of SAM records of the reads that contributed to file 3).

5) .one_or_both_unmapped.fq:
   an interleaved fastq file of read pairs where one or both the reads were
   not mapped. Useful to do a bin assembly.
''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('bam_in', help='Name of input bam file')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('--max_insert_size', type=int, help='Maximum allowed distance between two reads in a pair. This determines whether a read pair mapped to different contigs is used or not, by looking at their separation if the two sequences were joined together [%(default)s]', default=1000, metavar='INT')
options = parser.parse_args()

seq_lengths = sam.get_sequence_lengths(options.bam_in)

sam_reader = sam.file_reader(options.bam_in)
link_hits = {}
links = {}
unpaired_ends = {}
unpaired_ends_hits = {}
f_orphan_sam = utils.open_file_write(options.outprefix + '.orphaned_ends.sam')
one_or_both_unmapped = {}
f_one_or_both_unmapped = utils.open_file_write(options.outprefix + '.one_or_both_unmapped.fq')


for sam_record in sam_reader:
    if sam_record.is_paired() and not sam_record.is_proper_pair() and sam_record.mrname != '=':
        sam_record.id = sam_record.id.replace('_left','')
        sam_record.id = sam_record.id.replace('_right','')
        if sam_record.is_mapped() and sam_record.is_mate_mapped():
            if sam_record.id not in link_hits:
                link_hits[sam_record.id] = copy.copy(sam_record)
            else:
                (sam1, sam2) = get_ordered_sam_pair(sam_record, link_hits[sam_record.id])
                assert sam1.is_first_of_pair()
                assert sam2.is_second_of_pair()

                end_distance1 = get_distance_to_end(sam1, seq_lengths)
                end_distance2 = get_distance_to_end(sam2, seq_lengths)
                pair_insert_distance = end_distance1 + end_distance2
                if pair_insert_distance <= options.max_insert_size:
                    if sam1.rname < sam2.rname:
                        key = (sam1.rname, sam2.rname, sam1.query_strand(), sam2.query_strand())
                    else:
                        key = (sam2.rname, sam1.rname, sam2.query_strand(), sam1.query_strand())

                    links[key] = links.get(key, 0) + 1

                del link_hits[sam_record.id]
        elif sam_record.is_mapped() or sam_record.is_mate_mapped():
            if sam_record.id not in unpaired_ends_hits:
                unpaired_ends_hits[sam_record.id] = copy.copy(sam_record)
            else:
                (sam1, sam2) = get_ordered_sam_pair(sam_record, unpaired_ends_hits[sam_record.id])
                assert sam1.is_first_of_pair()
                assert sam2.is_second_of_pair()

                end_distance1 = get_distance_to_end(sam1, seq_lengths)
                end_distance2 = get_distance_to_end(sam2, seq_lengths)

                if 0 <= end_distance1 <= options.max_insert_size or 0 <= end_distance2 <= options.max_insert_size:
                    print(sam1, file=f_orphan_sam)
                    print(sam2, file=f_orphan_sam)
                    key = None

                    if 0 <= end_distance1:
                        key = (sam1.rname, sam1.query_strand())
                    elif 0 <= end_distance2:
                        key = (sam2.rname, sam2.query_strand())

                    assert key is not None
                    unpaired_ends[key] = unpaired_ends.get(key, 0) + 1

                del unpaired_ends_hits[sam_record.id]

        if sam_record.is_paired() and (not sam_record.is_mapped() or not sam_record.is_mate_mapped()):
            if sam_record.id not in one_or_both_unmapped:
                one_or_both_unmapped[sam_record.id] = sam_record.to_fastn()
            else:
                if sam_record.is_first_of_pair():
                    assert one_or_both_unmapped[sam_record.id].id.endswith('2')
                    print(sam_record.to_fastn(), file=f_one_or_both_unmapped)
                    print(one_or_both_unmapped[sam_record.id], file=f_one_or_both_unmapped)
                else:
                    assert one_or_both_unmapped[sam_record.id].id.endswith('1')
                    print(one_or_both_unmapped[sam_record.id], file=f_one_or_both_unmapped)
                    print(sam_record.to_fastn(), file=f_one_or_both_unmapped)

                del one_or_both_unmapped[sam_record.id]

utils.close(f_orphan_sam)
utils.close(f_one_or_both_unmapped)

# write links to file
f = utils.open_file_write(options.outprefix + '.links')
for k in sorted(links):
    print('\t'.join([str(x) for x in k]), links[k], sep='\t', file=f)
utils.close(f)

# write a 'matrix' spreadsheet
f = utils.open_file_write(options.outprefix + '.matrix.tsv')
counts = {}
all_seqs = set()
for l in links:
    if l[0] not in counts:
        counts[l[0]] = {}

    if l[1] not in counts[l[0]]:
        counts[l[0]][l[1]] = {}

    key = l[2] + l[3]
    counts[l[0]][l[1]][key] = counts[l[0]][l[1]].get(key, 0) + 1
    all_seqs.add(l[0])
    all_seqs.add(l[1])

all_seqs = sorted(list(all_seqs))

print('\t'.join(all_seqs), file=f)
for s1 in all_seqs:
    outstring = s1
    for s2 in all_seqs:
        if s1 in counts and s2 in counts[s1]:
            outstring += '\t' + '/'.join([x + str(counts[s1][s2][x]) for x in counts[s1][s2]])
        else:
            outstring += '\t.'

    print(outstring, file=f)

utils.close(f)


# write orphaned end reads
f = utils.open_file_write(options.outprefix + '.orphaned_ends')
for k in sorted(unpaired_ends):
    print('\t'.join([str(x) for x in k]), unpaired_ends[k], sep='\t', file=f)
utils.close(f)


if len(link_hits) + len(one_or_both_unmapped):
    print('Warning! Paired reads found without their mate...', file=sys.stderr)
    for h in link_hits:
        print('Unpaired:', h, sep='\t', file=sys.stderr)
    for h in one_or_both_unmapped:
        print('Unpaired:', h, sep='\t', file=sys.stderr)
