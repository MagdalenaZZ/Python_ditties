#!/usr/bin/env python3

import argparse
import pysam
import sys
from fastaq import *


parser = argparse.ArgumentParser(
    description = 'Makes read pairs for scaffolding from a BAM file that is sorted by read name',
    usage = '%(prog)s <in.bam> <reference.fa> <outprefix>')
parser.add_argument('bam_in', help='Input BAM file', metavar='in.bam')
parser.add_argument('ref_fa', help='Reference fasta file. Assumes reference.fa.fai exists', metavar='reference.fa')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('--pairs_per_read', type=int, help='Number of read pairs to make per pacbio linking read [%(default)s]', default=5)
parser.add_argument('--read_length', type=int, help='Length of reads [%(default)s]', default=250)
parser.add_argument('--insert', type=int, help='Insert size of reads [%(default)s]', default = 3000)
options = parser.parse_args()



class Hit:
    def __init__(self, sam, sam_reader, ref_lengths, ref_leeway=10, min_overlap=200, read_leeway=20):
        if sam.is_unmapped:
            self.is_unmapped = True
            return

        self.sam = sam
        self.is_unmapped = False
        self.qname = sam.qname
        self.qstart = sam.qstart
        self.qend = sam.qend
        self.q_length = sam.rlen
        self.rname = sam_reader.getrname(sam.tid)
        self.rstart = sam.pos
        self.rend = sam.aend - 1
        self.is_fwd = not sam.is_reverse
        self.r_length = ref_lengths[self.rname]
        self.set_overhangs_vars(ref_leeway=ref_leeway, min_overlap=min_overlap, read_leeway=read_leeway)


    def set_overhangs_vars(self, ref_leeway=10, min_overlap=200, read_leeway=20):
        self.overhangs_left_of_ref = self.rstart < ref_leeway and self.qstart >= min_overlap and (self.q_length - self.qend) < read_leeway
        self.overhangs_right_of_ref = (self.r_length - self.rend < ref_leeway) and (self.q_length - self.qend > min_overlap) and self.qstart < read_leeway
        self.read_start_overhangs = (self.is_fwd and self.overhangs_left_of_ref) or ((not self.is_fwd) and self.overhangs_right_of_ref)
        self.read_end_overhangs = (self.is_fwd and self.overhangs_right_of_ref) or ((not self.is_fwd) and self.overhangs_left_of_ref)
        self.overhangs = (self.overhangs_left_of_ref or self.overhangs_right_of_ref) and (self.read_start_overhangs or self.read_end_overhangs)


    def __str__(self):
        return '\t'.join([str(x) for x in [
            self.qname,
            self.qstart,
            self.qend,
            self.q_length,
            self.rname,
            self.rstart,
            self.rend,
            self.r_length,
            self.is_fwd,
            self.overhangs_left_of_ref,
            self.overhangs_right_of_ref
        ]])

    def reverse(self):
        self.qstart = self.q_length - self.qstart
        self.qend = self.q_length - self.qstart
        self.rstart = self.r_length - self.start
        self.rend = self.r_length - self.rend
        self.is_fwd = not self.is_fwd
        self.set_overhangs_vars()

    def links(self, other):
        if self.rname == other.rname \
          or self.qname != other.qname \
          or not self.overhangs \
          or not other.overhangs:
            return None

        if self.is_fwd and self.overhangs_right_of_ref:
            if (not other.is_fwd)  and other.overhangs_right_of_ref:
                return self.qname, self.q_length, self.rname, '+', self.rstart, '|', other.rname, '+', other.rstart, (self.q_length - other.qend) - self.qend
            elif other.is_fwd and other.overhangs_left_of_ref:
                return self.qname, self.q_length, self.rname, '+', self.rstart, '|', other.rname, '-', other.rend, other.qstart - self.qend

        if (not self.is_fwd) and self.overhangs_right_of_ref:
            if other.is_fwd and other.overhangs_right_of_ref:
                return self.qname, self.q_length, other.rname, '+', other.rstart, '|', self.rname, '+', self.rstart, (self.q_length - self.qend) - other.qend
            elif (not other.is_fwd) and other.overhangs_left_of_ref:
                return self.qname, self.q_length, self.rname, '+', self.rstart, '|', other.rname, '-', other.rend, (self.q_length - self.qend) - (other.q_length - other.qstart)
        if (not self.is_fwd) and self.overhangs_left_of_ref:
            if (not other.is_fwd) and other.overhangs_right_of_ref:
                return self.qname, self.q_length, self.rname, '-', self.rend, '|', other.rname, '+', other.rstart,  (other.q_length - other.qend) - (self.q_length - self.qstart)
            elif other.is_fwd and other.overhangs_left_of_ref:
                return self.qname, self.q_length, self.rname, '-', self.rend, '|', other.rname, '-', other.rend, other.qend - (self.q_length - self.qstart)

        if self.is_fwd and self.overhangs_left_of_ref:
            if other.is_fwd and other.overhangs_right_of_ref:
                return self.qname, self.q_length, self.rname, '-', self.rend, '|', other.rname, '+', other.rstart, self.qstart - other.qend
            elif (not other.is_fwd) and other.overhangs_left_of_ref:
                return self.qname, self.q_length, self.rname, '-', self.rend, '|', other.rname,'-', other.rend, self.qstart - (other.q_length - other.qstart)

        return None


def list_of_hits_to_link(l, leeway=10, min_overlap=200):
    if len(l) < 2:
        return None

    hits_at_end = []
    for hit in l:
        if hit.is_unmapped:
            continue

        if hit.sam.qlen >= min_overlap \
          and (hit.sam.rlen -  hit.sam.qlen) >= min_overlap \
          and ((hit.sam.pos < leeway) ^ (hit.r_length - hit.sam.aend < leeway)):
            hits_at_end.append(hit)

    if len(hits_at_end) == 2:
        return hits_at_end[0].links(hits_at_end[1])
    else:
        return None


ref_lengths = {}
tasks.lengths_from_fai(options.ref_fa + '.fai', ref_lengths)
links = {}

# find reads that can be used to make read pairs
current_hits = []

sam_reader = pysam.Samfile(options.bam_in, "rb")

for current_sam in sam_reader.fetch(until_eof=True):
    if len(current_hits) == 0 or current_sam.qname != current_hits[0].qname:
        if len(current_hits) > 0:
            link = list_of_hits_to_link(current_hits)

            if link is not None:
                links[link] = links.get(link, 0) + 1


            current_hits = []

    if not current_sam.is_unmapped:
        current_hits.append(Hit(current_sam, sam_reader, ref_lengths))


# load all ref seqs into memory - not efficient use of memory, so could
# probably be improved later
ref_seqs = {}
tasks.file_to_dict(options.ref_fa, ref_seqs)
insert_size = max([abs(t[9]) for t in links.keys()])
insert_size += 10 + options.pairs_per_read + 2 * options.read_length
print(insert_size, file=sys.stderr)

# make a pair of reads for each link
f1 = utils.open_file_write(options.outprefix + '.reads_1.fq')
f2 = utils.open_file_write(options.outprefix + '.reads_2.fq')
read_counter = 1

for l in links:
    refname1 = l[2]
    refname2 = l[6]
    assert refname1 in ref_seqs and refname2 in ref_seqs

    gap_length = l[9]
    outer_read_dist = int(0.5 * (insert_size - gap_length))
    assert outer_read_dist >= options.read_length + options.pairs_per_read

    orientation1 = l[3]
    orientation2 = l[7]

    def make_read(refseq, dist, orientation, length, offset, name):
        if orientation == '+':
            start = max(0, len(refseq) - dist - offset)
            end = min(start + length - 1, len(refseq) - offset - 1)
        elif orientation == '-':
            end = min(dist + offset, len(refseq) - offset -1)
            start = max(0, end - length + 1)

        return sequences.Fastq(name, refseq[start:end+1].upper(), 'I'*(end-start+1))

    reads1 = [make_read(ref_seqs[refname1], outer_read_dist, orientation1, options.read_length, i, str(read_counter + i) + '/1') for i in range(options.pairs_per_read)]
    reads2 = [make_read(ref_seqs[refname2], outer_read_dist, orientation2, options.read_length, i, str(read_counter + i) + '/2') for i in range(options.pairs_per_read)]
    read_counter += len(reads1)

    if orientation1 == '-':
        for read in reads1:
            read.revcomp()

    if orientation1 == orientation2:
        for read in reads2:
            read.revcomp()

    for i in range(len(reads1)):
        print(reads1[i], file=f1)
        print(reads2[i], file=f2)

utils.close(f1)
utils.close(f2)


