#!/usr/bin/env python3.3

import argparse
import os
import copy

import utils
import fastn
import nucmer
import genome_intervals

def update_perfect_contigs(nucmer_hit, ref_fasta, contigs):
    id = nucmer_hit.ref_name + ":" + str(nucmer_hit.ref_start) + '-' + str(nucmer_hit.ref_end)
    contig = fastn.Fasta('x', ref_fasta[nucmer_hit.ref_start-1:nucmer_hit.ref_end])
    contigs[(nucmer_hit.ref_name, nucmer_hit.ref_start, nucmer_hit.ref_end)] = contig


parser = argparse.ArgumentParser(
    description = 'Takes contigs and a reference sequence. Makes a new fasta file of the contigs, but they are now perfect sequences by using the reference instead',
    usage = '%(prog)s [options] <contigs.fa> <reference.fa> <outprefix>')
parser.add_argument('--min_seq_length', type=int, help='Minimum length of contig to output [%(default)s]', default=200)
parser.add_argument('--nucmer_options', help='Options when running nucmer [%(default)s]', default='')
parser.add_argument('contigs_fa', help='Name of contigs fasta file', metavar='contigs.fa')
parser.add_argument('ref_fa', help='Name of reference fasta file', metavar='reference.fa')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

ref_seqs = {}
fastn.file_to_dict(options.ref_fa, ref_seqs)

mummer_dir = os.path.join(os.path.expanduser('~mh12'), 'bin', 'MUMmer3.23')
nucmer_exe = os.path.join(mummer_dir, 'nucmer')
delta_filter = os.path.join(mummer_dir, 'delta-filter')
show_coords = os.path.join(mummer_dir, 'show-coords')

nucmer_out_prefix = options.outprefix + '.nucmer'
nucmer_out_delta = nucmer_out_prefix + '.delta'
nucmer_out_filter = nucmer_out_prefix + '.delta-filter'
nucmer_out_coords = nucmer_out_filter + '.coords'

# run nucmer of contigs vs ref
utils.syscall(' '.join([nucmer_exe, options.nucmer_options, '-p', nucmer_out_prefix, options.ref_fa, options.contigs_fa]))
utils.syscall(' '.join([delta_filter, '-i 98 -l 180 -q', nucmer_out_delta, '>', nucmer_out_filter]))
utils.syscall(' '.join([show_coords, '-dTlro', nucmer_out_filter, '>', nucmer_out_coords]))

# load hits into hash. key=ref_name, value=another hash with key=qry_name, value=list of hit positions in that ref seq
nucmer_hits = {}
contigs_to_print = {}

nucmer_reader = nucmer.file_reader(nucmer_out_coords)

for hit in nucmer_reader:
    if hit.ref_name not in nucmer_hits:
        nucmer_hits[hit.ref_name] = {}

    if hit.qry_name not in nucmer_hits[hit.ref_name]:
        nucmer_hits[hit.ref_name][hit.qry_name] = []

    nucmer_hits[hit.ref_name][hit.qry_name].append(genome_intervals.Interval(min(hit.ref_start, hit.ref_end), max(hit.ref_start, hit.ref_end)))

# merge all the overalpping hits for each list of hits corresponding to one contig
for ref_name, d in nucmer_hits.items():
    for qry_name, hits in d.items():
        genome_intervals.merge_overlapping_in_list(hits)

        for hit in hits:
            if hit.end - hit.start + 1 >= options.min_seq_length:
                if ref_name not in contigs_to_print:
                    contigs_to_print[ref_name] = []

                contigs_to_print[ref_name].append(copy.copy(hit))

# remove any contigs that are completely contained in another contig
for ref, l in contigs_to_print.items():
    genome_intervals.remove_contained_in_list(l)


# print the final perfect contigs
f_out = utils.open_file_write(options.outprefix + '.fa')
counter = 1
last_id = None
for ref_name in sorted(contigs_to_print):
    counter = 1

    for interval in contigs_to_print[ref_name]:
        id = ':'.join([str(x) for x in [ref_name, counter, interval.start, interval.end]])
        print(fastn.Fasta(id, ref_seqs[ref_name][interval.start - 1: interval.end]), file=f_out)
        counter += 1

utils.close(f_out)

