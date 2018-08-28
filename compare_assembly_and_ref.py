#!/usr/bin/env python3.3

import argparse
import os
import copy

import utils
import fastn
import nucmer
import genome_intervals
import external_progs

parser = argparse.ArgumentParser(
    description = 'Given an assembly and a reference, reports regions of the reference that were (un)assembled, found by running nucmer. Assumes all gaps in the reference got assembled',
    usage = '%(prog)s [options] <reference.fasta> <query.fasta> <outprefix>')
parser.add_argument('reference', help='Name of reference assembly fasta file')
parser.add_argument('assembly', help='Name of assembly fasta file')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('--nucmer_ops', help='String of nucmer options, in quotes [%(default)s]', default='--maxmatch')
parser.add_argument('--df_ops', help='String of delta-filter options, in quotes [%(default)s]', default='-i 95 -l 200 -m')
options = parser.parse_args()


# returns hashes of gap positions and sequence lengths from a fasta file
def get_gaps_and_lengths(infile):
    seq_reader = fastn.file_reader(infile)
    lengths = {}
    gaps = {}

    for seq in seq_reader:
        assert seq.id not in lengths
        lengths[seq.id] = len(seq)
        gaps[seq.id] = seq.gaps()

    return lengths, gaps

# returns union of nucmer hits in ref and query
def get_nucmer_hits(coords_file):
    qry_hits = {}
    ref_hits = {}

    nucmer_reader = nucmer.file_reader(coords_file)
    for hit in nucmer_reader:
        # nucmer hits are 1-based. INside the script, use 0-based.
        start, end = sorted([hit.ref_start - 1, hit.ref_end - 1])
        if hit.ref_name not in ref_hits:
            ref_hits[hit.ref_name] = []
        ref_hits[hit.ref_name].append(genome_intervals.Interval(start, end))

        start, end = sorted([hit.qry_start - 1, hit.qry_end - 1])
        if hit.qry_name not in qry_hits:
            qry_hits[hit.qry_name] = []
        qry_hits[hit.qry_name].append(genome_intervals.Interval(start, end))

    for l in ref_hits.values():
        genome_intervals.merge_overlapping_in_list(l)
    for l in qry_hits.values():
        genome_intervals.merge_overlapping_in_list(l)


    return ref_hits, qry_hits


def make_hits_union(ids, d1, d2):
    d = {}

    for id in ids:
        d[id] = []
        if id in d1:
            d[id].extend(d1[id])
        if id in d2:
            d[id].extend(d2[id])

        genome_intervals.merge_overlapping_in_list(d[id])

    return d

# returns the total length of intervals from dictionary
# where values are lists of intervals
def total_length_from_dict(d):
    return sum([genome_intervals.length_sum_from_list(x) for x in d.values()])


def print_dict_as_tsv(d, filename):
    f = utils.open_file_write(filename)

    for id in d:
        for interval in d[id]:
            print(id, interval.start+1, interval.end+1, sep='\t', file=f)

    utils.close(f)

# run nucmer
nucmer_outprefix = options.outprefix + '.nucmer'
nucmer_script = nucmer_outprefix + '.sh'
nucmer_coords = nucmer_outprefix + '.coords'
f = utils.open_file_write(nucmer_script)
print(external_progs.nucmer, options.nucmer_ops, '-p', nucmer_outprefix, options.reference, options.assembly, file=f)
print(external_progs.delta_filter, options.df_ops, nucmer_outprefix + '.delta >', nucmer_outprefix + '.filter', file=f)
print(external_progs.show_coords, '-dTlro', nucmer_outprefix + '.filter >', nucmer_coords, file=f)
utils.close(f)
utils.syscall('bash ' + nucmer_script)

# gather the results
ref_lengths, ref_gaps = get_gaps_and_lengths(options.reference)
assembly_lengths, assembly_gaps = get_gaps_and_lengths(options.assembly)
ref_hits, assembly_hits = get_nucmer_hits(nucmer_coords)

ref_hits_and_gaps = make_hits_union(ref_lengths.keys(), ref_gaps, ref_hits)
assembly_hits_and_gaps = make_hits_union(assembly_lengths.keys(), assembly_gaps, assembly_hits)

ref_bases = sum(ref_lengths.values())
assembly_bases = sum(assembly_lengths.values())
ref_bases_assembled = total_length_from_dict(ref_hits_and_gaps)
assembly_bases_in_ref = total_length_from_dict(assembly_hits_and_gaps)

ref_gaps_sum = total_length_from_dict(ref_gaps)
assembly_gaps_sum = total_length_from_dict(assembly_gaps)

# Print summary numbers
f = utils.open_file_write(options.outprefix + '.stats.tsv')
print('#ref_bases',
    'ref_gap_bases',
    'assembly_bases',
    'assembly_gap_bases',
    'ref_bases_assembled',
    'ref_bases_assembled_pc',
    'assembly_bases_match_ref',
    'assembly_bases_match_ref_pc',
    'assembly_bases_not_gap_no_hit',
    'assembly_bases_not_gap_no_hit_pc',
    sep='\t', file=f)

print(ref_bases,
    ref_gaps_sum,
    assembly_bases,
    assembly_gaps_sum,
    ref_bases_assembled,
    round(100 * ref_bases_assembled / ref_bases,2),
    assembly_bases_in_ref,
    round(100 * assembly_bases_in_ref / assembly_bases,2),
    assembly_bases - assembly_bases_in_ref,
    round(100 * (assembly_bases - assembly_bases_in_ref) / assembly_bases,2),
    sep='\t', file=f)
utils.close(f)

# print regions that were covered
print_dict_as_tsv(assembly_hits_and_gaps, options.outprefix + '.assembly_regions_hit_ref_or_gaps.tsv.gz')
print_dict_as_tsv(ref_hits_and_gaps, options.outprefix + '.ref_regions_covered_by_assembly_or_gaps.tsv.gz')

