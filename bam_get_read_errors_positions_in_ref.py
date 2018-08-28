#!/usr/bin/env python3.3

import argparse
import fastn
import sam
import utils

parser = argparse.ArgumentParser(
    description = 'Report positions in the reference where any read had an error (i.e. difference between read and reference)',
    usage = '%(prog)s [options] <in.bam> <reference.fasta> <outfile>')
parser.add_argument('bam_in', help='Name of input bam file')
parser.add_argument('fasta_in', help='Name of reference fasta file')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

sam_reader = sam.file_reader(options.bam_in)
errors = {}
ref_seqs = {}
fastn.file_to_dict(options.fasta_in, ref_seqs)

for sam_record in sam_reader:
    if sam_record.is_mapped():
        new_errors = sam_record.get_differences_from_ref(ref_seqs[sam_record.rname])
        if sam_record.rname not in errors:
            errors[sam_record.rname] = {}

        for e in new_errors:
            errors[sam_record.rname][e] = errors[sam_record.rname].get(e, 0) + 1

fout = utils.open_file_write(options.outfile)

for id in errors:
    for err in sorted(errors[id]):
        print(id, err[0]+1, '\t'.join([str(x) for x in err[1:]]), errors[id][err], sep='\t', file=fout)
utils.close(fout)
