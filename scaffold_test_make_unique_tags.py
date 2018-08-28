#!/usr/bin/env python3.3

import argparse
import os
import copy
import sys

import utils
import fastn
import sam
import external_progs
import nucmer
from scaffold_test_helper import *


def map_and_parse_sam(ref_index, query_fa, tags, qry_or_ref, ops, get_unique=True):
    samfile = ops.outprefix + '.maptags.sam'
    if get_unique:
        utils.syscall(external_progs.bowtie2_align + ' -f -x ' + ref_index + ' -U ' + query_fa + ' -S ' + samfile)
    else:
        utils.syscall(external_progs.bowtie2_align + ' -a --score-min L,0,0 -f -x ' + ref_index + ' -U ' + query_fa + ' -S ' + samfile)
    sam_reader = sam.file_reader(samfile)

    for sam_record in sam_reader:
        assert sam_record.id in tags
        if sam_record.is_mapped() and sam_record.tags['AS'][1] == 0:
            if (get_unique and (('XS' not in sam_record.tags or sam_record.tags['XS'][1] < 0))) \
               or not get_unique:
                if qry_or_ref == 'qry':
                    tags[sam_record.id].qry_hits.add(Hit(sam_record.rname, sam_record.pos, sam_record.query_strand()))
                elif qry_or_ref == 'ref':
                    tags[sam_record.id].ref_hits.add(Hit(sam_record.rname, sam_record.pos, sam_record.query_strand()))
                else:
                    print('Error parsing SAM', file=sys.stderr)
                    sys.exit(1)

    os.unlink(samfile)


def make_tag_from_fastn(seq, tag_length):
    if tag_length == 0 or len(seq) <= tag_length:
        return Tag(seq.id, 0, len(seq) - 1, seq.seq)
    else:
        left_coord = int(0.5 * len(seq) - 0.5 * tag_length)
        right_coord = left_coord + tag_length - 1
        return Tag(seq.id, left_coord, right_coord, seq[left_coord:right_coord + 1])

# returns hash: original seq id -> Tag for that sequence
def make_tags_file(filename, seqs_for_tagging, tag_length):
    tags = {}
    f = utils.open_file_write(filename)

    for id, seq in seqs_for_tagging.items():
        tag = make_tag_from_fastn(seq, tag_length)
        try:
            print(str(tag.to_fasta()), file=f)
        except:
            print('Error! str(tag.to_fasta()), tag=', tag)
            sys.exit(1)
        tags[seq.id] = copy.copy(tag)

    utils.close(f)
    return tags

parser = argparse.ArgumentParser(
    description = 'Makes tag sequences that are unique where possible to the contigs and the reference genome',
    usage = '%(prog)s [options] <in.fasta> <reference.fasta> <outprefix>')
parser.add_argument('--tag_lengths', help='comma-separated list of tag lengths to try [%(default)s]', default='50,100,200,400,600,1000,2000,5000')
parser.add_argument('--use_non_unique', action='store_true', help='If tag is not found for a sequence, use its whole length as a tag, even if the tag is not unique')
parser.add_argument('fasta_in', help='Name of input fasta file', metavar='in.fasta')
parser.add_argument('reference_fasta', help='Name of reference fasta file', metavar='in.fasta')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

ref_seqs = {}
fastn.file_to_dict(options.reference_fasta, ref_seqs)

seqs_to_be_tagged = {}
tags = {}
fastn.file_to_dict(options.fasta_in, seqs_to_be_tagged)
tags_fasta = options.outprefix + '.tmp.tags.fa'

external_progs.index_with_bowtie2(options.reference_fasta)
external_progs.index_with_bowtie2(options.fasta_in)

for tag_length in [int(x) for x in options.tag_lengths.split(',')]:
    print('---------------', tag_length, '------------')
    test_tags = make_tags_file(tags_fasta, seqs_to_be_tagged, tag_length)
    map_and_parse_sam(options.reference_fasta, tags_fasta, test_tags, 'ref', options)
    map_and_parse_sam(options.fasta_in, tags_fasta, test_tags, 'qry', options)

    for id, tag in test_tags.items():
        if tag.is_unique():
            tags[id] = copy.copy(tag)
            del seqs_to_be_tagged[id]

if options.use_non_unique and len(seqs_to_be_tagged):
    test_tags = make_tags_file(tags_fasta, seqs_to_be_tagged, 0)
    print(len(test_tags), 'non-unique')
    map_and_parse_sam(options.reference_fasta, tags_fasta, test_tags, 'ref', options, get_unique=False)
    map_and_parse_sam(options.fasta_in, tags_fasta, test_tags, 'qry', options, get_unique=False)

    for id, tag in test_tags.items():
        tags[id] = copy.copy(tag)


os.unlink(tags_fasta)


fastn.file_to_dict(options.fasta_in, seqs_to_be_tagged)

f_tags_tsv = utils.open_file_write(options.outprefix + '.tags.tsv')
f_tags_fa = utils.open_file_write(options.outprefix + '.uniquely-tagged.tags.fa')
f_unique = utils.open_file_write(options.outprefix + '.uniquely-tagged.contigs.fa')
f_bad = utils.open_file_write(options.outprefix + '.non-unique.fa')
for id, tag in sorted(tags.items()):
    if len(tag.qry_hits) == 1:
        print(tag, file=f_tags_tsv)
        print(tag.to_fasta(), file=f_tags_fa)
        print(seqs_to_be_tagged[id], file=f_unique)
    else:
        print(seqs_to_be_tagged[id], file=f_bad)


utils.close(f_tags_tsv)
utils.close(f_tags_fa)
utils.close(f_unique)
utils.close(f_bad)

