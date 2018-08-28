#!/usr/bin/env python3.3

import argparse
import os
import copy

import utils
import fastn
import sam
import external_progs
import nucmer



def map_and_parse_sam(ref_index, tags_fasta, tag_counts, log_fh):
    samfile = options.outprefix + '.maptags.sam'
    #utils.syscall('smalt map -d -1 -y 1 -f samsoft -o ' + samfile + ' ' + ref_smalt_index + ' ' + tags_fasta)
    utils.syscall(external_progs.bowtie2_align + ' -f -x ' + ref_index + ' -U ' + tags_fasta + ' -S ' + samfile)
    sam_reader = sam.file_reader(samfile)
    for sam_record in sam_reader:
        (contig_name, range) = sam_record.id.rsplit(':', 1)
        assert contig_name not in tag_counts
        if sam_record.is_mapped() \
            and sam_record.tags['AS'][1] == 0 \
            and ('XS' not in sam_record.tags or sam_record.tags['XS'][1] < 0):
                tag_counts[contig_name] = 1
        else:
            tag_counts[contig_name] = 2

    os.unlink(samfile)

def get_unique_tags(ref_index, tag_length, unique_tagged_seqs, untagged_seqs, unique_tags, log_fh, second_index=None):
    if len(untagged_seqs) == 0:
        return
    tags_fasta_fname = options.outprefix + '.tags.test.fa'
    seqs_sam = options.outprefix + '.seqs.bowtie2.sam'
    second_seqs_sam = options.outprefix + '.second_seqs.bowtie2.sam'
    fout_tags = utils.open_file_write(tags_fasta_fname)
    tags = {}
    tag_info = {}

    # make fasta file of tags
    for id, seq in untagged_seqs.items():
        tag = ''
        if len(seq) < tag_length:
            tag = fastn.Fasta(seq.id + ':1-' + str(len(seq)), seq.seq)
            tag_info[id] = [id, '1', str(len(seq)), tag.seq]
        else:
            left_coord = int(0.5 * len(seq) - 0.5 * tag_length)
            right_coord = left_coord + tag_length - 1
            tag = fastn.Fasta(seq.id + ':' + str(left_coord+1) + '-' + str(right_coord+1), seq[left_coord:right_coord + 1])
            tag_info[id] = [id, left_coord+1, right_coord+1, tag.seq]


        print(tag, file=fout_tags)
        tags[id] = copy.copy(seq)

    utils.close(fout_tags)


    # get the count of number of hits per tag from the results
    tag_counts = {}
    second_tag_counts = {}
    map_and_parse_sam(ref_index, tags_fasta_fname, tag_counts, log_fh)

    if second_index:
        map_and_parse_sam(second_index, tags_fasta_fname, second_tag_counts, log_fh)
        assert len(tag_counts) == len(second_tag_counts)

    # update the unique/non-unique tagged sequences
    for contig_name, hit_count in tag_counts.items():
        assert contig_name not in unique_tagged_seqs
        if second_index:
            second_hit_count = second_tag_counts[contig_name]
        else:
            second_hit_count = 1
        if hit_count == 1 == second_hit_count:
            unique_tagged_seqs[contig_name] = tags[contig_name]
            unique_tags.append(tag_info[contig_name])
            del untagged_seqs[contig_name]

    try:
        os.unlink(tags_fasta_fname)
    except:
        print('Error deleting file "' + tags_fasta_fname + '"', file=sys.stderr)
        sys.exit(1)


parser = argparse.ArgumentParser(
    description = 'Given a fasta file, finds a unique tag for each sequence. Can optionally give a second fasta file, and the script will check that the tag is unique in this second file as well.',
    usage = '%(prog)s [options] <in.fasta> <outprefix>')
parser.add_argument('--min_tag_length', type=int, help='Starting tag length [%(default)s]', default=50)
parser.add_argument('--max_tag_length', type=int, help='Max tag length [%(default)s]', default=1000)
parser.add_argument('--tag_step', type=int, help='Step size in tag length [%(default)s]', default=50)
parser.add_argument('--second_fasta', help='Name of a second fasta file to check for uniqueness of tags', default=None)
parser.add_argument('fasta_in', help='Name of input fasta file', metavar='in.fasta')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

untagged_seqs = {}
fastn.file_to_dict(options.fasta_in, untagged_seqs)

second_seqs = {}
unique_tags = []


seqs_index = options.outprefix + '.seqs.bowtie2.index'
#utils.syscall('smalt index -k 20 -s 10 ' + seqs_smalt_index + ' ' + options.fasta_in)
utils.syscall('bowtie2-build ' + options.fasta_in + ' ' + seqs_index)

if options.second_fasta:
    #second_seqs_smalt_index = options.outprefix + '.second_seqs_smalt_index'
    second_seqs_index = options.outprefix + '.second_seqs_bowtie2_index'
    #utils.syscall('smalt index -k 20 -s 10 ' + second_seqs_smalt_index + ' ' + options.second_fasta)
    utils.syscall('bowtie2-build ' + options.second_fasta + ' ' + second_seqs_index)
else:
    second_seqs_index=None


uniquely_tagged = {}
f_log = utils.open_file_write(options.outprefix + '.log')

for tag_length in range(options.min_tag_length, options.max_tag_length + 1, options.tag_step):
    get_unique_tags(seqs_index, tag_length, uniquely_tagged, untagged_seqs, unique_tags, f_log, second_seqs_index)
    print(tag_length,'unique:', len(uniquely_tagged), file=f_log)
    print(tag_length, 'nonunique', len(untagged_seqs), file=f_log)


non_unique_fa = options.outprefix + '.seqs-without-unique-tags.fa'
f = utils.open_file_write(non_unique_fa)
for id, seq in untagged_seqs.items():
    print(seq, file=f)
utils.close(f)

f = utils.open_file_write(options.outprefix + '.seqs-with-unique-tags.fa')
for id, seq in uniquely_tagged.items():
    print(seq, file=f)
utils.close(f)


for ext in ['1.bt2','2.bt2','3.bt2','4.bt2','rev.1.bt2','rev.2.bt2']:
    os.unlink(seqs_index + '.' + ext)

second_coords = {}
tag_counts = {}

if options.second_fasta:
    tags_tmp_fa = options.outprefix + '.tags.tmp.fa'
    f = utils.open_file_write(tags_tmp_fa)
    for t in unique_tags:
        #print('>' + t[0] + ':' + str(t[1]) + ':' + str(t[2]) + '\n' + t[3], file=f)
        print('>' + t[0] + '\n' + t[3], file=f)
    utils.close(f)
    samfile = options.outprefix + '.maptags.sam'
    #utils.syscall('smalt map -d -1 -y 1 -f samsoft -o ' + samfile + ' ' + second_seqs_smalt_index + ' ' + tags_tmp_fa)
    utils.syscall(external_progs.bowtie2_align + ' -f -x ' + second_seqs_index + ' -U ' + tags_tmp_fa + ' -S ' + samfile)
    os.unlink(tags_tmp_fa)
    sam_reader = sam.file_reader(samfile)
    for sam_record in sam_reader:
        if sam_record.is_mapped():
            tag_counts[sam_record.id] = tag_counts.get(sam_record.id, 0) + 1
            second_coords[sam_record.id] = [sam_record.rname, sam_record.pos+1]
        else:
            tag_counts[contig_name] = -1

    os.unlink(samfile)
    #os.unlink(second_seqs_smalt_index + '.smi')
    #os.unlink(second_seqs_smalt_index + '.sma')
    for ext in ['1.bt2','2.bt2','3.bt2','4.bt2','rev.1.bt2','rev.2.bt2']:
        os.unlink(second_seqs_index + '.' + ext)


f = utils.open_file_write(options.outprefix + '.tag.info.tsv')
for t in sorted(unique_tags):
    if options.second_fasta:
        if tag_counts[t[0]] == 1:
            t += second_coords[t[0]]
        elif tag_counts[t[0]] == -1:
            print('Warning: unmapped tag!', t[0], file=f_log)
        else:
            print('Warning: repetitive mapped tag!', t[0], file=f_log)
    print('\t'.join([str(x) for x in t]), file=f)

utils.close(f)
utils.close(f_log)
