#!/usr/bin/env python3.3

import argparse
import fastn
import utils
import sys
import os
import glob

parser = argparse.ArgumentParser(
    description = 'Assembles each read pair in a interleaved fastq file using cap3 (made for overlapping capillary reads). Assumes reads have suffixes /1 or /2 for fwd/rev',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Name of fastaq file')
parser.add_argument('outprefix', help='Prefix of output files. Writes a gzipped interleaved fastq file and a log file of read pairs that didn\'t produce 1 contig with cap3')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
f_log = utils.open_file_write(options.outprefix + '.log')
f_all_seqs = utils.open_file_write(options.outprefix + '.assembled_pairs.fq.gz')

for seq in seq_reader:
    id_prefix = seq.id[0:-2]
    files_prefix = options.outprefix + '.' + id_prefix + '.'
    (seq1, qual1) = seq.to_Fasta_and_qual()

    # make fasat and fasta.qual file for the next read pair
    try:
        next(seq_reader)
    except:
        print('Error getting mate for', seq1.id, file=sys.stderr)
        sys.exit(1)

    (seq2, qual2) = seq.to_Fasta_and_qual()
    assert seq1.id[0:-2] == seq2.id[0:-2]


    seq_file = files_prefix + 'tmp.read_pairs.fa'
    qual_file = files_prefix + 'tmp.read_pairs.fa.qual'
    f_seq = utils.open_file_write(seq_file)
    f_qual = utils.open_file_write(qual_file)

    print(seq1, file=f_seq)
    print(seq2, file=f_seq)

    print('>' + seq1.id, file=f_qual)
    print(' '.join([str(x) for x in qual1]), file=f_qual)
    print('>' + seq2.id, file=f_qual)
    print(' '.join([str(x) for x in qual2]), file=f_qual)

    utils.close(f_seq)
    utils.close(f_qual)
    utils.syscall_get_stdout('cap3 ' + seq_file)
    cap3_contigs = seq_file + '.cap.contigs'
    cap3_qual = seq_file + '.cap.contigs.qual'

    # check that cap3 output exactly 1 contig
    try:
        seq_count = fastn.count_sequences(cap3_contigs)
    except fastn.Error:
        pass

    if seq_count != 1:
        print(id_prefix, 'got', seq_count, 'contigs', file=f_log)
    else:
        # make a fastq file with the contig renamed to be the same as the reads
        cap3_fastq = files_prefix + 'tmp.cap3.fastq'
        fastn.fasta_to_fastq(cap3_contigs, cap3_qual, cap3_fastq)
        d = {}
        fastn.file_to_dict(cap3_fastq, d)
        fq = list(d.values())[0]
        fq.id = id_prefix + '.cap3.contig'
        print(fq, file=f_all_seqs)

    # tidy up files
    os.unlink(seq_file)
    os.unlink(qual_file)
    for fname in glob.glob(files_prefix + '*'):
        os.unlink(fname)

utils.close(f_log)
utils.close(f_all_seqs)

