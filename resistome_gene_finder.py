#!/usr/bin/env python3

import os
import argparse
import fastaq
import iva
import sys
import pysam
import copy

parser = argparse.ArgumentParser(
    description = 'Test resistome finder script',
    usage = '%(prog)s [options] <db.fa> <reads1.fq> <reads2.fq> <outprefix>')
parser.add_argument('db', help='FASTA file of reference genes')
parser.add_argument('reads_1', help='Name of fwd reads fastq file')
parser.add_argument('reads_2', help='Name of rev reads fastq file')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()





def faidx_reference(fname):
    if not os.path.exists(fname + '.fai'):
        iva.common.syscall('samtools faidx ' + fname)


def get_seq_from_fasta(infile, name, outfile):
    faidx_reference(infile)
    iva.common.syscall(' '.join([
        'samtools faidx',
        infile,
        name,
        '>', outfile
    ]))


def sam_to_fastq(s):
    name = s.qname
    if s.is_read1:
        name += '/1'
    elif s.is_read2:
        name += '/2'
    else:
        raise Error('Read', name, 'must be first of second of pair according to flag. Cannot continue')

    seq = fastaq.sequences.Fastq(name, s.seq.decode(), s.qual.decode())
    if s.is_reverse:
        seq.revcomp()

    return seq




def get_read_clustered_by_ref_seq(bam):
    clusters = {}
    sam_reader = pysam.Samfile(bam, "rb")
    sam1 = None

    for s in sam_reader.fetch(until_eof=True):
        if sam1 is None:
            sam1 = s
            continue

        ref_seqs = set()
        if not s.is_unmapped:
            ref_seqs.add(sam_reader.getrname(s.tid))
        if not sam1.is_unmapped:
            ref_seqs.add(sam_reader.getrname(sam1.tid))

        for ref in ref_seqs:
            if ref not in clusters:
                clusters[ref] = []

            read1 = sam_to_fastq(sam1)
            read2 = sam_to_fastq(s)
            if read1.id.endswith('/2'):
                read1, read2 = read2, read1
            clusters[ref].append((read1, read2))
                   
        sam1 = None 

    return clusters


def write_cluster_fastq_and_ref(clusters, reference, outprefix):
    for ref in read_clusters:
        directory = outprefix + '.' + ref
        os.mkdir(directory)
        f1 = fastaq.utils.open_file_write(os.path.join(directory, 'reads_1.fq'))
        f2 = fastaq.utils.open_file_write(os.path.join(directory, 'reads_2.fq'))
        for t in read_clusters[ref]:
            print(t[0], file=f1)
            print(t[1], file=f2)
        fastaq.utils.close(f1)
        fastaq.utils.close(f2)
        get_seq_from_fasta(reference, ref, os.path.join(directory, 'ref.fa'))


def assemble_and_map_back(root_dir, assembler='spades'):
    reads1 = os.path.join(root_dir, 'reads_1.fq')
    reads2 = os.path.join(root_dir, 'reads_2.fq')
    assembly_dir = os.path.join(root_dir, 'Assembly')
    if assembler == 'spades':
        cmd = ' '.join([
            '~/bin/SPAdes-3.1.1-Linux/bin/spades.py',
            '-1', reads1,
            '-2', reads2,
            '-o', assembly_dir
        ])
        iva.common.syscall(cmd, verbose=True)
        assembly_contigs = os.path.join('Assembly', 'scaffolds.fasta')
    else:
        print('Unknown assembler:', assembler, '- cannot continue', file=sys.stderr)
        sys.exit(1)

    
    contigs = os.path.join(root_dir, 'contigs.fasta')
    os.symlink(assembly_contigs, contigs)

    bam_prefix = contigs + '.smalt'

    iva.mapping.map_reads(
        reads1,
        reads2,
        contigs,
        bam_prefix, 
        index_k=15,
        index_s=3,
        minid=0.75,
        sort=True
    )

    os.unlink(bam_prefix + '.unsorted.bam')


faidx_reference(options.db)

# map reads - keep in read pair order to easliy extract pairs for each gene
bam_prefix = options.outprefix + '.smalt'
#iva.mapping.map_reads(
#    options.reads_1,
#    options.reads_2,
#    options.db,
#    bam_prefix, 
#    index_k=20,
#    index_s=5,
#    minid=0.75,
#    sort=False
#)


read_clusters = get_read_clustered_by_ref_seq(bam_prefix + '.bam')
write_cluster_fastq_and_ref(read_clusters, options.db, options.outprefix)


for gene in read_clusters:
    if 'strA' in gene or 'strB' in gene:
        assemble_and_map_back(
            options.outprefix + '.' + gene,
            assembler='spades'
        )



# checking assembly:
#  - existing contigs have fragment covereage across whole length
#  - no read pairs suggesting scaffolding between contigs (if assembly in >1 contig)

