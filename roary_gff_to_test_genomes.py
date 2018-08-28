#!/usr/bin/env python3

import copy
import random
import argparse
import pyfastaq


parser = argparse.ArgumentParser(
    description = 'Make simulated genomes for testing pan-genome pipelines',
    usage = '%(prog)s [options] <core> <genomes> <infile> <prefix of output files>')
parser.add_argument('--spacing_ns', type=int, help='Number of Ns to add between output genes [%(default)s]', default=50, metavar='INT')
parser.add_argument('--unique', type=int, help='Number of unique genes to add to each genome [%(default)s]', default=1, metavar='INT')
parser.add_argument('--step', type=int, help='Determines drop off of genes not present in all samples [%(default)s]', default=1, metavar='INT')
parser.add_argument('core', type=int, help='Number of core genes to put into each genome')
parser.add_argument('genomes', type=int, help='Number of genomes to output')
parser.add_argument('infile', help='Name of input GFF3 file')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()


def get_longest_seq(infile):
    all_seqs = {}
    pyfastaq.tasks.file_to_dict(options.infile, all_seqs)
    max_len = -1
    longest = None
    for seq in all_seqs.values():
        if max_len < len(seq):
            longest = seq.id
            max_len = len(seq)
   
    assert longest is not None
    return all_seqs[longest]


def strip_id(s):
    fields = s.rstrip().split(';')
    new_fields = [x for x in fields if not x.startswith('ID=')]
    assert len(fields) == len(new_fields) + 1
    return ';'.join(new_fields)


def get_cds_annotations(infile, seqname):
    annotations = []
    f = pyfastaq.utils.open_file_read(infile)
    for line in f:
        if line.startswith('##FASTA'):
            break

        if line.startswith('#'):
            continue

        data = line.rstrip().split('\t')
        if data[2] != 'CDS' or data[0] != seqname:
            continue

        data[-1] = 'gene.' + str(len(annotations) + 1) + ';' + strip_id(data[-1])
        annotations.append(data)

    pyfastaq.utils.close(f)
    return annotations


def cds_gene_seqs(annotations, seq):
    gene_seqs = []
    for a in annotations:
        assert seq.id == a[0]
        
        start = int(a[3]) - 1
        end = int(a[4]) - 1
        assert start < end
        assert end < len(seq)
        gene_seqs.append(pyfastaq.sequences.Fasta('gene.' + str(len(gene_seqs) + 1), seq.seq[start:end+1]))

    assert len(gene_seqs) == len(annotations)
    return gene_seqs


def make_fake_fasta(gene_seqs, indexes, name, spacing_ns):
    new_seq = []
    for i in indexes:
        new_seq.append(gene_seqs[i].seq)
        new_seq.append('N' * spacing_ns)
    new_seq = ''.join(new_seq[:-1])
    return pyfastaq.sequences.Fasta(name, new_seq)


def write_gff(gene_seqs, annotations, indexes, name, outfile, spacing_ns):
    assert max(indexes) < len(gene_seqs)
    assert len(annotations) == len(gene_seqs)
    f = pyfastaq.utils.open_file_write(outfile)
    print('##gff-version 3', file=f)
    fa = make_fake_fasta(gene_seqs, indexes, name, spacing_ns)
    print('##sequence-region', name, 1, len(fa), file=f)

    position = 1
    for i in indexes:
        a = copy.copy(annotations[i])
        a[3] = position
        a[4] = position + len(gene_seqs[i]) - 1
        a[-1] = 'ID=' + name + '.' + a[-1]
        print('\t'.join([str(x) for x in a]), file=f)
        position += len(gene_seqs[i]) + spacing_ns
        

    print('##FASTA', file=f)
    print(fa, file=f)
    pyfastaq.utils.close(f)


def make_indexes(opts):
    core_indexes = list(range(options.core))
    indexes = []
    free_index = options.core
    for i in range(opts.genomes):
        l = copy.copy(core_indexes)
        l += list(range(free_index, free_index + options.unique, 1))
        free_index += options.unique
        indexes.append(l)

    for i in range(opts.genomes - 1, 0, -opts.step):
        indexes_that_get_gene = random.sample(range(opts.genomes), i)
        for j in indexes_that_get_gene:
            indexes[j].append(free_index)
        free_index += 1

    return indexes


seq = get_longest_seq(options.infile)
annotations = get_cds_annotations(options.infile, seq.id)
gene_seqs = cds_gene_seqs(annotations, seq)
assert len(annotations) == len(gene_seqs)
indexes = make_indexes(options)
assert len(indexes) == options.genomes
f = pyfastaq.utils.open_file_write(options.outprefix + '_genome2gene')
print('genome', 'genes', sep='\t', file=f)

for i in range(options.genomes):
    genome = 'genome_' + str(i+1)
    print(genome, ','.join([str(x+1) for x in indexes[i]]), sep='\t', file=f)
    write_gff(
        gene_seqs,
        annotations,
        indexes[i],
        genome,
        options.outprefix + '_' + str(i+1) + '.gff',
        options.spacing_ns
    )


pyfastaq.utils.close(f)

genes_used = set()
for l in indexes:
    genes_used.update(l)

genes_used = sorted(list(genes_used))

f = pyfastaq.utils.open_file_write(options.outprefix + '_genes.fa')
for i in genes_used:
    print(gene_seqs[i], file=f)
