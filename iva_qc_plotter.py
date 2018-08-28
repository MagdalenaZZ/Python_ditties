#!/usr/bin/env python3

import sys
import os
import argparse
import fastaq
import iva

qc_stats = [
    'ref_EMBL_dir',
    'ref_EMBL_files',
    'ref_bases',
    'ref_sequences',
    'ref_bases_assembled',
    'ref_sequences_assembled',
    'ref_sequences_assembled_ok',
    'ref_bases_assembler_missed',
    'assembly_bases',
    'assembly_bases_in_ref',
    'assembly_contigs',
    'assembly_contigs_hit_ref',
    'assembly_bases_reads_disagree',
    'cds_number',
    'cds_assembled',
    'cds_assembled_ok',
    'gage_Missing_Reference_Bases',
    'gage_Missing_Assembly_Bases',
    'gage_Missing_Assembly_Contigs',
    'gage_Duplicated_Reference_Bases',
    'gage_Compressed_Reference_Bases',
    'gage_Bad_Trim',
    'gage_Avg_Idy',
    'gage_SNPs',
    'gage_Indels_<_5bp',
    'gage_Indels_>=_5',
    'gage_Inversions',
    'gage_Relocation',
    'gage_Translocation',
    'ratt_elements_found',
    'ratt_elements_transferred',
    'ratt_elements_transferred_partially',
    'ratt_elements_split',
    'ratt_parts_of_elements_not_transferred',
    'ratt_elements_not_transferred',
    'ratt_gene_models_to_transfer',
    'ratt_gene_models_transferred',
    'ratt_gene_models_transferred_partially',
    'ratt_exons_not_transferred_from_partial_matches',
    'ratt_gene_models_not_transferred',
    'longest_matching_contig',
]


class Assembly:
    def __init__(self, l):
        self.id, self.assembler, self.contigs_file, self.qc_prefix = l.rstrip().split()
        self.qc_stats_file = self.qc_prefix + '.stats.txt'
        self.nucmer_coords_file = self.qc_prefix + '.assembly_vs_ref.coords'
        self.qc = {x: -1 for x in qc_stats}
        self._load_qc_stats()


    def _load_qc_stats(self):
        f = fastaq.utils.open_file_read(self.qc_stats_file)
        for line in f:
            key, value = line.rstrip().split('\t')
            assert key in self.qc
            if key.startswith('ref_EMBL'):
                self.qc[key]=value
            elif '.' in value:
                self.qc[key] = float(value)
            elif value.isdigit():
                self.qc[key] = int(value)
                
        fastaq.utils.close(f)

        self.qc['longest_matching_contig'] = self._get_longest_matching_contig_length()


    def _get_longest_matching_contig_length(self):
        max_length = -1
        f = fastaq.utils.open_file_read(self.nucmer_coords_file)
        for line in f:
             if len(line) and line[0].isdigit():
                 max_length = max(max_length, int(line.split('\t')[8]))
        fastaq.utils.close(f)
        return max_length
                

    def to_R_data_line(self):
        return '\t'.join([
            self.id,
            self.assembler,
            '\t'.join([str(self.qc[x]) for x in qc_stats]) 
        ])


class Assemblies:
    def __init__(self, input_file, outprefix):
        self.outprefix = outprefix
        self.r_data_file = self.outprefix + '.dat'
        self.r_script = self.outprefix + '.R'
        self._read_input_file(input_file)
        self._get_assemblers() 
        self._filter_assemblies()
        

    def _read_input_file(self, filename):
        self.assemblies = {}
        f = fastaq.utils.open_file_read(filename)
        for line in f:
            a = Assembly(line)
            if a.id not in self.assemblies:
                self.assemblies[a.id] = []
            self.assemblies[a.id].append(a)
        fastaq.utils.close(f)


    def _get_assemblers(self):
        all_assembler_tuples = set()
        for assembly_id in self.assemblies:
            assemblers = tuple(sorted([x.assembler for x in self.assemblies[assembly_id]]))
            all_assembler_tuples.add(assemblers)

        biggest_tuple_size = max([len(x) for x in all_assembler_tuples])
        biggest_tuples = [x for x in all_assembler_tuples if len(x) == biggest_tuple_size]
        number_of_biggest = len(biggest_tuples)
        if number_of_biggest != 1:
            print('Error! conflicting assembler names. Cannot continue', file=sys.stderr)

        self.assemblers = set(biggest_tuples[0])


    def _filter_assemblies(self):
        to_delete = set()
        for assembly_id in self.assemblies:
            assemblers = set([x.assembler for x in self.assemblies[assembly_id]])
            if assemblers != self.assemblers:
                print('Ignoring ID', assembly_id, 'because not all assemblers present. Needed:', self.assemblers, ' got:', assemblers)
                to_delete.add(assembly_id)
        
        for assembly_id in to_delete:
            del self.assemblies[assembly_id]

        print(len(self.assemblies), 'samples remaining after filtering')


    def _write_R_data_file(self):
        f = fastaq.utils.open_file_write(self.r_data_file)
        print('ID', 'assembler', '\t'.join(qc_stats), sep='\t', file=f)
        for assembly_id in self.assemblies:
            for assembly in self.assemblies[assembly_id]:
                print(assembly.to_R_data_line(), file=f)
        fastaq.utils.close(f) 


    def _write_R_script(self):
        f = fastaq.utils.open_file_write(self.r_script)
        print('library(ggplot2)',
              'd = read.csv(file="' + self.r_data_file + r'''", sep="\t", header=T)''',
              'pre="' + self.outprefix + '"',
              r'''              
ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=longest_matching_contig)) +
    geom_hline(aes(yintercept=9000), color="blue", linetype="dashed") +
    xlab("Assembler") +
    ylab("Length of longest contig that matches reference")
ggsave(filename=paste(pre, "plot.longest_matching_contig.boxplot.pdf", sep=""))


ggplot() + geom_freqpoly(data=d, aes(x=longest_matching_contig, colour=assembler)) +
    xlab("Length of longest contig that matches reference") +
    ylab("Count")
ggsave(filename=paste(pre, "plot.longest_matching_contig.hist.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=100*ref_bases_assembler_missed/ref_bases)) +
    xlab("Assembler") +
    ylab("Per cent of reference missed, that had good read coverage")
ggsave(filename=paste(pre, "plot.percent_should_have_assembled.pdf",  sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler,y=100*gage_Missing_Reference_Bases/ref_bases)) +
    xlab("Assembler") +
    ylab("Per cent reference bases not assembled according to GAGE")
ggsave(filename=paste(pre, "plot.gage_percent_unassembled.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=gage_Duplicated_Reference_Bases/ref_bases)) +
    xlab("Assembler") +
    ylab("GAGE duplication rate")
ggsave(filename=paste(pre, "plot.gage_duplication_rate.pdf",  sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=gage_Avg_Idy)) +
    xlab("Assembler") +
    ylab("GAGE average identity between reference and assembly")
ggsave(filename=paste(pre, "plot.gage_average_identity.pdf",  sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=gage_Compressed_Reference_Bases/ref_bases)) +
    xlab("Assembler") +
    ylab("GAGE compressed reference bases")
ggsave(filename=paste(pre, "plot.gage_percent_compressed.pdf",  sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler,y=gage_Inversions+gage_Relocation+gage_Translocation)) +
    xlab("Assembler") +
    ylab("GAGE assembly errors")
ggsave(filename=paste(pre, "plot.gage_errors.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler,y=100*ratt_elements_transferred/ratt_elements_found)) +
    xlab("Assembler") +
    ylab("Per cent elements transferred by RATT")
ggsave(filename=paste(pre, "plot.ratt_elements_transferred.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler,y=100*ratt_gene_models_transferred/ratt_gene_models_to_transfer)) +
    xlab("Assembler") +
    ylab("Per cent genes transferred by RATT")
ggsave(filename=paste(pre, "plot.ratt_genes_transferred.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=100*assembly_bases_reads_disagree/assembly_bases)) +
    xlab("Assembler") +
    ylab("Per cent of assembly bases that disagree with mapped reads")
ggsave(filename=paste(pre,"plot.bases_disagree.pdf", sep=""))


ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=assembly_contigs)) +
    xlab("Assembler") +
    ylab("Number of contigs or scaffolds")
ggsave(filename=paste(pre,"plot.contig_count.pdf", sep=""))

ggplot() + geom_boxplot(data=d, mapping=aes(x=assembler, y=ref_sequences_assembled_ok)) +
    xlab("Assembler") +
    ylab("Number of segments assembled into 1 unique contig")
ggsave(filename=paste(pre,"plot.segments_assembled_well.pdf", sep=""))


''', sep='\n', file=f)
        fastaq.utils.close(f) 

    def make_plots(self):
        self._write_R_data_file()
        self._write_R_script()
        iva.common.syscall('R CMD BATCH ' + self.r_script)


parser = argparse.ArgumentParser(
    usage = '%(prog)s <infile> <outdir>')
parser.add_argument('infile', help='Name of input file')
parser.add_argument('outdir', help='Name of output directory')
options = parser.parse_args()



a = Assemblies(options.infile, os.path.join(options.outdir, 'out'))
os.mkdir(options.outdir)
a.make_plots()
