#!/usr/bin/env python3

import os
import argparse
import fastaq
import iva
import sys
import pysam
import copy
import subprocess


def syscall(cmd, allow_fail=False, verbose=False):
    if verbose:
        print('syscall:', cmd, flush=True)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        print('The following command failed with exit code', error.returncode, file=sys.stderr)
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(errors, file=sys.stderr)

        if allow_fail:
            return False, errors
        else:
            sys.exit(1)

    return True, None


class Clusters:
    def __init__(self, db_fasta, reads_1, reads_2, outdir, threads=1):
        self.db_fasta = os.path.abspath(db_fasta)
        self.reads_1 = os.path.abspath(reads_1)
        self.reads_2 = os.path.abspath(reads_2)
        self.outdir = os.path.abspath(outdir)
        self.bam_prefix = os.path.join(self.outdir, 'map_all_reads')
        self.bam = self.bam_prefix + '.bam'
        self.reads_by_clusters = {}
        self.clusters = []
        self.threads = threads
        self.catted_assemblies_fa = os.path.join(self.outdir, 'all_assemblies.fa')
        self.catted_genes_fa = os.path.join(self.outdir, 'all_genes.fa')
       

        try:
            os.mkdir(self.outdir)
        except:
            print('Error mkdir', self.outdir, file=sys.stderr)
            sys.exit(1)


    def _map_reads(self):
        iva.mapping.map_reads(
            self.reads_1,
            self.reads_2,
            self.db_fasta,
            self.bam_prefix, 
            index_k=20,
            index_s=5,
            minid=0.75,
            sort=False,
            threads = self.threads,
        )


    def _sam_to_fastq(self, s):
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


    def _get_reads_per_cluster(self):
        assert os.path.exists(self.bam)
        sam_reader = pysam.Samfile(self.bam, "rb")
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
                if ref not in self.reads_by_clusters:
                    self.reads_by_clusters[ref] = []

                read1 = self._sam_to_fastq(sam1)
                read2 = self._sam_to_fastq(s)
                if read1.id.endswith('/2'):
                    read1, read2 = read2, read1
                self.reads_by_clusters[ref].append((read1, read2))
                       
            sam1 = None 


    def _make_clusters(self):
        for ref in self.reads_by_clusters:
            directory = os.path.join(self.outdir, ref)
            cluster = Cluster(directory, ref, self.db_fasta, self.reads_by_clusters[ref], threads=self.threads)
            self.clusters.append(cluster)
            

    def _run_clusters(self):
        for cluster in self.clusters:
            cluster.run()


    def _make_act_files(self):
        files_to_cat = [x.final_assembly_fa for x in self.clusters if os.path.exists(x.final_assembly_fa)]
        if len(files_to_cat) > 0:
            fastaq.utils.syscall('cat ' + ' '.join(files_to_cat) + ' > ' + self.catted_assemblies_fa)

        files_to_cat = [x.gene_fa for x in self.clusters if os.path.exists(x.gene_fa)]
        if len(files_to_cat) > 0:
            fastaq.utils.syscall('cat ' + ' '.join(files_to_cat) + ' > ' + self.catted_genes_fa)


    def run(self):
        self._map_reads()
        self._get_reads_per_cluster()
        self._make_clusters()
        self._run_clusters()
        self._make_act_files()



class Cluster:
    def __init__(self, root_dir, gene_name, ref_fa, reads, threads=1):
        self.root_dir = os.path.abspath(root_dir)
        try:
            os.mkdir(self.root_dir)
        except:
            print('Error mkdir', self.root_dir, file=sys.stderr)
            sys.exit(1)

        self.gene_name = gene_name
        self.ref_fa_all = os.path.abspath(ref_fa)
        self.gene_fa = os.path.join(self.root_dir, 'gene.fa')
        self.reads = reads
        self.threads = threads
        self.reads1 = os.path.join(self.root_dir, 'reads_1.fq')
        self.reads2 = os.path.join(self.root_dir, 'reads_2.fq')
        self._write_gene_fasta()
        self._write_reads_fastqs()
        self.assembly_dir = os.path.join(self.root_dir, 'Assembly')
        self.assembly_scaffolds = os.path.join(self.root_dir, 'scaffolds.spades.fa')
        self.final_scaffolds = os.path.join(self.root_dir, 'scaffolds.final.fa')
        self.final_scaffolds_gapfilled = os.path.join(self.root_dir, 'scaffolds.final.gapfilled.fa')
        self.assembly_bam = os.path.join(self.root_dir, 'scaffolds.final.gapfilled.bam')
        self.gene_bam = os.path.join(self.root_dir, 'gene.bam')
        self.nucmer_coords = os.path.join(self.root_dir, 'nucmer.coords')
        self.final_assembly_fa = os.path.join(self.root_dir, 'assembly.fa')


    def _write_gene_fasta(self):
        if not os.path.exists(self.ref_fa_all + '.fai'):
            iva.common.syscall('samtools faidx ' + self.ref_fa_all)

        iva.common.syscall(' '.join([
            'samtools faidx',
            self.ref_fa_all,
            self.gene_name,
            '>', self.gene_fa
        ]))


    def _write_reads_fastqs(self):
        f1 = fastaq.utils.open_file_write(self.reads1)
        f2 = fastaq.utils.open_file_write(self.reads2)
        for t in self.reads:
            print(t[0], file=f1) 
            print(t[1], file=f2) 
        fastaq.utils.close(f1)
        fastaq.utils.close(f2)


    def _assemble(self):
        #if 'strA' not in self.gene_name:
        #    iva.common.syscall('touch ' + self.assembly_scaffolds)
        #    return

        cmd = ' '.join([
            '~mh12/bin/SPAdes-3.1.1-Linux/bin/spades.py',
            '--phred-offset 33',
            '-1', self.reads1,
            '-2', self.reads2,
            '-o', self.assembly_dir,
            '--threads', str(self.threads),
        ])
        self.assembled_ok, err = syscall(cmd, verbose=True, allow_fail=True)
        cwd = os.getcwd()
        os.chdir(self.root_dir)
        if self.assembled_ok:
            os.symlink(os.path.join('Assembly', 'scaffolds.fasta'), os.path.basename(self.assembly_scaffolds))
        else:
            with open('spades_errors', 'w') as f:
                print(err, file=f)
            f.close()
            
        os.chdir(cwd)

        
    def _scaffold_sopra(self):
        #if 'strA' not in self.gene_name:
        #    iva.common.syscall('touch ' + self.final_scaffolds)
        #    return

        cmd = ' '.join([
            '~mh12/git/perl/scaffold-wrapper-sopra.pl',
            self.assembly_scaffolds,
            'bowtie2',
            self.reads1,
            self.reads2,
            '10000',
            '500',
            os.path.join(self.root_dir, 'SOPRA')
        ])
        iva.common.syscall(cmd, verbose=True)
        cwd = os.getcwd()
        os.chdir(self.root_dir)
        os.symlink(os.path.join('SOPRA', 'scaffolds.fa'), os.path.basename(self.final_scaffolds))
        os.chdir(cwd)


    def _scaffold(self):
        #if 'strA' not in self.gene_name:
        #    iva.common.syscall('touch ' + self.final_scaffolds)
        #    return

        cmd = ' '.join([
            '~mh12/git/perl/scaffold-wrapper-sspace.pl',
            self.assembly_scaffolds,
            self.reads1,
            self.reads2,
            os.path.join(self.root_dir, 'SSPACE'),
            '350',
            '0.2',
        ])
        iva.common.syscall(cmd, verbose=True)
        cwd = os.getcwd()
        os.chdir(self.root_dir)
        os.symlink(os.path.join('SSPACE', 'scaffolds.fa'), os.path.basename(self.final_scaffolds))
        os.chdir(cwd)
        

    def _has_gaps(self, filename):
        seq_reader = fastaq.sequences.file_reader(filename) 
        for seq in seq_reader:
            if 'n' in seq.seq or 'N' in seq.seq:
                return True
        return False


    def _gap_fill(self):
        #if 'strA' not in self.gene_name:
        #    iva.common.syscall('touch ' + self.final_scaffolds)
        #    return
        cwd = os.getcwd()
        os.chdir(self.root_dir)

        if self._has_gaps(self.final_scaffolds):
            cmd = ' '.join([
                '~mh12/git/perl/gapFiller_wrapper.pl',
                self.final_scaffolds,
                os.path.join(self.root_dir, 'GapFiller'),
                self.reads1,
                self.reads2,
                '500',
                '0.4'
            ])
            iva.common.syscall(cmd, verbose=True)
            to_rename = os.path.join('GapFiller', 'standard_output.gapfilled.final.fa')
        else:
            to_rename = self.final_scaffolds
           
        seq_reader = fastaq.sequences.file_reader(to_rename)
        f = fastaq.utils.open_file_write(self.final_scaffolds_gapfilled)
        counter = 1
        for seq in seq_reader:
            seq.id = self.gene_name + '.' + str(counter)
            counter += 1
            print(seq, file=f)
        fastaq.utils.close(f)

        os.chdir(cwd)


    def _fix_orientation(self):
        if os.path.exists(self.final_scaffolds_gapfilled):
            f = fastaq.utils.open_file_write(self.final_assembly_fa)
            seq_reader = fastaq.sequences.file_reader(self.final_scaffolds_gapfilled)
            for seq in seq_reader:
                orfs = seq.all_orfs(min_length=50)
                if len(orfs):
                    orfs.sort(key=lambda t:len(t[0]))
                    if orfs[-1][1]:
                        seq.revcomp()
                print(seq, file=f)                

            fastaq.utils.close(f)


    def _map_reads(self, ref, out):
        iva.mapping.map_reads(
            self.reads1,
            self.reads2,
            ref,
            out[:-4], 
            index_k=9,
            index_s=2,
            minid=0.9,
            sort=True,
            threads = self.threads,
        )
        os.unlink(out[:-4] + '.unsorted.bam')
        syscall('samtools index ' + out)


    def _map_reads_to_gene(self):
        self._map_reads(self.gene_fa, self.gene_bam)

    def _map_reads_to_assembly(self):
        print('\n\nmap reads. self.final_scaffolds_gapfilled=', self.final_scaffolds_gapfilled, '    self.assembly_bam=', self.assembly_bam, '\n\n')
        self._map_reads(self.final_scaffolds_gapfilled,  self.assembly_bam)


    def _nucmer_assembly_v_gene(self):
        iva.mummer.run_nucmer(
            self.final_assembly_fa,
            self.gene_fa,
            self.nucmer_coords,
            min_id=95,
            min_length=50,
            breaklen=50
        )


    def _compare_assembly_with_gene(self):
        # load nucmer file into memory
        # work out how assembly compares to gene
        pass



    def run(self):
        self._map_reads_to_gene()
        self._assemble()
        if self.assembled_ok:
            self._scaffold()
            self._gap_fill()
            self._fix_orientation()
            self._map_reads_to_assembly()
            self._nucmer_assembly_v_gene()
            self._compare_assembly_with_gene()


parser = argparse.ArgumentParser(
    description = 'Test resistome finder script',
    usage = '%(prog)s [options] <db.fa> <reads1.fq> <reads2.fq> <outdir>')
parser.add_argument('db', help='FASTA file of reference genes')
parser.add_argument('reads_1', help='Name of fwd reads fastq file')
parser.add_argument('reads_2', help='Name of rev reads fastq file')
parser.add_argument('outdir', help='Output directory')
parser.add_argument('--threads', type=int, help='Number of threads for smalt and spades [%(default)s]', default=1)
options = parser.parse_args()

clusters = Clusters(
    options.db,
    options.reads_1,
    options.reads_2,
    options.outdir,
    threads=options.threads
)

clusters.run()
