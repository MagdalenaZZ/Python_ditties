import sys
import os
import re
import string
import copy

import utils
import genome_intervals

class Error (Exception): pass


# python 3's seek is glacially slow.  When we read a fasta file, we know
# we've reached the end of a sequence when we get a new line starting with
# '>'.  Instead of using seek and tell, we just remember the previous line
# of the file, for any given filehandle
previous_lines = {}


codon2aa = {
'GCA': 'A',
'GCC': 'A',
'GCG': 'A',
'GCT': 'A',
'AGA': 'R',
'AGG': 'R',
'CGA': 'R',
'CGC': 'R',
'CGG': 'R',
'CGT': 'R',
'AAC': 'N',
'AAT': 'N',
'GAC': 'D',
'GAT': 'D',
'TGC': 'C',
'TGT': 'C',
'GAA': 'E',
'GAG': 'E',
'CAA': 'Q',
'CAG': 'Q',
'GGA': 'G',
'GGC': 'G',
'GGG': 'G',
'GGT': 'G',
'CAC': 'H',
'CAT': 'H',
'ATA': 'I',
'ATC': 'I',
'ATT': 'I',
'TTA': 'L',
'TTG': 'L',
'CTA': 'L',
'CTC': 'L',
'CTG': 'L',
'CTT': 'L',
'AAA': 'K',
'AAG': 'K',
'ATG': 'M',
'TTC': 'F',
'TTT': 'F',
'CCA': 'P',
'CCC': 'P',
'CCG': 'P',
'CCT': 'P',
'AGC': 'S',
'AGT': 'S',
'TCA': 'S',
'TCC': 'S',
'TCG': 'S',
'TCT': 'S',
'ACA': 'T',
'ACC': 'T',
'ACG': 'T',
'ACT': 'T',
'TGG': 'W',
'TAC': 'Y',
'TAT': 'Y',
'GTA': 'V',
'GTC': 'V',
'GTG': 'V',
'GTT': 'V',
'TAA': '*',
'TAG': '*',
'TGA': '*'}

def file_reader(fname, read_quals=False):
    f = utils.open_file_read(fname)
    line = f.readline()
    if line.startswith('>'):
        seq = Fasta()
        previous_lines[f] = line
    elif line.startswith('@'):
        seq = Fastq()
        previous_lines[f] = line
    elif line == '':
        utils.close(f)
        return
    else:
        raise Error('Error determining file type from file "' + fname + '".  First line is:\n' + line.rstrip())

    while seq.get_next_from_file(f, read_quals):
        yield seq

    utils.close(f)

class Fasta:
    # this defines the line length when printing sequences
    line_length = 60

    def _get_id_from_header_line(self, line):
        if line.startswith('>'):
            return line.rstrip()[1:]
        else:
            raise Error('Error! expected line starting with ">", but got this:\n', line)


    def __init__(self, id_in=None, seq_in=None):
        self.id = id_in
        self.seq = seq_in

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.seq)

    def split_capillary_id(self):
        try:
            a = self.id.rsplit('.', 1)
            if a[1].startswith('p'):
                dir = 'fwd'
            elif a[1].startswith('q'):
                dir = 'rev'
            else:
                dir = 'unk'

            return {'prefix': a[0], 'dir': dir, 'suffix':a[1]}
        except:
            raise Error('Error in split_capillary_id() on ID', self.id)

    def strip_illumina_suffix(self):
        if self.id.endswith('/1') or self.id.endswith('/2'):
            self.id = self.id[:-2]

    def revcomp(self):
        self.seq = self.seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]


    def trim_Ns(self):
        self.seq = self.seq.strip('Nn')


    def replace_bases(self, old, new):
        self.seq = self.seq.replace(old, new)

    # returns list of gaps in the sequence (zero based) that are
    # at least min_length long
    def gaps(self, min_length = 1):
        gaps = []
        regex = re.compile('N+', re.IGNORECASE)
        for m in regex.finditer(self.seq):
             if m.span()[1] - m.span()[0] + 1 >= min_length:
                 gaps.append(genome_intervals.Interval(m.span()[0], m.span()[1] - 1))
        return gaps

    # Returns list of coords (zero based) of the contigs within the
    # sequence. i.e. everthing that's not a gap.
    def contig_coords(self):
        # contigs are the opposite of gaps, so work out the coords from the gap coords
        gaps = self.gaps()

        if len(gaps) == 0:
            return [genome_intervals.Interval(0, len(self) - 1)]

        coords = [0]
        for g in gaps:
            if g.start == 0:
                coords = [g.end + 1]
            else:
                coords += [g.start - 1, g.end + 1]

        if coords[-1] + 1 < len(self):
            coords.append(len(self) - 1)

        return [genome_intervals.Interval(coords[i], coords[i+1]) for i in range(0, len(coords)-1,2)]



    # Fills the object with the next sequence in the file.  Returns
    # True if this was successful, False if no more sequences in the file.
    # If reading a file of quality scores, set read_quals = True
    def get_next_from_file(self, f, read_quals=False):
        if f in previous_lines:
            if previous_lines[f] == None:
                self.id = self.seq = None
                return False
            else:
                self.id = self._get_id_from_header_line(previous_lines[f])
        else:
            line = '\n'
            while line == '\n':
                line = f.readline()
            self.id = self._get_id_from_header_line(line)

        self.seq = ''
        seq_lines = [] # much faster to store the seq lines in an array,
                       # then join at the end

        while 1:
            line = f.readline()

            if line.startswith('>'):
                previous_lines[f] = line.rstrip()
                break
            elif line == '':
                previous_lines[f] = None
                break
            else:
                 seq_lines.append(line.rstrip())

        if read_quals:
            self.seq = ' '.join(seq_lines)
        else:
            self.seq = ''.join(seq_lines)
        return True

    def __str__(self):
        if Fasta.line_length == 0:
            return '>' + self.id + '\n' + self.seq
        else:
            return '>' + self.id + '\n' + '\n'.join(self.seq[i:i+Fasta.line_length] for i in range(0, len(self), Fasta.line_length))

    def __getitem__(self, index):
        return self.seq[index]

    def trim(self, start, end):
        self.seq = self.seq[start:len(self.seq) - end]

    # qual_scores should be a list of quality scores
    def to_Fastq(self, qual_scores):
        if len(self) != len(qual_scores):
            raise Error('Error making Fastq from Fasta, lengths differ.', self.id)
        return Fastq(self.id, self.seq, ''.join([chr(max(0, min(x, 93)) + 33) for x in qual_scores]))

    def search(self, search_string):
        seq = self.seq.upper()
        search_string = search_string.upper()
        pos = 0
        found = seq.find(search_string, pos)
        hits = []

        while found != -1:
            hits.append((found, '+'))
            pos = found + 1
            found = seq.find(search_string, pos)


        pos = 0
        search_string = Fasta('x', search_string)
        search_string.revcomp()
        search_string = search_string.seq
        found = seq.find(search_string, pos)

        while found != -1:
            hits.append((found, '-'))
            pos = found + 1
            found = seq.find(search_string, pos)

        return hits

    def translate(self):
        return Fasta(self.id, ''.join([codon2aa.get(self.seq[x:x+3].upper(), 'X') for x in range(0, len(self)-1, 3)]))

class Fastq(Fasta):
    def __init__(self, id_in=None, seq_in=None, qual_in=None):
        super().__init__(id_in, seq_in)
        self.qual = qual_in
        if (not self.seq == self.qual == None) and len(self.qual) != len(self.seq):
            raise Error('Error constructing Fastq.  Mismatch in sequence and quality length\n' + str(self))


    def __str__(self):
        return '@' + self.id + '\n' + self.seq + '\n+\n' + self.qual
        return '>' + self.id + '\n' + self.seq

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_next_from_file(self, f, read_quals=False):
        if f in previous_lines:
            line = previous_lines[f]
            del previous_lines[f]
        else:
            line = f.readline()

        while line == '\n':
            line = f.readline()

        if not line:
            self = Fastq('', '', '')
            return False

        if not line.startswith('@'):
            raise Error('Error getting next sequence from fastq file.  Got line:\n' + line)

        self.id = line.rstrip()[1:]
        line = f.readline()
        if not line:
            raise Error('Error getting next sequence from fastq file, sequence has ID ' + self.id)

        self.seq = line.strip()

        line = f.readline()
        if not (line and line.startswith('+')):
            raise Error('Error getting next sequence from fastq file, no line starting with +,  sequence has ID ' + self.id)

        line = f.readline()
        if not line:
            raise Error('Error getting next sequence from fastq file, sequence has ID ' + self.id)

        self.qual = line.rstrip()
        return True

    def revcomp(self):
        super().revcomp()
        self.qual = self.qual[::-1]

    def trim(self, start, end):
        super().trim(start, end)
        self.qual = self.qual[start:len(self.qual) - end]

    def to_Fasta_and_qual(self):
        quals = [ord(x) - 33 for x in self.qual]
        return (Fasta(self.id, self.seq), quals)


    def trim_Ns(self):
        # get index of first base that is not an N
        i = 0
        while i < len(self) and self.seq[i] in 'nN':
            i += 1

        # strip off start of sequence and quality
        self.seq = self.seq[i:]
        self.qual = self.qual[i:]

        # strip the ends
        self.seq = self.seq.rstrip('Nn')
        self.qual = self.qual[:len(self.seq)]

    def translate(self):
        fa = super().translate()
        return Fastq(fa.id, fa.seq, 'I'*len(fa.seq))

def deinterleave(infile, outfile_1, outfile_2, fasta_out=False):
    seq_reader = file_reader(infile)
    f_1 = utils.open_file_write(outfile_1)
    f_2 = utils.open_file_write(outfile_2)
    for seq in seq_reader:
        if fasta_out:
            print(Fasta(seq.id, seq.seq), file=f_1)
        else:
            print(seq, file=f_1)
        try:
            next(seq_reader)
        except StopIteration:
            utils.close(f_1)
            utils.close(f_2)
            raise Error('Error getting mate for sequence. Cannot continue')
        if fasta_out:
            print(Fasta(seq.id, seq.seq), file=f_2)
        else:
            print(seq, file=f_2)

    utils.close(f_1)
    utils.close(f_2)

def interleave(infile_1, infile_2, outfile):
    seq_reader_1 = file_reader(infile_1)
    seq_reader_2 = file_reader(infile_2)
    f_out = utils.open_file_write(outfile)

    for seq_1 in seq_reader_1:
        try:
            seq_2 = next(seq_reader_2)
        except:
            utils.close(f_out)
            raise Error('Error getting mate for sequence', seq_1.id, ' ... cannot continue')

        print(seq_1, file=f_out)
        print(seq_2, file=f_out)

    try:
        seq_2 = next(seq_reader_2)
    except:
        seq_2 = None

    if seq_2 is not None:
        utils.close(f_out)
        raise Error('Error getting mate for sequence', seq_2.id, ' ... cannot continue')

    utils.close(f_out)

def reverse_complement(infile, outfile):
    seq_reader = file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.revcomp()
        print(seq, file=fout)

    utils.close(fout)

def trim(infile, outfile, start, end):
    seq_reader = file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.trim(start, end)
        if len(seq):
            print(seq, file=fout)

    utils.close(fout)

def fastq_to_mira_xml(infile, outfile):
    seq_reader = file_reader(infile)
    fout = utils.open_file_write(outfile)
    print('<?xml version="1.0"?>', '<trace_volume>', sep='\n', file=fout)

    for seq in seq_reader:
        print('    <trace>',
              '        <trace_name>' + seq.id + '</trace_name>',
              '        <clip_quality_right>' + str(len(seq)) + '</clip_quality_right>',
              '        <clip_vector_left>1</clip_vector_left>',
              '    </trace>', sep='\n', file=fout)


    print('</trace_volume>', file=fout)
    utils.close(fout)

def file_to_dict(infile, d):
    seq_reader = file_reader(infile)
    for seq in seq_reader:
        d[seq.id] = copy.copy(seq)

def lengths_from_fai(fai_file, d):
    f = utils.open_file_read(fai_file)
    for line in f:
        (id, length) = line.rstrip().split()[:2]
        d[id] = int(length)
    utils.close(f)


def split_by_base_count(infile, outfiles_prefix, max_bases, max_seqs=None):
    '''Splits a fasta/q file into separate files, file size determined by number of bases.

    Puts <= max_bases in each split file The eException is a single sequence >=max_bases
    is put in its own file.  This does not split sequences.
    '''
    seq_reader = file_reader(infile)
    base_count = 0
    file_count = 1
    seq_count = 0
    fout = None
    if max_seqs is None:
        max_seqs = float('inf')

    for seq in seq_reader:
        if base_count == 0:
            fout = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
            file_count += 1

        if base_count + len(seq) > max_bases or seq_count >= max_seqs:
            if base_count == 0:
                print(seq, file=fout)
                utils.close(fout)
            else:
                utils.close(fout)
                fout = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
                print(seq, file=fout)
                base_count = len(seq)
                file_count += 1
                seq_count = 1
        else:
            base_count += len(seq)
            seq_count += 1
            print(seq, file=fout)

    utils.close(fout)

def count_sequences(infile):
    '''Returns the number of sequences in a file'''
    seq_reader = file_reader(infile)
    n = 0
    for seq in seq_reader:
        n += 1
    return n


def fasta_to_fastq(fasta_in, qual_in, outfile):
    fa_reader = file_reader(fasta_in)
    qual_reader = file_reader(qual_in, read_quals=True)
    f_out = utils.open_file_write(outfile)

    for seq in fa_reader:
        qual = next(qual_reader)
        if seq.id != qual.id:
            raise Error('Mismatch in names from fasta and qual file', seq.id, qual.id)

        qual.seq = [int(x) for x in qual.seq.split()]
        print(seq.to_Fastq(qual.seq), file=f_out)

    utils.close(f_out)

def fastn_to_quasr_primers(infile, outfile):
    seq_reader = file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq2 = copy.copy(seq)
        seq2.revcomp()
        print(seq.seq, seq2.seq, sep='\t', file=f_out)

    utils.close(f_out)

def replace_bases(infile, outfile, old, new):
    seq_reader = file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.replace_bases(old, new)
        print(seq, file=f_out)

    utils.close(f_out)

