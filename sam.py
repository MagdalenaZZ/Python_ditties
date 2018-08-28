import cigar
import utils
import fastn

class Error (Exception): pass

def file_reader(fname):
    f = utils.open_file_read(fname)
    for line in f:
        if line.startswith('@'):
            continue

        yield SamRecord(line)

    utils.close(f)



class SamRecord:
    def __init__(self, line):
        # example line:
        # HS4_6280:2:1104:12102:124607  99  PyYM_01_v1  1   47  2S73M   =   362 438 TGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT HHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF AS:i:73
        try:
            (self.id, self.flag, self.rname, self.pos, self.mapq, self.cigar,
             self.mrname, self.mpos, self.isize, self.seq, self.qual, *self.tags_list) = line.rstrip().split('\t')

            self.pos = int(self.pos) - 1
            self.flag = int(self.flag)
            self.mapq = int(self.mapq)
            self.cigar = cigar.Cigar(self.cigar)
            self.mpos = int(self.mpos) - 1
            self.isize = int(self.isize)
            self.tags = {}
            for tag in self.tags_list:
                (tag, type, value) = tag.split(':', 2)
                if type == 'i':
                    value = int(value)
                self.tags[tag] = (type, value)
        except:
            raise Error('Error reading this sam line:\n' + line)


    def __str__(self):
        return '\t'.join([str(x) for x in
            [self.id, self.flag, self.rname, self.pos+1, self.mapq, self.cigar,
             self.mrname, self.mpos+1, self.isize, self.seq, self.qual] + self.tags_list])

    def get_differences_from_ref(self, ref):
        if self.is_mapped():
            return self.cigar.get_differences_from_ref(self.seq, ref, self.pos)
        else:
            return []

    def is_paired(self):
        return self.flag & 0x0001 != 0

    def is_proper_pair(self):
        return self.flag & 0x002 != 0

    def is_mapped(self):
        return self.flag & 0x0004 == 0

    def is_mate_mapped(self):
        return self.flag & 0x0008 == 0

    def is_forward_strand(self):
        return self.flag & 0x0010 == 0

    def is_first_of_pair(self):
        return self.flag & 0x0040

    def is_second_of_pair(self):
        return self.flag & 0x0080

    def is_duplicate(self):
        return self.flag & 0x0400 != 0

    def query_strand(self):
        if self.flag & 0x0010:
            return '-'
        else:
            return '+'

    def ref_hit_end_position(self):
        if self.is_mapped():
            return self.pos + self.cigar.ref_hit_length() - 1
        else:
            raise Error('Cannot get end of hit position for unmapped read!\n' + str(self))

    def to_fastn(self):
        if self.qual == '*':
            seq = fastn.Fasta(self.id, self.seq)
        else:
            seq = fastn.Fastq(self.id, self.seq, self.qual)

        if self.query_strand() == '-':
            seq.revcomp()

        if self.is_first_of_pair():
            seq.id += '/1'
        elif self.is_second_of_pair():
            seq.id += '/2'

        return seq




def get_sequence_lengths(fname):
    lengths = {}
    f = utils.open_file_read(fname)
    for line in f:
        if not line.startswith('@'):
            break
        elif line.startswith('@SQ'):
            try:
                l = line.rstrip().split('\t')[1:]
                d = {x[:2]:x[3:] for x in l}
                lengths[d['SN']] = int(d['LN'])
            except:
                raise Error('Error getting length from line of BAM header\n' + line)

    utils.close(f)
    return lengths

