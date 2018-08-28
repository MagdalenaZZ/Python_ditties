import fastn
import utils

class Error (Exception): pass

def file_reader(fname):
    f = utils.open_file_read(fname)

    for line in f:
        yield BlastHit(line)

    utils.close(f)



class BlastHit:
    def __init__(self, line):
        try:
            l = line.rstrip().split('\t')
            self.qry_name = l[0]
            self.ref_name = l[1]
            self.percent_identity = float(l[2])
            self.alignment_length = int(l[3])
            self.mismatches = int(l[4])
            self.gap_openings = int(l[5])
            self.qry_start = int(l[6])
            self.qry_end = int(l[7])
            self.ref_start = int(l[8])
            self.ref_end = int(l[9])
            self.e_value = float(l[10])
            self.bit_score = float(l[11])

            if len(l) == 12:
                self.qry_length = None
                self.ref_length = None
            elif len(l) == 14:
                self.qry_length = int(l[12])
                self.ref_length = int(l[13])
            else:
                raise Error('Error reading this blast line:\n' + line)
        except:
            raise Error('Error reading this blast line:\n' + line)


    def __str__(self):
        s =  '\t'.join(str(x) for x in
            [self.qry_name,
             self.ref_name,
             '{:.2f}'.format(self.percent_identity),
             self.alignment_length,
             self.mismatches,
             self.gap_openings,
             self.qry_start,
             self.qry_end,
             self.ref_start,
             self.ref_end,
             self.e_value,
             '{:g}'.format(self.bit_score)])

        if self.qry_length is not None:
            s += '\t' + str(self.qry_length) + '\t' +  str(self.ref_length)

        return s

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def add_sequence_lengths(self, ref_lengths, qry_lengths):
        try:
            self.ref_length = ref_lengths[self.ref_name]
            self.qry_length = qry_lengths[self.qry_name]
        except:
            raise Error('Error adding lengths\n', self)


def add_sequence_lengths(infile, ref_fai, qry_fai, outfile):
    ref_lengths = {}
    qry_lengths = {}

    fastn.lengths_from_fai(ref_fai, ref_lengths)
    fastn.lengths_from_fai(qry_fai, qry_lengths)

    f = utils.open_file_write(outfile)
    blast_reader = file_reader(infile)
    for hit in blast_reader:
        hit.add_sequence_lengths(ref_lengths, qry_lengths)
        print(hit, file=f)
    utils.close(f)

