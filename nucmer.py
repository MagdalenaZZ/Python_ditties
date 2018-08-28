import utils

class Error (Exception): pass

def file_reader(fname):
    f = utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    utils.close(f)



class NucmerHit:
    def __init__(self, line):
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]

        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0])
            self.ref_end = int(l[1])
            self.qry_start = int(l[2])
            self.qry_end = int(l[3])
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.strand = int(l[10])
            self.ref_name = l[11]
            self.qry_name = l[12]

            if len(l) == 14:
                self.tag = l[13][1:-1]
            else:
                self.tag = None
        except:
            raise Error('Error reading this nucmer line:\n' + line)


    def __str__(self):
        s = '\t'.join(str(x) for x in
            [self.ref_start,
            self.ref_end,
            self.qry_start,
            self.qry_end,
            self.hit_length_ref,
            self.hit_length_qry,
            '{0:.2f}'.format(self.percent_identity),
            self.ref_length,
            self.qry_length,
            self.frame,
            self.strand,
            self.ref_name,
            self.qry_name])

        if self.tag is not None:
            s += '\t[' + self.tag + ']'

        return s

