import utils
import fastn

class Error (Exception): pass

def file_reader(fname):
    f = utils.open_file_read(fname)

    for line in f:
        yield MpileupLine(line)

    utils.close(f)


class MpileupLine:
    def __init__(self, line):
        # MAL1   3       N       3       AA^]A   HHH
        try:
            l = line.rstrip().split('\t')
            self.id = l[0]
            self.position = int(l[1])
            self.ref_base = l[2]
            self.depth = int(l[3])
            self.qry_bases = l[4]
            self.qry_quals = l[5]
        except:
            raise Error('Error reading this mpileup line:\n' + line)


    def __str__(self):
        s = '\t'.join(str(x) for x in
            [self.id,
            self.position,
            self.ref_base,
            self.depth,
            self.qry_bases,
            self.qry_quals])

        return s


def get_mean_coverage_per_seq(mpileup_file, fai_file):
    mean_cov_per_seq = {}
    seq_lengths = {}
    fastn.lengths_from_fai(fai_file, seq_lengths)
    depth_sum = 0
    current_id = None
    missed_seqs = set(seq_lengths.keys())

    reader = file_reader(mpileup_file)
    for mp in reader:
        if current_id != mp.id:
            if current_id is not None:
                mean_cov_per_seq[current_id] = 1.0 * depth_sum / seq_lengths[current_id]

            current_id = mp.id
            depth_sum = 0
            missed_seqs.remove(current_id)

        depth_sum += mp.depth

    mean_cov_per_seq[current_id] = 1.0 * depth_sum / seq_lengths[current_id]
    for id in missed_seqs:
        mean_cov_per_seq[id] = 0

    return mean_cov_per_seq

