#!/usr/bin/env python3.3

import fastn

class Hit:
    def __init__(self, id, start, strand):
        self.id = id
        self.start = start
        self.strand = strand

    def __str__(self):
        return self.id + ':::' + str(self.start + 1) + ':::' + self.strand

    def __lt__(self, other):
        return self.id < other.id or (self.id == other.id and self.start < other.start)


class Tag:
    def __init__(self, id, start, end, seq):
        self.id = id
        self.start = start
        self.end = end
        self.seq = seq
        self.qry_hits = set()
        self.ref_hits = set()

    def __str__(self):
        def hits2string(s):
            if len(s):
                return ' '.join([str(x) for x in s])
            else:
                return '*'

        return '\t'.join([self.id,
                          str(self.start + 1),
                          str(self.end + 1),
                          self.seq,
                          hits2string(self.qry_hits),
                          hits2string(self.ref_hits)])

    def __lt__(self, other):
        return self.id < other.id or (self.id == other.id and self.start < other.start) or (self.id == other.id and self.start == other.start and self.end < other.end)

    def to_fasta(self):
        return fastn.Fasta(self.id, self.seq)

    def is_unique(self):
        return len(self.qry_hits) == len(self.ref_hits) == 1

    def is_unique_in_qry(self):
        return len(self.qry_hits) == 1

