#!/usr/bin/env python3.3

import fastn
import utils
import argparse

parser = argparse.ArgumentParser(
    description = 'Generates every possible sequence combination as FASTA. e.g. [ACT]GT would give AGT, CGT and TGT.',
    usage = '%(prog)s [options] <sequence> <outfile>')
parser.add_argument('--prefix', help='Prefix to add to evey output sequence [none]', default = '')
parser.add_argument('sequence', help='Sequence. Put bases with more than one possibility in square brackets.')
parser.add_argument('outfile', help='Name of fasta output file')
options = parser.parse_args()

class Seq:
    def __init__(self, s):
        self.bases = []
        in_list = False
        self.first_to_increment = None
        for b in s:
            if not in_list:
                if b == '[':
                    in_list = True
                    self.bases.append([])
                    if self.first_to_increment is None:
                        self.first_to_increment = len(self.bases) - 1
                else:
                    self.bases.append([b])
            else:
                if b == ']':
                    in_list = False
                else:
                    self.bases[-1].append(b)

        self.positions = [0] * len(self.bases)
        self.all_incremented = False

    def __repr__(self):
        return str(self.bases) + '\n' + str(self.positions)

    def to_fasta(self, id):
        seq = ''.join([self.bases[i][self.positions[i]] for i in range(len(self.positions))])
        return fastn.Fasta(id, seq)

    def increment(self):
        for i in reversed(range(len(self.positions))):
            if len(self.bases[i]) == 1:
                continue
            elif self.positions[i] + 1 < len(self.bases[i]):
                self.positions[i] += 1
                break
            else:
                 self.positions[i] = 0
                 if i == self.first_to_increment:
                     self.all_incremented = True
                     break

    def reset(self):
        self.positions = [0] * len(self.bases)
        self.all_incremented = False

    def all_combinations_found(self):
        return self.all_incremented


s = Seq(options.sequence)
f = utils.open_file_write(options.outfile)
i = 1

while not s.all_combinations_found():
    print(s.to_fasta(options.prefix + str(i)), file=f)
    i += 1
    s.increment()

utils.close(f)

