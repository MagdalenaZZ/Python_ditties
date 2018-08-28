import utils

# Class for reading GenomeDiff file. See
# http://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_format.html

class Error (Exception): pass

mutation_types = set(['SNP', 'SUB', 'DEL', 'INS', 'MOB', 'AMP', 'CON', 'INV'])

class Mutation:
    def __init__(self, s):
        fields = s.rstrip().split('\t')

        if fields[0] not in mutation_types:
            raise Error("Error constructing mutation object from line: " + s)

        self.type = fields[0]
        self.id = fields[1]
        self.parent_ids = fields[2]
        self.seq_id = fields[3]
        self.position = int(fields[4])
        self.strand = '.'
        self.optional = {}

        if self.type == 'SNP':
            self.optional['new_seq'] = fields[5]
        elif self.type == 'SUB':
            self.optional['size'] = int(fields[5])
            self.optional['new_seq'] = fields[6]
        elif self.type == 'DEL':
            self.optional['size']= int(fields[5])
        elif self.type == 'INS':
            self.optional['new_seq'] = fields[5]
        elif self.type == 'MOB':
            self.optional['repeat_name'] = fields[5]
            if fields[6] == '1':
                self.strand = '+'
            elif fields[6] == '-1':
                self.strand = '-'
            else:
                raise Error("Error getting strand " + line)
            self.optional['duplication_size'] = int(fields[7])
        elif self.type == 'AMP':
            self.optional['size'] = int(fields[5])
            self.optional['new_copy_number'] = int(fields[6])
        elif self.type == 'CON':
            self.optional['size'] = int(fields[5])
            self.optional['region'] = fields[6]
        elif self.type == 'INV':
            self.optional['size'] = int(fields[5])


    def to_gff(self):
        l = [self.seq_id, 'x', self.type, self.position]
        if 'size' in self.optional:
            l.append(self.position + self.optional['size'] - 1)
        else:
            l.append(self.position)

        l.extend(['.', self.strand, '.'])
        extra = [k + '=' + str(self.optional[k]) for k in sorted(self.optional)]
        l.append(';'.join(extra))
        return('\t'.join([str(x) for x in l]))


class GenomeDiff:
    def __init__(self, filename):
        f = utils.open_file_read(filename)

        self.version = None
        self.mutations = {}  # (seq name, pos) -> [list of mutations]

        for line in f:
            # first line should define that this is a genome diff file
            if self.version is None:
                if not line.startswith('#=GENOME_DIFF'):
                    raise Error("Error. first line of file '" + filename + "' should start with: #=GENOME_DIFF")

                self.version = line.rstrip().split()[-1]
                continue

            # for now, ignore the rest of the metadata
            if line.startswith('#'):
                continue

            fields = line.rstrip().split('\t')

            if fields[0] in mutation_types:
                mutation = Mutation(line)
                self.mutations[mutation.seq_id, mutation.position] = mutation

        utils.close(f)


    def write_gff(self, filename):
        # sort the output by reference name then position
        f = utils.open_file_write(filename)

        for k in sorted(self.mutations.keys()):
            print(self.mutations[k].to_gff(), file=f)

        utils.close(f)
