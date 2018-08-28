import re

class Error (Exception): pass

class Operation:
    def __init__(self, letter, number):
        self.operator = letter
        self.number = int(number)

    def __str__(self):
        return str(self.number) + self.operator


class Cigar:
    def __init__(self, string_in=None):
        self.operations = []

        if string_in is not None:
            try:
                numbers = re.split('[A-Z]', string_in)[:-1]
                letters = re.split('[0-9]*', string_in)[1:]
                assert(len(numbers) == len(letters))
            except:
                raise Error('Error making Cigar object from this string: "' +  string_in + '"')

            for i in range(len(letters)):
                self.operations.append(Operation(letters[i], numbers[i]))

    def __str__(self):
        return ''.join([str(x) for x in self.operations])

    def reverse(self):
        self.operations.reverse()

    def ref_hit_length(self):
        return sum([o.number for o in self.operations if o.operator in ['M','N','D']])

    def read_hit_length(self):
        return sum([o.number for o in self.operations if o.operator in ['M','I']])

    def soft_clipped_bases(self):
        return sum([o.number for o in self.operations if o.operator == 'S'])

    def get_differences_from_ref(self, read, ref, ref_start=0):
        diffs = []
        ref_pos = ref_start
        read_pos = 0
        for o in self.operations:
            if o.operator == 'M':
                for i in range(o.number):
                    if read[read_pos] != ref[ref_pos]:
                        diffs.append((ref_pos, 'S', ref[ref_pos] + '/' + read[read_pos], 1))
                    ref_pos += 1
                    read_pos += 1
            elif o.operator == 'S':
                read_pos += o.number
            elif o.operator == 'I':
                diffs.append((ref_pos, 'I', read[read_pos:read_pos + o.number], o.number))
                read_pos += o.number
            elif o.operator == 'D':
                diffs.append((ref_pos, 'D', ref[ref_pos:ref_pos + o.number], o.number))
                ref_pos += o.number

        return diffs
