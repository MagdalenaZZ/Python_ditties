#!/usr/bin/env python3.3

import sys
import filecmp
import os
import copy
sys.path.insert(1, '..')
import fastn
import unittest
import utils
import genome_intervals


class TestFasta(unittest.TestCase):
    def setUp(self):
        self.fasta = fastn.Fasta('ID', 'ACGTA')


    def test_init(self):
        '''__init__ should get the ID and sequence correctly'''
        self.assertEqual(self.fasta.id, 'ID')
        self.assertEqual(self.fasta.seq, 'ACGTA')


    def test_get_next_from_file(self):
        '''get_next_from_file() should read seqs from OK, including weirdness in file'''
        f_in = utils.open_file_read('fastn_unittest.fa')
        fa = fastn.Fasta()
        counter = 1

        while fa.get_next_from_file(f_in):
            self.assertEqual(fa, fastn.Fasta(str(counter), 'ACGTA'))
            counter += 1

        utils.close(f_in)

    def test_getitem(self):
        '''getitem() should return the right subsequence'''
        seq = 'AACGTGTCA'
        fa = fastn.Fasta('x', seq)
        self.assertEqual(seq[1], fa[1])
        self.assertEqual(seq[0:2], fa[0:2])
        self.assertEqual(seq[1:], fa[1:])

    def test_len(self):
        '''len() should return the length of the sequence'''
        self.assertEqual(5, len(self.fasta))

    def test_print_line_length(self):
        '''__str__ should be formatted correctly with the right number of chars per line of sequence'''
        line_lengths = [0, 3]
        correct_files = ['fastn_unittest_one-per-line.fa', 'fastn_unittest_3-per-line.fa']

        for i in range(len(line_lengths)):
            seq_reader = fastn.file_reader('fastn_unittest_one-per-line.fa')
            fastn.Fasta.line_length = line_lengths[i]
            tmp_out = 'tmp.line_length_test.fa'
            f = utils.open_file_write(tmp_out)
            for s in seq_reader:
                print(s, file=f)
            utils.close(f)
            self.assertTrue(filecmp.cmp(correct_files[i], tmp_out))
            os.unlink(tmp_out)

        fastn.Fasta.line_length = 60

    def test_strip_illumina_suffix(self):
        '''Check that /1 and /2 removed correctly from IDs'''
        seqs = [fastn.Fasta('name/1', 'A'),
                fastn.Fasta('name/2', 'A'),
                fastn.Fasta('name', 'A'),
                fastn.Fasta('name/1/2', 'A'),
                fastn.Fasta('name/2/1', 'A'),
                fastn.Fasta('name/3', 'A')]

        correct_names = ['name', 'name', 'name', 'name/1', 'name/2', 'name/3']

        for seq in seqs:
            seq.strip_illumina_suffix()

        for i in range(len(seqs)):
            self.assertEqual(seqs[i].id, correct_names[i])

    def test_revcomp(self):
        '''revcomp() should correctly reverse complement a sequence'''
        fa = fastn.Fasta('ID', 'ACGTNacgtn')
        fa.revcomp()
        self.assertEqual(fa, fastn.Fasta('ID', 'nacgtNACGT'))

    def test_gaps(self):
        '''gaps() should find the gaps in a sequence correctly'''
        test_seqs = [fastn.Fasta('ID', 'ACGT'),
                     fastn.Fasta('ID', 'NACGT'),
                     fastn.Fasta('ID', 'NACGTN'),
                     fastn.Fasta('ID', 'ANNCGT'),
                     fastn.Fasta('ID', 'NANNCGTNN')]

        correct_gaps = [[],
                        [genome_intervals.Interval(0, 0)],
                        [genome_intervals.Interval(0, 0), genome_intervals.Interval(5, 5)],
                        [genome_intervals.Interval(1, 2)],
                        [genome_intervals.Interval(0, 0), genome_intervals.Interval(2, 3), genome_intervals.Interval(7, 8)]]

        for i in range(len(test_seqs)):
            gaps = test_seqs[i].gaps()
            self.assertListEqual(correct_gaps[i], gaps)

    def test_contig_coords(self):
        '''contig_coords() should get the coords of all contigs in a sequence correctly'''
        test_seqs = [fastn.Fasta('ID', 'ACGT'),
                     fastn.Fasta('ID', 'NACGT'),
                     fastn.Fasta('ID', 'NNACGT'),
                     fastn.Fasta('ID', 'ACGTN'),
                     fastn.Fasta('ID', 'ACGTNN'),
                     fastn.Fasta('ID', 'NANNCGT'),
                     fastn.Fasta('ID', 'ANNCGTNNAAAAA')]

        correct_coords = [[genome_intervals.Interval(0,3)],
                         [genome_intervals.Interval(1, 4)],
                         [genome_intervals.Interval(2, 5)],
                         [genome_intervals.Interval(0, 3)],
                         [genome_intervals.Interval(0, 3)],
                         [genome_intervals.Interval(1, 1), genome_intervals.Interval(4,6)],
                         [genome_intervals.Interval(0, 0), genome_intervals.Interval(3, 5), genome_intervals.Interval(8, 12)]]

        for i in range(len(test_seqs)):
            gaps = test_seqs[i].contig_coords()
            self.assertListEqual(correct_coords[i], gaps)


    def test_trim_Ns(self):
        '''trim_Ns() should do the right trimming of a sequence'''
        fa = fastn.Fasta('ID', 'ANNANA')
        test_seqs = [fastn.Fasta('ID', 'ANNANA'),
                     fastn.Fasta('ID', 'NANNANA'),
                     fastn.Fasta('ID', 'NANNANAN'),
                     fastn.Fasta('ID', 'ANNANAN'),
                     fastn.Fasta('ID', 'NNNNNNANNANAN'),
                     fastn.Fasta('ID', 'NNANNANANn')]

        for s in test_seqs:
            s.trim_Ns()
            self.assertEqual(fa, s)

    def test_replace_bases(self):
        '''Check that bases get replaced correctly'''
        fa = fastn.Fasta('X', 'AUCGTUUACT')
        fa.replace_bases('U', 'T')
        self.assertEqual(fa, fastn.Fasta('X', 'ATCGTTTACT'))

    def test_search_string(self):
        '''Check that search_string() finds all the hits'''
        fa = fastn.Fasta('X', 'AAA')
        hits = fa.search('G')
        self.assertTrue(len(hits) == 0)
        hits = fa.search('AAA')
        self.assertListEqual(hits, [(0, '+')])
        hits = fa.search('AA')
        self.assertListEqual(hits, [(0, '+'), (1, '+')])
        hits = fa.search('TTT')
        self.assertListEqual(hits, [(0, '-')])

    def test_to_Fastq(self):
        '''Check to_Fastq converts OK, including out of range quality scores'''
        fa = fastn.Fasta('X', 'AAAAA')
        quals = [-1, 0, 40, 93, 94]
        self.assertEqual(fastn.Fastq('X', 'AAAAA', '!!I~~'), fa.to_Fastq(quals))


    def test_translate(self):
        '''Test nucleatide -> amino acid conversion works on Fasta'''
        fa = fastn.Fasta('ID', 'GCAGCCGCGGCTAGAAGGCGACGCCGGCGTAACAATGACGATTGCTGTGAAGAGCAACAGGGAGGCGGGGGTCACCATATAATCATTTTATTGCTACTCCTGCTTAAAAAGATGTTCTTTCCACCCCCGCCTAGCAGTTCATCCTCGTCTACAACCACGACTTGGTACTATGTAGTCGTGGTTTAATAGTGA')

        self.assertEqual(fastn.Fasta('ID', 'AAAARRRRRRNNDDCCEEQQGGGGHHIIILLLLLLKKMFFPPPPSSSSSSTTTTWYYVVVV***'), fa.translate())


    def test_split_capillary_id(self):
        '''Tests that we get information from a sanger capillary read name OK'''
        ids = ['abcde.p1k', 'abcde.x.p1k', 'abcde.p1ka', 'abcde.q1k', 'abcde.w2k']
        expected = [{'prefix': 'abcde', 'dir': 'fwd', 'suffix': 'p1k'},
                    {'prefix': 'abcde.x', 'dir': 'fwd', 'suffix': 'p1k'},
                    {'prefix': 'abcde', 'dir': 'fwd', 'suffix': 'p1ka'},
                    {'prefix': 'abcde', 'dir': 'rev', 'suffix': 'q1k'},
                    {'prefix': 'abcde', 'dir': 'unk', 'suffix': 'w2k'}]

        for i in range(len(ids)):
            fa = fastn.Fasta(ids[i], 'A')
            self.assertEqual(fa.split_capillary_id(), expected[i])


class TestFastq(unittest.TestCase):
    def setUp(self):
        self.fastq = fastn.Fastq('ID', 'ACGTA', 'IIIII')

    def test_init(self):
        '''__init__ should get the ID, sequence and quality correctly'''
        self.assertEqual(self.fastq.id, 'ID')
        self.assertEqual(self.fastq.seq, 'ACGTA')
        self.assertEqual(self.fastq.qual, 'IIIII')

    def test_init_length_mismatch(self):
        '''__init__ should raise an error when length of seq and quality not the same'''
        with self.assertRaises(fastn.Error):
            fastn.Fastq('X', 'A', 'II')

    def test_get_next_from_file(self):
        '''get_next_from_file() should read seqs from OK, and raise error at badly formatted file'''
        bad_files = ['fastn_unittest_fail_no_AT.fq',
                     'fastn_unittest_fail_no_seq.fq',
                     'fastn_unittest_fail_no_plus.fq',
                     'fastn_unittest_fail_no_qual.fq']

        for fname in bad_files:
            f_in = utils.open_file_read(fname)
            fq = fastn.Fastq()
            with self.assertRaises(fastn.Error):
                while fq.get_next_from_file(f_in):
                    pass

            utils.close(f_in)

        fname = 'fastn_unittest_good_file.fq'
        try:
            f_in = open(fname)
        except IOError:
            print("Error opening '" + fname + "'", file=sys.stderr)
            sys.exit(1)

        fq = fastn.Fastq()
        while fq.get_next_from_file(f_in):
            self.assertEqual(fq, fastn.Fastq('ID', 'ACGTA', 'IIIII'))
        utils.close(f_in)

    def test_revcomp(self):
        '''revcomp() should correctly reverse complement a sequence'''
        fq = fastn.Fastq('ID', 'ACGTNacgtn', '1234567890')
        fq.revcomp()
        self.assertEqual(fq, fastn.Fastq('ID', 'nacgtNACGT', '0987654321'))

    def test_trim_Ns(self):
        '''trim_Ns() should do the right trimming of a fastq sequence'''
        fq = fastn.Fastq('ID', 'ANNANA', '111111')
        test_seqs = [fastn.Fastq('ID', 'ANNANA', '111111'),
                     fastn.Fastq('ID', 'NANNANA', '1111111'),
                     fastn.Fastq('ID', 'NANNANAN', '11111111'),
                     fastn.Fastq('ID', 'ANNANAN', '1111111'),
                     fastn.Fastq('ID', 'NNNNNNANNANAN', '1111111111111'),
                     fastn.Fastq('ID', 'NNANNANANn', '1111111111')]

        for s in test_seqs:
            s.trim_Ns()
            self.assertEqual(fq, s)

    def test_trim(self):
        '''trim() should trim the right number of bases off start and end'''
        fq = fastn.Fastq('ID', '1234567890', '1234567890')
        fq.trim(0, 0)
        self.assertEqual(fq, fastn.Fastq('ID', '1234567890', '1234567890'))

        fq = fastn.Fastq('ID', '1234567890', '1234567890')
        fq.trim(1, 0)
        self.assertEqual(fq, fastn.Fastq('ID', '234567890', '234567890'))

        fq = fastn.Fastq('ID', '1234567890', '1234567890')
        fq.trim(0, 1)
        self.assertEqual(fq, fastn.Fastq('ID', '123456789', '123456789'))

        fq = fastn.Fastq('ID', '1234567890', '1234567890')
        fq.trim(2, 2)
        self.assertEqual(fq, fastn.Fastq('ID', '345678', '345678'))

    def test_to_Fasta_and_qual(self):
        '''Check to_Fasta_and_qual converts quality scores correctly'''
        fq = fastn.Fastq('ID', 'ACGT', '>ADI')
        (fa, qual) = fq.to_Fasta_and_qual()
        self.assertEqual(fa, fastn.Fasta('ID', 'ACGT'))
        self.assertListEqual(qual, [29, 32, 35, 40])


    def test_translate(self):
        '''Test nucleatide -> amino acid conversion works on Fasta'''
        fq = fastn.Fastq('ID', 'GCAGCCGCGGCTAGAAGGCGACGCCGGCGTAACAATGACGATTGCTGTGAAGAGCAACAGGGAGGCGGGGGTCACCATATAATCATTTTATTGCTACTCCTGCTTAAAAAGATGTTCTTTCCACCCCCGCCTAGCAGTTCATCCTCGTCTACAACCACGACTTGGTACTATGTAGTCGTGGTTTAATAGTGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII')

        self.assertEqual(fastn.Fastq('ID', 'AAAARRRRRRNNDDCCEEQQGGGGHHIIILLLLLLKKMFFPPPPSSSSSSTTTTWYYVVVV***', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'), fq.translate())

class TestFileReader(unittest.TestCase):
    def test_file_reader_fasta(self):
        '''file_reader should iterate through a fasta file correctly'''
        reader = fastn.file_reader('fastn_unittest.fa')
        counter = 1
        for seq in reader:
            self.assertEqual(seq, fastn.Fasta(str(counter), 'ACGTA'))
            counter += 1

    def test_file_reader_fastq(self):
        '''file_reader should iterate through a fastq file correctly'''
        reader = fastn.file_reader('fastn_unittest_good_file.fq')
        for seq in reader:
            self.assertEqual(seq, fastn.Fastq('ID', 'ACGTA', 'IIIII'))

class TestDeinterleave(unittest.TestCase):
    def test_deinterleave(self):
        '''deinterleave should deal with an interleaved file correctly'''
        tmp_1 = 'tmp.deinterleaved_1.fa'
        tmp_2 = 'tmp.deinterleaved_2.fa'
        fastn.deinterleave('fastn_unittest_interleaved.fa', tmp_1, tmp_2)
        self.assertTrue(filecmp.cmp('fastn_unittest_deinterleaved_1.fa', tmp_1))
        self.assertTrue(filecmp.cmp('fastn_unittest_deinterleaved_2.fa', tmp_2))

        with self.assertRaises(fastn.Error):
            fastn.deinterleave('fastn_unittest_interleaved_bad.fa', tmp_1, tmp_2)
        os.unlink(tmp_1)
        os.unlink(tmp_2)

class TestInterleave(unittest.TestCase):
    def test_interleave(self):
        '''Check that interleave works as expected'''
        tmp = 'tmp.interleaved.fa'
        fastn.interleave('fastn_unittest_deinterleaved_1.fa', 'fastn_unittest_deinterleaved_2.fa', tmp)
        self.assertTrue(filecmp.cmp('fastn_unittest_interleaved.fa', tmp))

        with self.assertRaises(fastn.Error):
            fastn.interleave('fastn_unittest_deinterleaved_bad_1.fa', 'fastn_unittest_deinterleaved_bad_2.fa', tmp)
        os.unlink(tmp)

class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        '''reverse_complement should correctly reverse complement each seq in a file'''
        tmp = 'tmp.revcomp.fa'
        fastn.reverse_complement('fastn_unittest.fa', tmp)
        self.assertTrue(filecmp.cmp('fastn_unittest_revcomp.fa', tmp))
        os.unlink(tmp)

class TestTrim(unittest.TestCase):
    def test_trim(self):
        '''trim should correctly trim each seq in a file'''
        tmp = 'tmp.trim.fq'
        fastn.trim('fastn_unittest_untrimmed.fq', tmp, 2, 1)
        self.assertTrue(filecmp.cmp('fastn_unittest_trimmed.fq', tmp))
        os.unlink(tmp)

class TestFastqToMiraXml(unittest.TestCase):
    def test_fastq_to_mira_xml(self):
        '''check that fastq_to_mira_xml makes the correct xml file from a fastq file'''
        tmp = 'tmp.mira.xml'
        fastn.fastq_to_mira_xml('fastn_unittest_good_file.fq', tmp)
        self.assertTrue(filecmp.cmp('fastn_unittest_good_file_mira.xml', tmp))
        os.unlink(tmp)

class TestFileToDict(unittest.TestCase):
    def test_file_to_dict(self):
        '''check file_to_dict9 fills dictionary correctly'''
        d_test = {}
        d = {}
        fastn.file_to_dict('fastn_unittest.fa', d_test)
        for i in range(1,5):
            d[str(i)] = fastn.Fasta(str(i),'ACGTA')

        self.assertSequenceEqual(d_test.keys(),d.keys())
        for i in range(1,5):
            key = str(i)
            self.assertEqual(d_test[key].id, d[key].id)
            self.assertEqual(d_test[key].seq, d[key].seq)


class TestLengthsFromFai(unittest.TestCase):
    def test_lengths_from_fai(self):
        '''Check lengths_from_fai gets the length of each seq OK'''
        d = {}
        lengths = {str(x):x for x in range(1,5)}
        fastn.lengths_from_fai('fastn_unittest_fai_test.fa.fai', d)
        self.assertSequenceEqual(d.keys(), lengths.keys())
        for i in d:
            self.assertEqual(int(i), d[i])


class TestSplitByBasecount(unittest.TestCase):
    def test_split_by_base_count(self):
        '''Check that fasta/q files get split by base count correctly'''
        infile = 'fastn_unittest_split_test.fa'
        outprefix = 'fastn_unittest_split_test.fa.test'
        length2files = {3: ['1','2','3'], 4: ['1', '2', '3'], 6: ['1', '2']}
        for l in length2files:
            fastn.split_by_base_count(infile, outprefix, l)
            for x in range(len(length2files[l])):
                file_index = str(length2files[l][x])
                fname = outprefix + '.' + file_index
                self.assertTrue(filecmp.cmp(fname, infile + '.' + str(l) + '.' + file_index))
                os.unlink(fname)

        # check that limiting the numebr of files works
        fastn.split_by_base_count(infile, outprefix, 6, 2)
        for i in range(1,4):
            test_file = outprefix + '.' + str(i)
            self.assertTrue(filecmp.cmp(test_file, 'fastn_unittest_split_test.fa.6.limit2.' + str(i)))
            os.unlink(test_file)



class TestCountSequences(unittest.TestCase):
    def test_count_sequences(self):
        '''Check that count_sequences does as expected'''
        self.assertEqual(2, fastn.count_sequences('fastn_unittest_good_file.fq'))
        self.assertEqual(4, fastn.count_sequences('fastn_unittest.fa'))
        self.assertEqual(0, fastn.count_sequences('empty_file'))

class TestFastaToFastq(unittest.TestCase):
    def test_fasta_to_fastq(self):
        '''Check fasta_to_fastq converts files as expected'''
        fastn.fasta_to_fastq('fastn_unittest.fa', 'fastn_unittest.fa.qual', 'tmp.fq')
        self.assertTrue(filecmp.cmp('fastn_unittest.fasta_to_fastq.fq', 'tmp.fq'))
        os.unlink('tmp.fq')

class TestFastnToQuasrPrimers(unittest.TestCase):
    def test_fastn_to_quasr_primers(self):
        '''Check that fasta file gets converted to QUASR sequence file'''
        fastn.fastn_to_quasr_primers('fastn_unittest_fastn_to_quasr_primers.fa', 'tmp.primers')
        self.assertTrue(filecmp.cmp('fastn_unittest_fastn_to_quasr_primers.expected', 'tmp.primers'))
        os.unlink('tmp.primers')

class TestFastnReplaceBases(unittest.TestCase):
    def test_fastn_replace_bases(self):
        '''Check that fasta file gets all bases replaced OK'''
        tmpfile = 'tmp.replace_bases.fa'
        fastn.replace_bases('fastn_unittest_fastn_replace_bases.fa', tmpfile, 'T', 'X')
        self.assertTrue(filecmp.cmp('fastn_unittest_fastn_replace_bases.expected.fa', tmpfile))
        os.unlink(tmpfile)

if __name__ == '__main__':
    unittest.main()

