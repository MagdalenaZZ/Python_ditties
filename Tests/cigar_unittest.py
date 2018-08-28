#!/usr/bin/env python3.3

import sys
sys.path.insert(1, '..')
import cigar
import unittest
import fastn

class TestCigar(unittest.TestCase):
    def test_init(self):
        '''__init__ should do the expected with a cigar string'''
        with self.assertRaises(cigar.Error):
            cigar.Cigar("5M2")
        with self.assertRaises(cigar.Error):
            cigar.Cigar("H5M2X")

        test_cig = '5S2M3I1M1D4H'
        self.assertEqual(test_cig, str(cigar.Cigar(test_cig)))


    def test_ref_hit_length(self):
        '''Check that ref_hit_length() returns the right number'''
        self.assertEqual(cigar.Cigar("10M").ref_hit_length(), 10)
        self.assertEqual(cigar.Cigar("1S9M").ref_hit_length(), 9)
        self.assertEqual(cigar.Cigar("1S7M1I1M").ref_hit_length(), 8)
        self.assertEqual(cigar.Cigar("1S7M1D1M").ref_hit_length(), 9)

    def test_read_hit_length(self):
        '''Check that read_hit_length() returns the right number'''
        self.assertEqual(cigar.Cigar("10M").read_hit_length(), 10)
        self.assertEqual(cigar.Cigar("1S9M").read_hit_length(), 9)
        self.assertEqual(cigar.Cigar("1S7M1I1M").read_hit_length(), 9)
        self.assertEqual(cigar.Cigar("1S7M1D1M").read_hit_length(), 8)

    def test_reverse(self):
        '''Test that reverse works as expected'''
        c = cigar.Cigar('1S10M1I5M2S')
        c.reverse()
        self.assertEqual(str(c), '2S5M1I10M1S')

    def test_get_differences_from_ref(self):
        '''check test_get_differences_from_ref finds the correct differences'''
        ref = fastn.Fasta('ID', 'ACGTACGTACGT')
        c = cigar.Cigar("12M")

        pairs_to_check = [(cigar.Cigar("12M"), 'ACGTACGTACGT'),
                          (cigar.Cigar("12M"), 'AGGTACGTACGT'),
                          (cigar.Cigar("1S12M"), 'AAGGTACGTACGT'),
                          (cigar.Cigar("1S12M1S"), 'AAGGTACGTACGTA'),
                          (cigar.Cigar("1M1I10M"), 'AiCGTACGTACGT'),
                          (cigar.Cigar("3M1I3M1D3M"), 'AGGiTACTACGT'),
                          (cigar.Cigar("2S3M1I3M1D3M5S"), 'ssAGGiTACTACGTsssss')]
        correct_answers = [[],
                           [(1, 'S', 'C/G', 1)],
                           [(1, 'S', 'C/G', 1)],
                           [(1, 'S', 'C/G', 1)],
                           [(1, 'I', 'i', 1)],
                           [(1, 'S', 'C/G', 1), (3, 'I', 'i', 1), (6, 'D', 'G', 1)],
                           [(1, 'S', 'C/G', 1), (3, 'I', 'i', 1), (6, 'D', 'G', 1)]]

        for i in range(len(pairs_to_check)):
            self.assertListEqual(pairs_to_check[i][0].get_differences_from_ref(pairs_to_check[i][1], ref), correct_answers[i])


    def test_soft_clipped_bases(self):
        '''Check that sort_clipped_bases() counts the right number of bases'''
        self.assertEqual(cigar.Cigar("10M").soft_clipped_bases(), 0)
        self.assertEqual(cigar.Cigar("1S10M").soft_clipped_bases(), 1)
        self.assertEqual(cigar.Cigar("10M2S").soft_clipped_bases(), 2)
        self.assertEqual(cigar.Cigar("2S10M3S").soft_clipped_bases(), 5)

if __name__ == '__main__':
    unittest.main()
