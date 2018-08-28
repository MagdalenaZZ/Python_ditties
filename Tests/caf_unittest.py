#!/usr/bin/env python3.3

import sys
import os
import copy
sys.path.insert(1, '..')
import fastn
import unittest
import utils
import caf


class TestCaf(unittest.TestCase):
    def test_get_next_from_file(self):
        '''get_next_from_file() should read caf records from file correctly'''

        f_in = utils.open_file_read('caf_unittest.caf')

        c = caf.Caf()
        c.get_next_from_file(f_in)

        e = caf.Caf()
        e.id = 'pknbac5b2Aa01.p1k'
        seq = ''.join(['NGGAGAGACTCGGACTAGTTCTACACCCTCACACCTTTGTCCTAAACCTTGAATCTAAGT'
                       'CCTAACACCCTGACACCTTTGTCCTAAGCCCGGAATCTAACTTCTAGCACCCCTACGACC',
                       'CTTATTCCTAAACCCAGAATCTGACTATTGACACCCCTACAACCCTAATTCCAACACCCT',
                       'TACAACCTTCATTCCAACACCGCAACAACCTTCATTCCAGCACCCCAACAACCTTCATTC',
                       'CAACACCCCAAACAACATCATTCCAACACCCCAAACAACATCATTCCAACACCCCAAACA',
                       'ACATCATTCCAACACCCCAAACAACATCATTCCAACACGGCAACAACATCATTCGAACAC',
                       'CCCTACAACATCATTCCAGCACCCCAACAACCTCCCTGCGAAACCCCGAATCCGAATTTT',
                       'GACACCCCTACAACCTTATTCTGACACCCCCAACAAACTTTCTCTAACACCCCAACAACG',
                       'TGACTACTAATACACCTAAAACCTTACTCCTAAACCCGGAATCCGACTTCTAATACCGCA',
                       'ACAACCTTCATTCCTAAACCCGGAATCTGAACCCTGAACCATTAAAACATAAAACGTGGA',
                       'AAATGAACCCCTGAACCATGAAAACCGTGAAAACCTATAACTTGGACCATGAACCTCTCA',
                       'ACCCCGAAATATGAGAACTTTGGAAACCCTAAATTTTGGGAAAACTCCTTTTTTTTTTTT',
                       'TTATTGTACATCCTGTGCGATGGTATACATTTTGGCGAATGCAAAAGAATTAGCATATAT',
                       'ATATGTGTAGGTCTTTGTGATGGTCAGGGGGGAGATCGACTAGGGTGTAGGTCTTTGTGA',
                       'TGGTCAAGGGAGATGGGCCAAAGGGAAGTCGGACAAGGTGAGATGGGCCAAGGAGATGGG',
                       'CCTAGGGTGGATGGGACAAGGGTGGATGGTCAGAGGTGGATGGTCAAGGGTGGATGGTCA',
                       'AGGATGAATGGGCAAGGGAGATGGGCAAAGTAGATGGGCAAGGGTGGATGGACAAGGTGG',
                       'ATGGCCAAAGTGGATGGCAAGGAGGATGGCCCAGGTAATAGGCAAGGAAATGGCCAGGTG',
                       'GATGGACCAGGTGGTGCCCTAATGGAGGCAGGGTGAAGTCCAGGAGGAGGCCCAGGAAAA',
                       'GGCCCAGAGAAACCCAAGGAAAGGCCCAGGGGGTGGGACAGGGGAAGCGCCAAGGGATGC',
                       'CAAGGTGGGGGCCAGAAAATAGCCCAGAAAAGGCCAAAATAAGCCAAGAAAAGCCCCAGA',
                       'AAACCCAAGAAA'])

        quals = [4, 4, 4, 4, 6, 6, 8, 6, 6, 6, 6, 10, 12, 11, 13, 13, 20, 19, 9, 10, 9, 9, 9, 19, 19, 34, 34, 39, 35, 35, 35, 37, 35, 34, 26, 26, 16, 17, 11, 21, 21, 32, 35, 37, 37, 32, 45, 23, 17, 17, 18, 27, 29, 32, 35, 32, 32, 32, 32, 39, 35, 35, 35, 35, 35, 37, 42, 31, 31, 14, 13, 13, 25, 25, 35, 40, 33, 29, 23, 23, 15, 25, 24, 35, 35, 35, 35, 23, 36, 18, 18, 23, 28, 33, 29, 29, 32, 32, 32, 32, 35, 35, 32, 35, 35, 32, 35, 35, 44, 44, 37, 35, 28, 26, 24, 19, 23, 30, 33, 40, 32, 32, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 50, 50, 50, 37, 30, 30, 27, 27, 21, 21, 21, 29, 26, 29, 29, 23, 23, 28, 37, 37, 50, 50, 40, 35, 35, 32, 32, 32, 35, 44, 37, 35, 35, 35, 35, 35, 35, 32, 32, 35, 35, 35, 35, 44, 42, 42, 41, 41, 41, 41, 41, 42, 41, 41, 41, 41, 41, 41, 44, 44, 42, 42, 42, 42, 42, 35, 37, 35, 35, 33, 37, 37, 44, 44, 44, 41, 42, 50, 42, 42, 42, 44, 44, 50, 50, 44, 44, 44, 44, 44, 44, 50, 50, 44, 44, 44, 44, 44, 41, 42, 44, 42, 42, 42, 44, 44, 42, 42, 41, 41, 41, 42, 44, 50, 50, 50, 44, 44, 44, 44, 44, 37, 37, 37, 37, 39, 41, 41, 44, 44, 44, 44, 47, 47, 44, 44, 44, 43, 43, 42, 42, 37, 37, 37, 41, 41, 42, 44, 44, 44, 44, 44, 42, 42, 42, 41, 41, 41, 44, 44, 44, 46, 42, 41, 37, 37, 37, 37, 37, 41, 42, 35, 35, 35, 35, 35, 35, 35, 42, 41, 42, 44, 50, 42, 42, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 37, 37, 41, 44, 44, 47, 37, 37, 33, 33, 33, 27, 27, 37, 37, 47, 47, 47, 47, 47, 50, 44, 44, 42, 50, 35, 35, 35, 42, 42, 44, 50, 50, 50, 42, 42, 42, 42, 35, 35, 37, 42, 50, 44, 44, 44, 44, 44, 44, 47, 47, 47, 47, 44, 50, 44, 44, 44, 44, 47, 47, 44, 47, 50, 50, 50, 48, 37, 17, 17, 13, 22, 22, 35, 36, 42, 42, 35, 35, 35, 37, 37, 42, 50, 35, 35, 35, 35, 37, 37, 35, 35, 33, 33, 33, 33, 42, 42, 42, 41, 41, 41, 41, 41, 41, 50, 37, 44, 44, 44, 42, 37, 37, 21, 21, 21, 33, 33, 42, 50, 50, 44, 44, 44, 44, 44, 44, 44, 42, 37, 44, 44, 44, 44, 42, 42, 42, 42, 42, 44, 44, 44, 50, 50, 44, 44, 44, 37, 37, 35, 33, 33, 21, 21, 33, 33, 33, 41, 42, 42, 41, 44, 44, 44, 44, 42, 42, 42, 42, 44, 41, 44, 37, 42, 37, 41, 41, 42, 42, 50, 50, 44, 44, 44, 44, 42, 42, 27, 33, 27, 33, 33, 37, 37, 50, 35, 35, 35, 37, 37, 44, 44, 50, 44, 44, 44, 37, 37, 35, 31, 31, 37, 37, 44, 44, 44, 44, 50, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 37, 37, 28, 28, 23, 28, 26, 33, 33, 33, 29, 29, 29, 33, 35, 46, 33, 23, 23, 26, 33, 33, 50, 44, 37, 37, 30, 37, 37, 42, 50, 30, 30, 30, 37, 37, 23, 28, 15, 15, 11, 15, 27, 37, 33, 37, 26, 26, 28, 37, 42, 48, 48, 37, 23, 23, 23, 31, 31, 33, 23, 23, 24, 31, 31, 31, 33, 24, 25, 21, 21, 21, 28, 31, 33, 42, 42, 42, 42, 44, 44, 44, 44, 30, 23, 16, 10, 10, 16, 24, 33, 24, 24, 24, 30, 33, 36, 42, 42, 44, 44, 42, 39, 39, 33, 46, 27, 28, 28, 33, 33, 37, 37, 37, 22, 22, 17, 19, 19, 33, 31, 33, 27, 27, 18, 18, 24, 29, 32, 33, 35, 33, 40, 40, 37, 34, 27, 27, 14, 14, 13, 13, 18, 12, 20, 25, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 47, 56, 56, 56, 47, 42, 42, 42, 42, 27, 23, 15, 11, 11, 27, 33, 42, 42, 33, 24, 10, 10, 10, 13, 15, 18, 14, 14, 14, 14, 14, 25, 30, 33, 30, 21, 21, 27, 27, 22, 15, 13, 22, 22, 19, 19, 15, 15, 11, 10, 17, 27, 27, 21, 15, 18, 13, 13, 16, 22, 24, 37, 31, 40, 40, 37, 47, 40, 37, 27, 27, 24, 24, 17, 20, 13, 10, 10, 11, 11, 14, 12, 19, 10, 10, 12, 14, 11, 10, 10, 10, 10, 15, 11, 15, 15, 25, 12, 12, 8, 8, 10, 17, 10, 21, 21, 8, 8, 8, 10, 10, 19, 25, 21, 19, 10, 10, 8, 9, 10, 12, 14, 17, 24, 22, 16, 16, 10, 9, 8, 14, 12, 12, 9, 9, 9, 9, 9, 9, 13, 19, 15, 18, 22, 22, 15, 15, 15, 15, 9, 10, 9, 8, 8, 9, 10, 14, 10, 10, 19, 15, 12, 9, 15, 4, 4, 4, 8, 8, 10, 12, 9, 8, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 16, 10, 10, 10, 8, 7, 7, 7, 7, 7, 13, 20, 19, 15, 15, 10, 10, 8, 8, 10, 10, 10, 15, 10, 8, 8, 9, 8, 9, 10, 11, 10, 8, 8, 8, 8, 8, 4, 8, 4, 7, 7, 9, 13, 16, 11, 10, 12, 11, 13, 8, 8, 8, 8, 8, 9, 10, 9, 9, 9, 8, 8, 8, 12, 8, 9, 9, 11, 10, 10, 7, 7, 9, 7, 8, 9, 11, 10, 9, 10, 9, 10, 7, 7, 7, 9, 8, 8, 10, 8, 8, 4, 7, 4, 4, 4, 4, 4, 8, 7, 7, 8, 9, 9, 7, 7, 9, 9, 9, 9, 8, 7, 7, 7, 7, 10, 10, 7, 8, 8, 9, 10, 10, 10, 14, 13, 9, 8, 7, 7, 7, 6, 6, 7, 7, 6, 6, 6, 6, 6, 8, 15, 10, 8, 8, 8, 8, 6, 7, 6, 6, 6, 6, 7, 7, 7, 8, 7, 4, 4, 4, 6, 6, 6, 6, 7, 13, 7, 7, 8, 8, 8, 7, 7, 7, 7, 7, 7, 8, 9, 7, 7, 7, 6, 6, 6, 6, 6, 7, 9, 7, 7, 7, 8, 10, 8, 8, 8, 8, 9, 6, 6, 6, 6, 6, 6, 6, 7, 7, 10, 9, 9, 9, 8, 8, 8, 8, 8, 7, 7, 7, 6, 6, 6, 7, 6, 6, 7, 9, 7, 7, 11, 6, 6, 7, 6, 6, 8, 7, 7, 8, 8, 10, 8, 8, 8, 6, 6, 7, 6, 6, 6, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 7, 7, 7, 7, 7, 12, 9, 14, 10, 10, 10, 10, 8, 9, 8, 8, 8, 7, 7, 7, 7, 7, 13, 7, 7, 6, 6, 6, 6, 6, 6, 8, 7, 7, 7, 6, 6, 6, 6, 6, 8, 8, 8, 9, 7, 7, 7, 8, 9, 7, 6, 6, 6, 6, 6, 6, 8, 8, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 8, 9, 8, 8, 8, 8, 8, 8, 8, 12, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 6, 9, 6, 6, 6, 7, 7, 7, 7, 7, 8, 10, 10, 14, 9, 12, 7, 7, 7, 4, 4]
        e.seq = fastn.Fasta(e.id, seq)
        e.seq = e.seq.to_Fastq(quals)
        e.insert_min = 2000
        e.insert_max = 4000
        e.ligation = '96781'
        e.clone = 'pknbac5b2'
        e.clip_start = 23
        e.clip_end = 789

        self.assertEqual(c, e)

        c.get_next_from_file(f_in)

        e = caf.Caf()
        e.id = 'pknbac5b2Aa02.p1k'

        seq = ''.join(['AAAGACATACGACCTTTTTTTTTTTCGATAACAAAGGGTATCCTTTCACCAGAAAAAAAA',
                       'AAAGAACATTCTTCTTTTTTCTTGAAGAACATACATTCTTTTTTTTATTTTATTTTTTTT',
                       'TTTCGACCCCTCAGTGTTGTGGTAGCATGATGTGTTGGACTTGAATGGTATATGTATTGA',
                       'TTGTTTCGTTCGTTATGTAATTTCCGGTTTTTCCCCGTGGCATCCGGATAGTGTATAGTA',
                       'TCCGGTCCCTGTGTTCAAAAAGTTTTTCCTTTTCCCCTTAAAGCAACTGAAGTTAAACCC',
                       'TGAACCTTACTACTGAACCCGGAATTTGACTTCTAAAACCCTGAAGAATGATTCCTATAA',
                       'CCCTAAAAAATCCAACCTAAAACATCCAAACTGAACCATAGAACCTTCCTCCTAAACCCG',
                       'GAATCTATGTTCTAACACCCTGACATCTTTGTCCTAAACCCTGAATCTAAGTTCTAACAT',
                       'CCTGACAACTCTCCCTCCTAAACCCGGAATCTAAATTCGTACACCCTGACACCTCCCCCC',
                       'TAAACCCGGAATCCGCATTCTAACACCCTGACAATTTCCTCCTGAAAAGCGGAATCTGAC',
                       'TTCTAACACCCTGACACCTTTGTCCTGAACCCGGAATCTAAGTTCTTACACCCGGACACC',
                       'TCCCTCCTAAATCCGGAATCTAAGTTCTAACACCCTCACACCTTTGTCCTAAACCTTGAA',
                       'TCTAAGTCCTAACACCCTGACACCTTTGTCCTAAGCCCGGAATCTAACTTCTAGCACCCC',
                       'TACGACCCTTATTCCTAAACCCAGAATCTGACTATTGACACCCCTACAACCCTAATTCCA',
                       'ACACCCTTACAACCTTCATTCCAACACCGCAACAACCTTCATTCCAGCACCCCTACAACT',
                       'TCATTCCTACACCCCAAACAACATCATCCCTACACCCCAAACAACATCATTCCTACACCC',
                       'CAAACACATCATCCAACACCCCATAACACATCATTCCAACACGGCAACAACATCATTCGA',
                       'AACACCCCTACAAATCATTGCAGCACCCCCACTACCTCCCTGCGTATACCCGTATTCGAA',
                       'ATTTTGACACCCCTACTACCTTTATCTGACACCCCCAAAAAACTCCTCTTAAACCCAACA',
                       'AGGGGACTATAATACCCCTAAAACTTTATCTTAACCGGAATCCGAATTCTATACCGAAAA',
                       'AACTTCTTTCCTAACCGGGATCTGTACCCCGAACTTTTAAAATTAAAGGGGAAATGAACC',
                       'CCTGACCAGATAACGGGAAACCTTTATTGTGACAGGAACTCCTACCGCAATATGAAAATT',
                       'GGACCCCAAATTTGGGAAACCCCTTTT'])


        quals = [9, 9, 6, 4, 4, 4, 4, 7, 6, 6, 8, 6, 6, 6, 7, 7, 14, 8, 8, 8, 10, 17, 21, 12, 9, 10, 10, 9, 11, 8, 9, 11, 11, 21, 12, 15, 15, 21, 24, 33, 32, 35, 29, 29, 22, 22, 15, 29, 25, 26, 18, 18, 18, 31, 31, 47, 56, 56, 56, 42, 36, 44, 28, 28, 28, 39, 33, 35, 30, 36, 33, 35, 35, 36, 35, 37, 42, 35, 35, 31, 29, 26, 26, 20, 33, 15, 22, 22, 29, 29, 32, 35, 35, 36, 35, 35, 42, 42, 37, 37, 42, 47, 47, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 44, 47, 47, 47, 47, 47, 35, 30, 30, 23, 24, 30, 45, 37, 37, 37, 35, 23, 23, 11, 11, 13, 23, 31, 21, 21, 19, 20, 23, 29, 23, 20, 16, 16, 30, 29, 29, 28, 28, 28, 24, 24, 17, 29, 29, 33, 33, 35, 31, 37, 18, 15, 12, 16, 16, 23, 27, 24, 32, 29, 32, 32, 24, 26, 29, 37, 29, 30, 35, 35, 33, 35, 35, 31, 33, 31, 31, 35, 31, 31, 31, 31, 27, 33, 33, 42, 35, 37, 37, 21, 21, 21, 21, 37, 37, 50, 50, 50, 50, 50, 33, 33, 18, 16, 15, 25, 19, 20, 33, 33, 33, 35, 35, 33, 33, 33, 18, 18, 18, 33, 24, 33, 33, 33, 27, 33, 33, 33, 33, 33, 22, 33, 33, 33, 24, 24, 21, 24, 24, 31, 31, 11, 11, 11, 31, 33, 44, 44, 37, 42, 42, 47, 44, 44, 44, 44, 44, 44, 44, 47, 50, 50, 42, 42, 42, 41, 42, 42, 47, 47, 37, 37, 27, 33, 33, 33, 33, 35, 35, 42, 41, 37, 37, 44, 50, 50, 33, 33, 27, 33, 37, 42, 42, 42, 41, 41, 33, 33, 27, 27, 33, 33, 37, 50, 35, 35, 35, 35, 35, 35, 35, 42, 35, 37, 35, 37, 35, 41, 37, 42, 42, 42, 42, 50, 50, 50, 42, 35, 33, 33, 21, 21, 16, 23, 19, 27, 27, 33, 35, 41, 50, 37, 35, 35, 42, 50, 50, 50, 44, 44, 44, 50, 42, 42, 37, 37, 35, 35, 35, 44, 44, 50, 50, 41, 37, 37, 37, 37, 35, 35, 35, 37, 37, 37, 44, 37, 37, 33, 33, 22, 33, 37, 35, 33, 33, 21, 21, 21, 33, 33, 41, 41, 44, 44, 44, 44, 44, 50, 50, 44, 44, 37, 50, 33, 33, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 50, 50, 50, 50, 50, 44, 47, 44, 44, 48, 33, 21, 24, 21, 33, 33, 35, 37, 50, 50, 37, 35, 35, 50, 50, 56, 50, 50, 50, 50, 48, 33, 27, 33, 27, 33, 33, 44, 50, 50, 42, 37, 35, 42, 42, 50, 50, 50, 44, 44, 44, 33, 33, 18, 19, 18, 33, 33, 42, 35, 35, 44, 44, 44, 50, 44, 44, 44, 50, 50, 44, 44, 37, 37, 33, 33, 35, 35, 35, 35, 35, 37, 50, 37, 27, 27, 24, 37, 33, 35, 35, 37, 35, 37, 37, 46, 33, 24, 24, 21, 33, 33, 39, 42, 42, 44, 50, 50, 56, 50, 50, 37, 35, 35, 33, 37, 33, 33, 35, 35, 35, 35, 33, 33, 33, 33, 27, 27, 27, 37, 37, 44, 37, 41, 41, 41, 50, 46, 33, 24, 24, 16, 31, 19, 27, 31, 37, 37, 44, 44, 44, 37, 50, 23, 23, 22, 29, 31, 33, 23, 23, 23, 23, 23, 28, 25, 33, 26, 26, 22, 28, 37, 42, 44, 42, 42, 44, 44, 44, 44, 46, 33, 16, 19, 14, 27, 31, 42, 50, 50, 50, 44, 44, 44, 50, 50, 26, 26, 21, 28, 31, 29, 29, 26, 26, 26, 30, 30, 39, 27, 37, 26, 30, 30, 42, 42, 42, 36, 33, 29, 33, 33, 33, 20, 21, 23, 17, 23, 31, 36, 42, 43, 56, 56, 47, 47, 42, 42, 33, 33, 29, 29, 23, 31, 25, 26, 26, 26, 30, 30, 36, 27, 33, 28, 31, 33, 35, 44, 33, 33, 28, 33, 35, 44, 48, 48, 48, 42, 47, 42, 42, 42, 48, 44, 44, 37, 34, 34, 44, 48, 42, 37, 34, 42, 48, 33, 33, 34, 30, 30, 33, 33, 40, 30, 37, 28, 28, 26, 27, 27, 25, 19, 16, 25, 29, 40, 31, 27, 15, 18, 13, 25, 27, 40, 40, 33, 40, 33, 33, 33, 40, 37, 23, 12, 12, 17, 11, 10, 15, 15, 13, 13, 13, 18, 27, 23, 28, 28, 28, 28, 37, 28, 32, 26, 23, 26, 26, 19, 29, 25, 24, 25, 24, 15, 15, 15, 12, 17, 24, 24, 21, 21, 21, 25, 22, 29, 25, 22, 21, 24, 25, 17, 17, 14, 14, 12, 14, 19, 24, 18, 18, 14, 21, 11, 15, 10, 15, 18, 22, 27, 25, 25, 29, 29, 29, 25, 26, 25, 21, 22, 25, 22, 22, 18, 15, 15, 15, 25, 19, 25, 25, 16, 24, 24, 20, 20, 22, 20, 15, 10, 10, 10, 12, 13, 20, 20, 12, 14, 14, 12, 12, 12, 15, 15, 15, 18, 18, 11, 10, 11, 11, 10, 10, 14, 15, 18, 18, 19, 17, 12, 11, 10, 10, 20, 15, 19, 24, 24, 24, 23, 15, 13, 7, 6, 6, 6, 6, 6, 12, 13, 12, 9, 8, 10, 10, 9, 6, 6, 6, 6, 10, 10, 13, 15, 15, 15, 15, 17, 9, 9, 9, 9, 9, 11, 11, 9, 7, 7, 7, 6, 4, 4, 6, 9, 9, 8, 8, 8, 10, 9, 8, 7, 7, 7, 7, 7, 9, 13, 10, 10, 10, 15, 12, 9, 9, 9, 15, 19, 15, 15, 11, 7, 7, 7, 7, 7, 7, 8, 8, 19, 10, 10, 10, 12, 12, 19, 11, 15, 18, 11, 14, 9, 9, 6, 6, 6, 6, 6, 6, 6, 8, 11, 20, 13, 17, 14, 14, 9, 9, 10, 17, 11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 11, 10, 12, 11, 10, 12, 9, 12, 8, 8, 8, 9, 12, 12, 8, 11, 7, 8, 8, 8, 8, 11, 9, 8, 6, 4, 4, 4, 6, 6, 7, 10, 10, 12, 9, 7, 7, 6, 6, 6, 6, 8, 6, 9, 10, 13, 8, 11, 8, 7, 7, 8, 7, 7, 7, 10, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 14, 10, 8, 12, 8, 8, 6, 7, 9, 8, 7, 6, 8, 7, 4, 4, 7, 7, 6, 7, 6, 6, 6, 6, 8, 11, 8, 8, 8, 8, 12, 10, 12, 11, 11, 11, 10, 12, 10, 7, 7, 9, 4, 4, 8, 6, 6, 6, 6, 6, 6, 7, 10, 7, 7, 7, 7, 7, 9, 9, 9, 7, 7, 7, 6, 6, 6, 7, 7, 7, 10, 11, 9, 7, 6, 6, 8, 6, 6, 8, 8, 8, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 12, 10, 15, 15, 16, 7, 7, 6, 6, 6, 7, 6, 6, 6, 6, 6, 8, 7, 7, 8, 7, 8, 7, 7, 9, 8, 7, 7, 8, 8, 9, 7, 6, 7, 6, 9, 6, 7, 11, 7, 7, 11, 8, 8, 7, 10, 8, 9, 8, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 8, 8, 7, 7, 6, 9, 7, 6, 6, 6, 6, 6, 6, 6, 8, 7, 7, 7, 7, 7, 7, 7, 10, 12, 19, 13, 13, 10, 9, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 6, 6, 6, 6, 6, 4, 6, 4, 4, 4, 4, 4, 4, 6, 6, 6, 7, 7, 7, 7, 8, 7, 6, 6, 6, 6, 6, 6, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 11, 12, 9]

        e.seq = fastn.Fasta(e.id, seq)
        e.seq = e.seq.to_Fastq(quals)
        e.insert_min = 2000
        e.insert_max = 4000
        e.ligation = '96781'
        e.clone = 'pknbac5b2'
        e.clip_start = 33
        e.clip_end = 848
        self.assertEqual(c, e)

        utils.close(f_in)


if __name__ == '__main__':
    unittest.main()
