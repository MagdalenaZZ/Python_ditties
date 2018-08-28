#!/usr/bin/env python3.3

import sys
import os
sys.path.insert(1, '..')
import blast
import filecmp
import unittest
import utils

class TestBlast(unittest.TestCase):
    def test_init_BlastHit(self):
        '''Check that init followed by str gives the same as the original input'''

        testlines = ['Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.01\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7\t100\t200']

        for t in testlines:
            self.assertEqual(t, str(blast.BlastHit(t)))

        bad_testlines = ['Contig1\t1.1.B1.36.cap3_contig\tX\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51.1\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51\tX\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51\t1\tA\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7\t100']

        for l in bad_testlines:
            with self.assertRaises(blast.Error):
                self.assertEqual(l, str(blast.BlastHit(l)))

    def test_add_sequence_lengths(self):
        '''Check that add_sequence_lengths() does the job and fails when it should'''
        b = blast.BlastHit('qry\tref\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7')
        b_with_length = blast.BlastHit('qry\tref\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7\t999\t1000')
        d_ref = {}
        d_qry = {}
        with self.assertRaises(blast.Error):
            b.add_sequence_lengths(d_ref, d_qry)

        d_qry['qry'] = 999
        with self.assertRaises(blast.Error):
            b.add_sequence_lengths(d_ref, d_qry)

        d_ref['ref'] = 1000
        with self.assertRaises(blast.Error):
            b.add_sequence_lengths(d_qry, d_ref)

        b.add_sequence_lengths(d_ref, d_qry)
        self.assertEqual(b, b_with_length)


class TestFileReader(unittest.TestCase):
    def test_file_reader(self):
        '''file_reader should iterate through a blast file correctly'''
        tmp_out = 'blast_unittest.m8.tmp'

        for f in ['blast_unittest.m8', 'blast_unittest.m8.with_lengths']:
            blast_reader = blast.file_reader('blast_unittest.m8')
            fout = utils.open_file_write(tmp_out)
            for hit in blast_reader:
                print(hit, file=fout)
            utils.close(fout)
            self.assertTrue(filecmp.cmp('blast_unittest.m8', tmp_out))
            os.unlink(tmp_out)


class TestAddSequenceLengths(unittest.TestCase):
    def test_add_sequence_lengths(self):
        '''check add_sequence_lengths() works as expected'''
        tmp_out = 'blast_unittest.m8.with_lengths.tmp'
        blast.add_sequence_lengths('blast_unittest.m8', 'blast_unittest.m8.fai', 'blast_unittest.m8.fai', tmp_out)
        self.assertTrue(filecmp.cmp('blast_unittest.m8.with_lengths', tmp_out))
        os.unlink(tmp_out)


if __name__ == '__main__':
    unittest.main()
