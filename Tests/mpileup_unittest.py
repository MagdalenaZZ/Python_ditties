#!/usr/bin/env python3.3

import sys
import os
sys.path.insert(1, '..')
import filecmp
import unittest
import mpileup
import utils

class TestMpileupLine(unittest.TestCase):
    def test_init_SamRecord(self):
        '''__init__ for MpileupLine should do as expected'''
        testline = 'MAL1:78:189868-191088\t1\tN\t1\t^]A\tH'

        self.assertEqual(testline, str(mpileup.MpileupLine(testline)))

        bad_testlines = ['MAL1:78:189868-191088\t1\tN\t1\t^]A',
                         'MAL1:78:189868-191088\tX\tN\t1\t^]A\tH',
                         'MAL1:78:189868-191088\t1\tN\tX\t^]A\tH']

        for l in bad_testlines:
            with self.assertRaises(mpileup.Error):
                x = mpileup.MpileupLine(l)


class TestFileReader(unittest.TestCase):
    def test_file_reader_mpileup(self):
        '''file_reader should iterate through a pileup file correctly'''
        tmp_out = 'tmp.mpileup'
        fout = utils.open_file_write(tmp_out)
        mpileup_reader = mpileup.file_reader('mpileup_unittest.mpileup')
        for mp in mpileup_reader:
            print(mp, file=fout)
        utils.close(fout)
        self.assertTrue(filecmp.cmp('mpileup_unittest.mpileup', tmp_out))
        os.unlink(tmp_out)



class TestGetMeanCoveragePerSeq(unittest.TestCase):
    def test_get_mean_covereage_per_seq(self):
        '''Check get_mean_coverage_per_seq does as expected, including missing seqs'''
        found_covs = mpileup.get_mean_coverage_per_seq('mpileup_unittest.mpileup', 'mpileup_unittest.mpileup.fai')
        expected_covs = {'MAL1': 1.0 * (1+2+3+3+3+5+6) / 7,
                         'MAL2': 1.0 * (6+6+6+7+7) / 10,
                         'MAL3': 1.0 * (7+8+8+8+9+9+9) / 30,
                         'MAL4': 0.0}

        self.assertSequenceEqual(found_covs.keys(), expected_covs.keys())
        for chr in expected_covs:
            self.assertAlmostEqual(found_covs[chr], expected_covs[chr])


if __name__ == '__main__':
    unittest.main()
