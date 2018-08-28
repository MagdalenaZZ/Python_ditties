#!/usr/bin/env python3.3

import sys
import os
sys.path.insert(1, '..')
import nucmer
import filecmp
import unittest
import utils

class TestSam(unittest.TestCase):
    def test_init_NucmerHit(self):
        '''Check that init followed by str gives the same as the original input'''

        testline = '1162\t25768\t24536\t4\t24607\t24533\t99.32\t640851\t24536\t1\t-1\tref1\tNODE_25757_length_24482_cov_18.920391\t[CONTAINS]'
        testline_no_tags = '1162\t25768\t24536\t4\t24607\t24533\t99.32\t640851\t24536\t1\t-1\tref1\tNODE_25757_length_24482_cov_18.920391'

        self.assertEqual(testline, str(nucmer.NucmerHit(testline)))
        self.assertEqual(testline_no_tags, str(nucmer.NucmerHit(testline_no_tags)))

        bad_testlines = ['1162\t25768\t24536\t4\t24607\t24533\t99.32\t640851\t24536\t1\t-1\tref1']

        for l in bad_testlines:
            with self.assertRaises(nucmer.Error):
                self.assertEqual(l, str(nucmer.NucmerHit(l)))


class TestFileReader(unittest.TestCase):
    def test_file_reader(self):
        '''file_reader should iterate through a nucmer file correctly'''
        tmp_out = 'nucmer_unittest.coords.tmp'
        fout = utils.open_file_write(tmp_out)
        nucmer_reader = nucmer.file_reader('nucmer_unittest.coords')
        for hit in nucmer_reader:
            print(hit, file=fout)
        utils.close(fout)
        self.assertTrue(filecmp.cmp('nucmer_unittest.coords.out', tmp_out))
        os.unlink(tmp_out)

if __name__ == '__main__':
    unittest.main()
