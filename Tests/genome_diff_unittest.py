#!/usr/bin/env python3.3

import sys
import os
import filecmp
sys.path.insert(1, '..')
import genome_diff
import unittest


class TestGenomeDiff(unittest.TestCase):
    def test_write_gff(self):
        '''Check conversion to gff file works'''
        gff_out = 'tmp.gff'
        gd = genome_diff.GenomeDiff('genome_diff_unittest.gd')
        gd.write_gff(gff_out)
        self.assertTrue(filecmp.cmp('genome_diff_unittest.gff', gff_out))
        os.unlink(gff_out)


if __name__ == '__main__':
    unittest.main()
