#!/usr/bin/env python3.3

import sys
import os
sys.path.insert(1, '..')
import filecmp
import unittest
import utils
import external_progs

class TestBowtie2(unittest.TestCase):
    def test_bowtie2_index(self):
        '''Check bowtie2 index'''
        ref = 'external_progs_bowtie2_ref.fa'

        self.assertFalse(external_progs.is_bowtie2_indexed(ref))
        external_progs.index_with_bowtie2(ref)
        self.assertTrue(external_progs.is_bowtie2_indexed(ref))
        external_progs.clean_bowtie2_index(ref)
        self.assertFalse(external_progs.is_bowtie2_indexed(ref))


if __name__ == '__main__':
    unittest.main()
