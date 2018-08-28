#!/usr/bin/env python3.3

import sys
import os
sys.path.insert(1, '..')
import sam
import filecmp
import unittest
import utils
import fastn

class TestSam(unittest.TestCase):
    def test_init_SamRecord(self):
        '''__init__ for SamRecord should do as expected'''
        testline = 'HS4_6280:2:1104:12102:124607\t99\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1'
        testline_no_tags = 'HS4_6280:2:1104:12102:124607\t99\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF'

        self.assertEqual(testline, str(sam.SamRecord(testline)))
        self.assertEqual(testline_no_tags, str(sam.SamRecord(testline_no_tags)))

        bad_testlines = ['HS4_6280:2:1104:12102:124607\t99\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT',
            'HS4_6280:2:1104:12102:124607\tAA\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1']

        for l in bad_testlines:
            with self.assertRaises(sam.Error):
                self.assertEqual(l, str(sam.SamRecord(l)))

    def test_flag_operations(self):
        '''Check all the methods involving flags do the right thing'''
        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t16\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertTrue(s.is_mapped())
        self.assertFalse(s.is_forward_strand())
        self.assertFalse(s.is_paired())
        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t4\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertFalse(s.is_mapped())


        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t147\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertTrue(s.is_proper_pair())
        self.assertTrue(s.is_paired())

        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t81\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertFalse(s.is_proper_pair())


        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t87\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertTrue(s.is_first_of_pair())
        self.assertFalse(s.is_second_of_pair())

        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t147\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertFalse(s.is_first_of_pair())
        self.assertTrue(s.is_second_of_pair())
        self.assertTrue(s.is_mate_mapped())

        self.assertEqual(s.query_strand(), '-')

        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t0\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertEqual(s.query_strand(), '+')
        self.assertTrue(s.is_forward_strand())

        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t95\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertFalse(s.is_mate_mapped())

        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t95\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertFalse(s.is_duplicate())
        s = sam.SamRecord('HS4_6280:2:1104:12102:124607\t1089\tPyYM_01_v1\t1\t47\t2S73M\t=\t362\t438\tTGTTAAAAATATCATTTATATAATATAATTAAAATTATTTATTTTTAGATATTATAATATTATGAATAATAGTAT\tHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGFHHHHHHHEHHHHHHHFCHFHHHHGFHCHHFHFAFFEFECF@BF\tAS:i:73\tNM:i:1')
        self.assertTrue(s.is_duplicate())

    def test_to_fastn(self):
        '''Check conversion to fastq with to_fastq()'''
        sams = [sam.SamRecord('ID\t0\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\tIIIII'),
                sam.SamRecord('ID\t16\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\tIIIII'),
                sam.SamRecord('ID\t65\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\tIIIII'),
                sam.SamRecord('ID\t129\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\tIIIII'),
                sam.SamRecord('ID\t0\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\t*'),
                sam.SamRecord('ID\t16\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\t*')]
        seqs = [fastn.Fastq('ID', 'ACGTA', 'IIIII'),
                fastn.Fastq('ID', 'TACGT', 'IIIII'),
                fastn.Fastq('ID/1', 'ACGTA', 'IIIII'),
                fastn.Fastq('ID/2', 'ACGTA', 'IIIII'),
                fastn.Fasta('ID', 'ACGTA'),
                fastn.Fasta('ID', 'TACGT')]

        for i in range(len(sams)):
            self.assertEqual(seqs[i], sams[i].to_fastn())

    def test_ref_hit_end_position(self):
        '''Check that we get the end position of a hit correctly'''
        s_good = sam.SamRecord('ID\t0\tref\t1\t47\t2S10M1I10M1D10M4H\t=\t362\t438\tACGTA\tIIIII')
        s_bad = sam.SamRecord('ID\t4\tref\t1\t47\t2S73M\t=\t362\t438\tACGTA\tIIIII')
        self.assertTrue(s_good.ref_hit_end_position(), 30)

        with self.assertRaises(sam.Error):
            s_bad.ref_hit_end_position()


class TestFileReader(unittest.TestCase):
    def test_file_reader_sam(self):
        '''file_reader should iterate through a BAM file correctly'''
        tmp_sam_out = 'tmp.sam'
        fout = utils.open_file_write(tmp_sam_out)
        sam_reader = sam.file_reader('sam_unittest.bam')
        for sam_record in sam_reader:
            print(sam_record, file=fout)
        utils.close(fout)
        self.assertTrue(filecmp.cmp('sam_unittest.sam', tmp_sam_out))
        os.unlink(tmp_sam_out)

class TestGetSequenceLengths(unittest.TestCase):
    def test_get_sequence_lengths(self):
        '''Check that we can read the sequence lengths from a BAM file header OK'''
        lengths = {
            'PyYM_01_v1': 639781,
            'PyYM_02_v1': 819697,
            'PyYM_03_v1': 746845,
            'PyYM_04_v1': 899110,
            'PyYM_05_v1': 1045897,
            'PyYM_06_v1': 1045168,
            'PyYM_07_v1': 931358,
            'PyYM_08_v1': 1563612,
            'PyYM_09_v1': 1783079,
            'PyYM_10_v1': 1760775,
            'PyYM_11_v1': 1889938,
            'PyYM_12_v1': 1882022,
            'PyYM_13_v1': 2639166,
            'PyYM_14_v1': 2614191,
            'PyYM_API_v1': 29736,
            'PyYM_bin_v1': 1734795,
            'PyYM_MIT_v1': 6512}

        d = sam.get_sequence_lengths('sam_unittest.bam')
        self.assertSequenceEqual(d.keys(), lengths.keys())
        for key in d:
            self.assertEqual(d[key], lengths[key])


if __name__ == '__main__':
    unittest.main()
