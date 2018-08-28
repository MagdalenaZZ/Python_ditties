#!/usr/bin/env python3.3

import sys
import os
import filecmp
sys.path.insert(1, '..')
import genome_intervals
import unittest


class TestIntervals(unittest.TestCase):
    def test_init(self):
        '''Throw error if try to construct genome_interval from a non-int, or end<start'''
        with self.assertRaises(genome_intervals.Error):
            genome_intervals.Interval('a', 1)
        with self.assertRaises(genome_intervals.Error):
            genome_intervals.Interval(1, 'a')
        with self.assertRaises(genome_intervals.Error):
            genome_intervals.Interval('a', 'a')
        with self.assertRaises(genome_intervals.Error):
            genome_intervals.Interval(3, 2)

    def test_comparisons(self):
        '''< and <= should work as expected'''
        self.assertTrue(genome_intervals.Interval(1,2) < genome_intervals.Interval(2,2))
        self.assertTrue(genome_intervals.Interval(1,2) <= genome_intervals.Interval(2,2))
        self.assertFalse(genome_intervals.Interval(2,2) <= genome_intervals.Interval(1,2))
        self.assertFalse(genome_intervals.Interval(2,2) < genome_intervals.Interval(1,2))
        self.assertFalse(genome_intervals.Interval(2,2) < genome_intervals.Interval(2,2))

    def test_len(self):
        self.assertEqual(len(genome_intervals.Interval(1,2)), 2)
        self.assertEqual(len(genome_intervals.Interval(1,1)), 1)
        self.assertEqual(len(genome_intervals.Interval(10,20)), 11)

    def test_intersects(self):
        '''Intersection of two intervals should do the right thing'''
        a = genome_intervals.Interval(5, 10)
        no_intersect = [genome_intervals.Interval(3, 4),
                        genome_intervals.Interval(11,20)]
        intersect = [genome_intervals.Interval(3,5),
                     genome_intervals.Interval(3,6),
                     genome_intervals.Interval(9,12),
                     genome_intervals.Interval(10,12),
                     genome_intervals.Interval(6,7),
                     genome_intervals.Interval(1,20)]

        for i in no_intersect:
            self.assertFalse(a.intersects(i), 'shouldn\'t intersect: ' + str(a) + ', ' + str(i))

        for i in intersect:
            self.assertTrue(a.intersects(i), 'should intersect: ' + str(a) + ', ' + str(i))

    def test_contains(self):
        '''Check that contains() works as expected'''
        a = genome_intervals.Interval(5, 10)
        not_contained = [genome_intervals.Interval(1,2),
                         genome_intervals.Interval(4,5),
                         genome_intervals.Interval(4,10),
                         genome_intervals.Interval(4,11),
                         genome_intervals.Interval(5,11),
                         genome_intervals.Interval(1,2),
                         genome_intervals.Interval(9,11),
                         genome_intervals.Interval(10,11),
                         genome_intervals.Interval(11,20)]


        contained = [genome_intervals.Interval(5,5),
                     genome_intervals.Interval(5,10),
                     genome_intervals.Interval(6,7),
                     genome_intervals.Interval(6,10),
                     genome_intervals.Interval(10,10)]

        for i in not_contained:
            self.assertFalse(a.contains(i), 'shouldn\'t contain: ' + str(a) + ', ' + str(i))

        for i in contained:
            self.assertTrue(a.contains(i), 'should contain: ' + str(a) + ', ' + str(i))

    def test_union(self):
        '''Union should either return None or the correct union'''
        a = genome_intervals.Interval(5, 10)
        b = genome_intervals.Interval(8, 15)
        c = genome_intervals.Interval(12, 20)
        d = genome_intervals.Interval(21,22)
        self.assertEqual(a.union(c), None)
        self.assertEqual(c.union(a), None)
        self.assertEqual(a.union(b), genome_intervals.Interval(5,15))
        self.assertEqual(b.union(a), genome_intervals.Interval(5,15))
        self.assertEqual(c.union(d), genome_intervals.Interval(12,22))
        self.assertEqual(d.union(c), genome_intervals.Interval(12,22))

    def test_union_flll_gap(self):
        '''union_fill_gap() should ignore intersections and return the maximum range of coords'''
        a = genome_intervals.Interval(5, 10)
        b = genome_intervals.Interval(8, 15)
        c = genome_intervals.Interval(12, 20)
        d = genome_intervals.Interval(21,22)
        self.assertEqual(a.union_fill_gap(c), genome_intervals.Interval(5,20))
        self.assertEqual(c.union_fill_gap(a), genome_intervals.Interval(5,20))
        self.assertEqual(a.union_fill_gap(b), genome_intervals.Interval(5,15))
        self.assertEqual(b.union_fill_gap(a), genome_intervals.Interval(5,15))
        self.assertEqual(c.union_fill_gap(d), genome_intervals.Interval(12,22))
        self.assertEqual(d.union_fill_gap(c), genome_intervals.Interval(12,22))


    def test_intersection(self):
        '''Intersection should either return None or the correct intersection'''
        a = genome_intervals.Interval(5, 10)
        b = genome_intervals.Interval(8, 15)
        c = genome_intervals.Interval(12, 20)
        self.assertEqual(a.intersection(c), None)
        self.assertEqual(a.intersection(b), genome_intervals.Interval(8,10))

class Test_intersection(unittest.TestCase):
    def test_intersection(self):
        '''intersection() should correctly intersect two lists of intervals'''
        a = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(10,20),
             genome_intervals.Interval(51,52),
             genome_intervals.Interval(54,55),
             genome_intervals.Interval(57,58)]

        b = [genome_intervals.Interval(5,6),
             genome_intervals.Interval(9,11),
             genome_intervals.Interval(13,14),
             genome_intervals.Interval(17,18),
             genome_intervals.Interval(20,25),
             genome_intervals.Interval(50,60)]

        i = [genome_intervals.Interval(10,11),
             genome_intervals.Interval(13,14),
             genome_intervals.Interval(17,18),
             genome_intervals.Interval(20,20),
             genome_intervals.Interval(51,52),
             genome_intervals.Interval(54,55),
             genome_intervals.Interval(57,58)]

        self.assertSequenceEqual(genome_intervals.intersection(a,b), i)
        self.assertSequenceEqual(genome_intervals.intersection(b,a), i)


class Test_merge_overlapping_in_list(unittest.TestCase):
    def test_merge_overlapping_in_list(self):
        '''merge_overlapping_in_list() merges correctly'''
        a = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(51,60),
             genome_intervals.Interval(10,20),
             genome_intervals.Interval(20,30),
             genome_intervals.Interval(20,30),
             genome_intervals.Interval(29,50),
             genome_intervals.Interval(65,70)]

        b = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(10,60),
             genome_intervals.Interval(65,70)]

        genome_intervals.merge_overlapping_in_list(a)
        self.assertSequenceEqual(a, b)

class Test_remove_contained_in_list(unittest.TestCase):
    def test_remove_contained_in_list(self):
        '''test_remove_contained_in_list removes the right elements of list'''
        a = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(4,4),
             genome_intervals.Interval(4,5),
             genome_intervals.Interval(5,6),
             genome_intervals.Interval(7,9),
             genome_intervals.Interval(8,10),
             genome_intervals.Interval(9,11),
             genome_intervals.Interval(20,25),
             genome_intervals.Interval(20,24),
             genome_intervals.Interval(20,26),
             genome_intervals.Interval(30,38),
             genome_intervals.Interval(30,37),
             genome_intervals.Interval(30,36),
             genome_intervals.Interval(30,35),
             genome_intervals.Interval(30,35),
             genome_intervals.Interval(32,33),
             genome_intervals.Interval(38,50),
             genome_intervals.Interval(65,70),
             genome_intervals.Interval(67,70)]

        b = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(4,5),
             genome_intervals.Interval(5,6),
             genome_intervals.Interval(7,9),
             genome_intervals.Interval(8,10),
             genome_intervals.Interval(9,11),
             genome_intervals.Interval(20,26),
             genome_intervals.Interval(30,38),
             genome_intervals.Interval(38,50),
             genome_intervals.Interval(65,70)]

        genome_intervals.remove_contained_in_list(a)
        self.assertSequenceEqual(a, b)

class Test_length_sum_from_list(unittest.TestCase):
    def test_length_sum_from_list(self):
        '''Test that total length of intervals is summed correctly'''
        a = [genome_intervals.Interval(1,2),
             genome_intervals.Interval(4,5),
             genome_intervals.Interval(10,19)]

        self.assertEqual(14, genome_intervals.length_sum_from_list(a))


if __name__ == '__main__':
    unittest.main()
