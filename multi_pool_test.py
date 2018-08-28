#!/usr/bin/env python3

import multiprocessing
import fastaq
import time


class Test:
    def __init__(self, outprefix, threads):
        self.outprefix = outprefix
        self.threads = threads


    def func(self, p):
        print('start', p)
        f = fastaq.utils.open_file_write(self.outprefix + '.' + p)
        for i in range(20):
            print(p, i, file=f)
            time.sleep(0.2)
        fastaq.utils.close(f)
        print('finished', p)


    def multi_test(self):
        p = ['one', 'two', 'three']
        pool = multiprocessing.Pool(self.threads)
        pool.map(self.func, p)
        pool.close()
        pool.join()


t = Test('out', 2)

t.multi_test()





print('All done!')

