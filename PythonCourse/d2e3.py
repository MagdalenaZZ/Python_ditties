


def fib(n):
        """Assumes n is an integer >=0
        Returns Fibonacci of n"""

        if n==0 or n==1:
                return 1
        else:
            return fib(n-1) + fib(n-2)


def testFib():
    for i in range(16):
        print 'Fibonacci of', i, '=', fib(i)

# run testfib
#if __name__ == "__main__":
#    testFib()



