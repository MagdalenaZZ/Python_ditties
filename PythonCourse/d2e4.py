
import datetime
import d2e3

def fastFib(n, memo={}):


    if n == 0 or n == 1:
        return 1
    
    try:
        return memo[n]
    except KeyError:
        memo[n] = fastFib(n-1,memo) + fastFib(n-2,memo)
        return memo[n]


def testFib():
    for i in [0,1,5,20]:
        print 'Fibonacci of', i, '=', fastFib(i)

def testFibnodic():
    for i in [0,1,5,20]:
        print 'Fibonacci of', i, '=', d2e3.fib(i)


a = datetime.datetime.now().time()
testFib()
b = datetime.datetime.now().time()
testFibnodic()
c = datetime.datetime.now().time()


a =float(str(a).split(":")[2])
b =float(str(b).split(":")[2])
c =float(str(c).split(":")[2])

print 'First took', b-a
print 'Second took', c-b

