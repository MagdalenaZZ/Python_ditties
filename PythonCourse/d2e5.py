import random
from math import pi

def inCircle():
    # return the next random floating point number between 0 and 1
    x = random.random()
    y = random.random()

    #(x2 + y2)^0.5 less than 1
    if (x*x + y*y) **0.5  <= 1.0:
            return True
    else:
            return False



# I think numneedles is set to 1000, unless given a value, in which case it uses that instead
def estPi(numNeedles = 1000):
    # this isn't even needed? maybe used by random.random()
    #random.seed(123)
    print numNeedles
    count = 0 
    # from 1 to numNeedles, can be range too
    for Needles in xrange(1, numNeedles+1,1):
        # evaluates as true or false
        if inCircle():
                count +=1
    # the value of count is now a number between 0 and 1000, most likely around 7500
    print count, Needles

    # The estimation is 4 times (7500 divided by Needles, which is now = numNeedles 
    est = 4*(count/float(Needles))


    # print the estimation, the diff between est and pi, and the numNeedles
    print 'Est. = %7.5f, Diff.= %7.5f, Needles = %d' % (est, (est-pi),numNeedles)
    return est

#print(inCircle())

#est = estPi()
#print (est)

