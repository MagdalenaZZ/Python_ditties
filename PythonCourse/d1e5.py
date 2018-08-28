# finding the cube root of the perfect cube

#import sys

# Find the cube root of a perfect cube

x =  raw_input('Enter an integer: ')


# function isdigit will deal with non-integers
if x.isdigit() :
    x = int(x)
else:
    print "Not an integer"
    exit()


ans = 0

# step up ans by 1
while ans*ans*ans < abs(x):
    ans = ans+1

if ans*ans*ans != abs(x):
    print x, "is not a perfect cube"

else:
    if x < 0:
        ans = -ans
    print "Cube root of ", x, "is", ans
print ("Done")




#sys.exit()



