import sys

# Find the cube root of a perfect cube

x = int( raw_input('Enter an integer: '))
ans = 0

# step up ans by 1 
while ans*ans*ans < abs(x):
    ans = ans+1
    print (ans , "is now")

    
    # if it is the same as cube, great

if ans*ans*ans == abs(x):
    print ("Cube root of ", x, "is", ans)






sys.exit()


if ans*ans*ans != abs(x):
    print (x, "is not a perfect cube")

else:
    if x < 0:
        ans = -ans
    print ("Cube root of ", x, "is", ans)
print ("Done")



