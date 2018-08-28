#!/usr/bin/python3

import sys, getopt



######################################################################

# reads in arguments

def main(argv):
    number = ''
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"h:i:n:o:",["num=","ifile=","ofile="])
        #print("Opts:", opts)
    except getopt.GetoptError:
        print ('test.py -n <number> -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -n <number> -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-n", "--num"):
            number = arg
 
    print ('Input file is ', inputfile)
    print ('Output file is ', outputfile)
    print ('Number is ', number)
    return number
    
######################################################################

# make the input string a number if possible 

def make_number(x):
    #print(x)

    try:
        x = int(x)
        print ("Is an integer")
    except ValueError as e:
        print("Not an integer then. ")

        # how about float
        try:
            x=float(x)
        except ValueError as e:
            print("Not a float either. Try inputting a number instead of  ",x)
            sys.exit()

    return x
    print("Let's get on")

######################################################################


# make float integer

def float_2_int(x):

    # string split and multiply with the number of zeros

    #words = text.split(",")
    whole, small = str(x).split(".")
    lens = int(len(str(small)))
    zero = str("0")
    multiplier = int(str(1)+lens*zero) 
    print(lens, multiplier)
    newval  = int(multiplier * x)
    return newval,multiplier

######################################################################


def perfect_cube(x):
    ans = 0

# step up ans by 1
    while ans < (x/3)+0.5:
        #print ans 
        if ans**3==x:
            #print ("Perfect cube",ans,x)
            return ans
        else:
            ans +=1

    return "No perfect cube exists"


######################################################################

# allows you to both define a module and run it

if __name__ == "__main__":
    number = main(sys.argv[1:])


######################################################################
######################################################################



# check the number is a number, and change float to val
number = make_number(number)

#print("Num:",number, type(number))


# if it is a float fix it, print results
if   type(number)==float:
    # make int 
    new, mult = float_2_int(number)
    #print("Num:",number, type(number), new, mult)
    #number=int(new)
    res = perfect_cube(new)

    if type(res)==str:
        print ('\n',res, 'of', number )
    else:
        res = int(perfect_cube(new))/mult*100
        print ('\nThe prefect cube of {0:.3f} is {1}'.format(res,number))

else:
    #print("Is integer ",number)
    res = perfect_cube(number)

    if type(res)==str:
        print ('\n',res, 'of', number )
    else:
        print ('\nThe prefect cube of {0} is {1}'.format(perfect_cube(number), number))








