# numpy is a numerical library that provides additional (and very useful) functionality. Basically,
# its a bunch of code that already exists and we're reusing it. This next statement lets is telling
# Python to import two functions, exponential and normal, from the random module of numpy.
from numpy.random import exponential, normal

# pull 10000 random values from an exponential distribution
exponential_values = exponential(size=10000) 

# pull 20000 random values from a normal distribution with a mean of 1 and std of 2
normal_values = normal(loc=1, scale=2, size=20000)

# Create the histograms. The histogram function returns information about the histogram
# which we're storing in the variable hist_info simply so it doesn't clutter the output.
# We are passing in a few different pieces of information into these functions:
# color    : the color of the histogram
# alpha    : the level of transparency, where 0.0 is completely transparent and 1.0 is opaque
# bins     : the number of bins to use in the histogram
# histtype : the type of the histogram (in this case, we are not showing the boundaries on the bars) 
hist_info = hist(exponential_values, color='red', alpha=0.5, bins=100, histtype='stepfilled')
hist_info = hist(normal_values, color='blue', alpha=0.5, bins=100, histtype='stepfilled')

# let's label our plot
title("A pretty picture")
ylabel("Number of values")
xlabel("The actual value!")
legend(["Exponential", "Normal"])


