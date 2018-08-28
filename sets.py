

# find overlaps

a = set(["Jake", "John", "Eric"])
b = set(["John", "Jill"])

print(a.intersection(b))
print(b.intersection(a))



#To find out which members attended only one of the events, use the "symmetric_difference" method:

a = set(["Jake", "John", "Eric"])
b = set(["John", "Jill"])

print(a.symmetric_difference(b))
print(b.symmetric_difference(a))


# To find out which members attended only one event and not the other, use the "difference" method:



# To receive a list of all participants, use the "union" method

a = set(["Jake", "John", "Eric"])
b = set(["John", "Jill"])

print(a.union(b))




