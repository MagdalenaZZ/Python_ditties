


def selSort(L):
    """Assumes that L is a list of elements that can be compared using > """

    suffixStart=0

    while suffixStart != len(L):
    # Look at each element in suffix

        for i in xrange(suffixStart, len(L)):
            if L[i] < L[suffixStart]:
                # swap position of elements
                L[i], L[suffixStart] = L[suffixStart] , L[i]
        suffixStart +=1






