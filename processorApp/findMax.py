def findMax(x):
    maxi = x[0]
    index = 0
    for i in range(len(x)):
        if x[i] >= maxi:
            maxi = x[i]
            index = i

    return index, maxi
