

def mySaveMat(M,fName):
    f = open(fName,'w+')

    l, w = M.shape

    f.write('%d\t%d\n'%(l,w))
    for i in range(l):
        for j in range(w):
            f.write('%15.22g\n'%(M[i,j]))

    f.close()

def myReadMat(fName):
    with open(fName,'r') as infile:
        
        l , w = [int(p) for p in infile.readline().split()]
        M = numpy.zeros(l,w)
        for i in range(l):
            for j in range(w):
                M[i,j] = infile.readline().astype('double')

    return M


