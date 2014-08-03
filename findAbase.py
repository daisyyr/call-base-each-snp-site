#!/usr/lib/python

'''
find A base and basecount in the result file of callbase_HTSeq.py
return a tuple(base, basecount)
'''

def findreadA(line):
    baselist = ['A', 'T', 'C', 'G']
    if line.startswith('ref-name'):
        print line
    else:
        Adict = {}
        basecount = line.split()[3:7]
        for count, base in zip(basecount, baselist):
            Adict[count] = base
            neworder = sorted(Adict, key=lambda x:int(x), reverse=True)
        if line.split()[2] == 'A':
            if neworder[0] != basecount[0]:
                return Adict[neworder[0]], int(neworder[0])
            else:
                return Adict[neworder[1]], int(neworder[1])
        if line.split()[2] == 'T':
            if neworder[0] != basecount[1]:
                return Adict[neworder[0]], int(neworder[0])
            else:
                return Adict[neworder[1]], int(neworder[1])
        if line.split()[2] == 'C':
            if neworder[0] != basecount[2]:
                return Adict[neworder[0]], int(neworder[0])
            else:
                return Adict[neworder[1]], int(neworder[1])
        if line.split()[2] == 'G':
            if neworder[0] != basecount[3]:
                return Adict[neworder[0]], int(neworder[0])
            else:
                return Adict[neworder[1]], int(neworder[0])

