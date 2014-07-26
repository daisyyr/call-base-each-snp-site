#!/usr/lib/python
import os
from optparse import OptionParser

msg_usage = 'usage: %prog [-P] path [-O] output'
descr ='''Using split_vcf.py will get several results. This script will
concatenate those results by SNP sites coordinate.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-P', '--path', dest = 'path',
                     help = 'Those splited results in which directory?')
optparser.add_option('-O', '--output', dest = 'resultname',
                     help = 'Input the result file name.')
options, args = optparser.parse_args()

def concatenate(path):
    outfiles = [os.path.join(path, i) for i in os.listdir(path) if i.endswith('.out')]
    sort_key = lambda d : int(d.split('-')[-1].split('.')[0]) #sorted by num
    newfilelist = sorted(outfiles, key=sort_key) # file name sorted
    repeatline = set()
    resultname = os.path.join(path, 'concatenate.result')
    f0 = open(resultname, 'w')
    for i in newfilelist:
        f1 = open(i)
        for j in f1:
            if j not in repeatline:
                f0.write(j)
                repeatline.add(j)
            else :
                pass
        f1.close()
    f0.close()

if __name__ == '__main__':
    p = options.path
    o = options.resultname
    concatenate(p, o)

