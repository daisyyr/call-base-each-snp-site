#!/usr/lib/python
from optparse import OptionParser
from HTSeq import VCF_Reader

msg_usage = 'usage: %prog [-F] vcffile [-M] mode'
descr ='''Cause callbase.py have to pass bam file each line and bam
file always very big. I can't tackle this situation in mutiple threads by python,
so I wanna split vcf file first, then run the callbase.py with defferent vcf
file simutaneously.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-F', '--file', dest = 'filename',
                     help = 'input the vcf file name.')
optparser.add_option('-M', '--mode', dest = 'splitmode',
                     help = 'which mode you wanna choose,\
only bychrom supported now')
options, args = optparser.parse_args()

def split_by_chrom(vcffile):
    '''split your vcf file by chromsome name.
    '''
    chrom_set = set()
    VF = VCF_Reader(vcffile)
    for i in VF:
        chrom_set.add(i.chrom) #get the possible chrom names
    f0 = open(vcffile, 'r')
    for chr in chrom_set:
        f1 = open(chr, 'w')
        for i in f0:
            if i.startswith(chr+'\t'):
                f1.write(i)
        f1.close()
        f0.seek(0, 0)
    f0.close()

if __name__ == '__main__':
    f = options.filename
    m = options.splitmode
    if m == 'bychrom':
        split_by_chrom(f)
