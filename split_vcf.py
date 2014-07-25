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
                     help = "which mode do you wanna choose,bychrom or bynumberoffiles.\
                     \tbychrom: splited by differnt chromosome.\t\
bynumberoffiles: how many files you wanna to get.")
optparser.add_option('-N', '--numberoffiles', dest = 'numberoffiles',
                     type = 'int',
                     help = 'input how many files you wanna to get.')
options, args = optparser.parse_args()

def split_bychrom(vcffile):
    '''split your vcf file by chromsome name.
    '''
    chrom_set = set()
    VF = VCF_Reader(vcffile)
    for i in VF:
        chrom_set.add(i.chrom) #get the possible chrom names
    f0 = open(vcffile, 'r')
    for chr in chrom_set:
        f1 = open(chr+'.vcf', 'w')
        for i in f0:
            if i.startswith(chr+'\t'):
                f1.write(i)
        f1.close()
        f0.seek(0, 0)
    f0.close()

def split_bynumberoffiles(vcffile, n):
    '''n means the number of splited files.If total_lines
    in vcf file divide n without remainder, you will get
    n files.But if not, you will get n+1 files because you
    need one more file to save those remind lines
    '''
    F_prefix = '.'.join(vcffile.split('.')[0:-1])
    f0 = open(vcffile, 'r')
    vcflist = [i for i in f0 if not i.startswith('#')]
    f0.close()
    total_lines = len(vcflist)
    remainder = total_lines % n
    LIEF = total_lines / n #lines in each file
    if remainder == 0:
        for i in range(n):
            f1 = open(F_prefix+'-'+str(i+1)+'.vcf', 'w')
            start = i * LIEF
            end = (i+1) * LIEF
            f1.write(''.join(vcflist[start:end]))
            f1.close()
    else:
        for i in range(n):
            f1 = open(F_prefix+'-'+str(i+1)+'.vcf', 'w')
            start = i * LIEF
            end = (i+1) * LIEF
            f1.write(''.join(vcflist[start:end]))
            f1.close()
        f2 = open(F_prefix+'-'+str(n+1)+'.vcf', 'w')
        f2.write(''.join(vcflist[end:]))
        f2.close()

if __name__ == '__main__':
    f = options.filename
    m = options.splitmode
    n = options.numberoffiles
    if m == 'bychrom':
        split_bychrom(f)
    elif m == 'bynumberoffiles':
        split_bynumberoffiles(f, n)
    else:
        print 'please choose the split mode. -h for more help'
