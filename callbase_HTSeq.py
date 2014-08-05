#!/usr/lib/python
#-*- coding:utf8 -*-
from optparse import OptionParser
import itertools
from HTSeq import VCF_Reader
from HTSeq import BAM_Reader

msg_usage = 'usage: %prog [-B] bamfile [-S] SNPsites file [-O] outputfile'
descr ='''This script is used to call base in each SNP sites. Because
SNP sites in vcf files from freebayes, GATK, specially samtools always contain
wrong information on DP, AO and RO iterms. In order to evaluate Callers accuracy,
and for the accuracy of the next study about ASE,this script born. You should
have a bam file and a file containing SNP sites, each line of standard file
like this: 'Chr1\t65348', just like the first two columns in vcf files, you can
also use the vcf file which only the first two columns will be used.The bam file
which is binary file of sam file. Currently, the bamfile must have been splited N
caigar cause RNAseq data.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-B', '--bam', dest = 'bamname',
                     help = 'input the bam file name.')
optparser.add_option('-S', '--SNPsites', dest = 'snpsites',
                     help = 'input the vcf file name.')
optparser.add_option('-O', '--out', dest = 'output',
                     help = 'output file name.')
options, args = optparser.parse_args()

class ParseSNPsitesLine():
    def __init__(self, line):
        self.chrom = line.split()[0]
        self.pos = int(line.split()[1])

def dict_read(aln_obj):
    '''parse aligned object, then generate a dictionary containing
    coordinate, 0-based, and base. return this dictionary'''
    readlist = [i for i in aln_obj.read_as_aligned.seq]
    start = aln_obj.iv.start
    upper = 0
    for c in aln_obj.cigar:
        type = c.type #M, D, I, H
        len = c.size
        if type == 'M': upper += len
        elif type == 'D':
            for i in range(len):
                readlist.insert(upper, None)
            upper += len
        if type == 'I': del readlist[upper:upper+len]
    read_dict = {}
    for i, j in enumerate(readlist):
        read_dict[i+start] = j
    return read_dict

def callbase(bamfile, snpsites, out):
    BF = BAM_Reader(bamfile) #Bam File
    SF = open(snpsites, 'r') #SNPsites File
    RF = open(out, 'w') #Result File
    RF.write('ref-name\t' + 'position\t' + 'A\t' + 'T\t' + 'C\t' \
+ 'G\t' + 'N\t' + 'None\t' + 'others\n')
#the first line of the RF, position is what in vcf file,1-based
    for i in SF:
        line = ParseSNPsitesLine(i)
        vcf_position = line.pos -1 # change vcf_position into 0-based
        vcf_refname = line.chrom
        print 'processing: %s......'%(vcf_refname + '  ' + str(vcf_position+1))
        At, Tt, Ct, Gt, Nt, nonet, othert = 0, 0, 0, 0, 0, 0, 0
        for j in BF:
            bam_refname = j.iv.chrom
            bam_start = j.iv.start
            bam_end = j.iv.end - 1
            if (vcf_refname == bam_refname and bam_start <= vcf_position <= bam_end):
                readdict = dict_read(j) # j is aln object
                yourbase = readdict[vcf_position]
                if yourbase == 'A': At += 1
                elif yourbase == 'T': Tt += 1
                elif yourbase == 'C': Ct += 1
                elif yourbase == 'G': Gt += 1
                elif yourbase == 'N': Nt += 1
                elif yourbase == None: nonet += 1
                else: othert += 1
        RF.write(vcf_refname + '\t' + str(vcf_position+1) + '\t' + '\t'\
 + str(At) + '\t' + str(Tt) + '\t' + str(Ct) + '\t' + \
 str(Gt) + '\t' + str(Nt) + '\t' + str(nonet) + '\t' + str(othert) + '\n')
    RF.close()

if __name__ == '__main__':
#    import datetime
    b = options.bamname
    s = options.snpsites
    o = options.output
#    st = datetime.datetime.now()
    callbase(b, s, o)
#    et = datetime.datetime.now()
#    duration = et - st
#    print duration
