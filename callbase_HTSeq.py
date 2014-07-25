#!/usr/lib/python
#-*- coding:utf8 -*-
from optparse import OptionParser
import itertools
from HTSeq import VCF_Reader
from HTSeq import BAM_Reader
from HTSeq import SAM_Reader

msg_usage = 'usage: %prog [-S] bamfile [-V] vcffile [-O] outputfile'
descr ='''This script is used to call base in each SNP site in vcffile.
Because vcf files from freebayes, GATK, specially samtools always contain
wrong information on AO and RO iterms. In order to evaluate Callers accuracy,
this script born. You should have a bam file and a vcf file which offer the
SNP sites. The bam file which is binary file of sam file. Currently, the bam
file must have been splited N caigar cause RNAseq data.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-B', '--bam', dest = 'bamname',
                     help = 'input the sam file name.')
optparser.add_option('-V', '--vcf', dest = 'vcfname',
                     help = 'input the vcf file name.')
optparser.add_option('-O', '--out', dest = 'output',
                     help = 'output file name.')
options, args = optparser.parse_args()


def dict_read(aln_obj):
    '''parse aligned object, then generate a dictionary containing
    coordinate and base. return this dictionary'''
    readlist = [i for i in aln_obj.read_as_aligned.seq]
    start = aln_obj.iv.start
    upper = 0
    for c in aln_obj.cigar:
        type = c.type #M, D, I, H
        len = c.size
        if type == 'M':
            upper += len
        if type == 'D':
            for i in range(len):
                readlist.insert(upper, None)
            upper += len
        if type == 'I':
            del readlist[upper:upper+len]
    read_dict = {}
    for i, j in enumerate(readlist):
        read_dict[i+start] = j
    return read_dict

def callbase(bamfile, vcffile, out):
    BF = BAM_Reader(bamfile) #bam file
    VF = VCF_Reader(vcffile) #vcf file
    RF = open(out, 'w') #result file
    RF.write('ref-name\t' + 'position\t' + 'RO\t' + 'A\t' + 'T\t' + 'C\t' \
+ 'G\t' + 'N\t' + 'None\t' + 'others\n')
#the first line of the RF, position is what in vcf file,1-based
    begin = 0
    for i in VF:
        vcf_refname = i.chrom
        vcf_position = i.pos.pos -1 # change vcf_position into 0-based
        refbase = i.ref
        print i.chrom + '-' + str(i.pos.pos)
        begin, At, Tt, Ct, Gt, Nt, nonet, othert = 0, 0, 0, 0, 0, 0, 0, 0
        for ofset, j in enumerate(itertools.islice(BF, begin, None)):
            bam_refname = j.iv.chrom
            bam_start = j.iv.start
            bam_end = j.iv.end - 1
            if (vcf_refname == bam_refname and vcf_position >= bam_start
                and vcf_position <= bam_end):
                begin += ofset
                readdict = dict_read(j) # j is aln object
                yourbase = readdict[vcf_position]
                if yourbase == 'A':
                    At += 1
                elif yourbase == 'T':
                    Tt += 1
                elif yourbase == 'C':
                    Ct += 1
                elif yourbase == 'G':
                    Gt += 1
                elif yourbase == 'N':
                    Nt += 1
                elif yourbase == None:
                    nonet += 1
                else:
                    othert += 1
        RF.write(vcf_refname + '\t' + str(vcf_position+1) + '\t' + refbase + '\t'\
 + str(At) + '\t' + str(Tt) + '\t' + str(Ct) + '\t' + \
 str(Gt) + '\t' + str(Nt) + '\t' + str(nonet) + '\t' + str(othert) + '\n')
    RF.close()

if __name__ == '__main__':
    import sys
    import datetime

    b = options.bamname
    v = options.vcfname
    o = options.output
    st = datetime.datetime.now()
    callbase(b, v, o)
    et = datetime.datetime.now()
    duration = et - st
    print duration
