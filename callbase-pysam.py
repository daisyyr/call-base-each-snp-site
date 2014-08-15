#!/usr/lib/python
#-*- coding:utf-8 -*-
from optparse import OptionParser
from pysam import Samfile

msg_usage = 'usage: %prog [-B] bamfile [-S] SNPsites file [-O] outputfile'
descr ='''This script is used to call base in each SNP sites. Because
SNP sites in vcf files from freebayes, GATK, specially samtools always contain
wrong information on DP, AO and RO iterms. In order to evaluate Callers accuracy,
and for the accuracy of the next study about ASE,this script born. You should
have a bam file and a file containing SNP sites, each line of standard file
like this: 'Chr1\t65348', just like the first two columns in vcf files, you can
also use the vcf file which only the first two columns will be used.The bam file
which is binary file of sam file.
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
        self.Rbase = line.split()[3]
        self.Abase = line.split()[4]

def callbase(bamfile, snpsites, out):
    BF = Samfile(bamfile, 'rb') #open your bam file
    SF = open(snpsites, 'r')    #the file contain snp sites info
    RF = open(out, 'w')         #resulte file
    RF.write('ref_name\tpos\tRbase\tAbase\tA\tT\tC\tG\tN\tothers\n')
    for i in SF:
        if i.startswith('#'):
            continue
        else:
            line = ParseSNPsitesLine(i)
            vcf_pos = line.pos-1 #change 1-base to 0-based
            vcf_refname = line.chrom
            print 'processing: %s %s...'%(vcf_refname, str(vcf_pos))
            At, Tt, Ct, Gt, Nt, othert = 0, 0, 0, 0, 0, 0
            for i in BF.pileup(vcf_refname, vcf_pos, vcf_pos+1):
                if i.pos == vcf_pos:
                    vcf_Rbase = line.Rbase
                    vcf_Abase = line.Abase
                    for j in i.pileups:
                        yourbase = j.alignment.seq[j.qpos]
                        if yourbase == 'A': At += 1
                        elif yourbase == 'T': Tt += 1
                        elif yourbase == 'C': Ct += 1
                        elif yourbase == 'G': Gt += 1
                        elif yourbase == 'N': Nt += 1
                        else: othert += 1
        RF.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(vcf_refname, \
str(vcf_pos+1), vcf_Rbase, vcf_Abase, str(At), str(Tt), str(Ct), str(Gt), \
str(Nt), str(othert)))
    BF.close()

if __name__ == '__main__':
    import datetime
    b = options.bamname
    s = options.snpsites
    o = options.output
    st = datetime.datetime.now()
    callbase(b, s, o)
    et = datetime.datetime.now()
    duration = et - st
    print duration
