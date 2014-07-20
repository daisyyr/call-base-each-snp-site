#!/usr/lib/python
#-*- coding:utf8 -*-
from optparse import OptionParser

msg_usage = 'usage: %prog [-S] samfile [-V] vcffile [-O] outputfile'
descr ='''This script is used to call base in each SNP site in vcffile.
Because vcf files from freebayes, GATK, specially samtools always contain
wrong information on AO and RO iterms. In order to evaluate Callers accuracy,
this script born. You should have a sam file and a vcf file which offer the
SNP sites. The sam file usually come from bam file which is binary file of sam
file. Currently, the bam file have been splited N caigar cause RNAseq data.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-S', '--sam', dest = 'samname',
                     help = 'input the sam file name.')
optparser.add_option('-V', '--vcf', dest = 'vcfname',
                     help = 'input the vcf file name.')
optparser.add_option('-O', '--out', dest = 'output',
                     help = 'output file name.')
options, args = optparser.parse_args()

def SNP_sites(vcffile):
    '''Get snp sites information from a vcf file, return a list contain those
    informations'''
    f0 = open(vcffile, 'r')
    global snp_site_list
    snp_site_list = []
    for i in f0:
        if i.startswith('#'):
            pass
        else:
            j = i.split()[0] + '-' + i.split()[1] + '-' + i.split()[3]
            snp_site_list.append(j)
    f0.close()
    return snp_site_list

def parseCIGAR(cigar):
    '''parse the cigar field and calculate the read truth length.
    return two dictionaries about cigar info.'''
    ligation = ''
    dic1 = {}
    dic2 = {}
    for i in cigar:
        if i.isdigit():
            ligation += i
        if i.isalpha():
            if i in dic1:
                dic1[i] = dic1[i] + '-' + ligation
                ligation = ''
            else:
                dic1[i] = ligation
                ligation = ''
    for i in cigar:
        if i.isdigit():
            ligation += i
        if i.isalpha():
            if i in dic2:
                dic2[i] += int(ligation)
                ligation = ''
            else:
                dic2[i] = int(ligation)
                ligation = ''
    return dic1,dic2

def parse_read(dic1, CIGARorder, readlist, start_site):
    '''assign each read base a cordinate number based on refsequence.
    return a dictionary {cordinate:base}
    '''
    countM = CIGARorder.count('M')
    countI = CIGARorder.count('I')
    countD = CIGARorder.count('D')
    Mn, In, Dn, upper = 0
    for i in CIGARorder:
        if i == 'M' and Mn < countM:
            upper += int(dic1[i].split('-')[Mn])
            Mn += 1
        if i == 'I' and In < countI:
            Ilens = int(dic1[i].split('-')[In])
            del readlist[upper:upper+Ilens]
            In += 1
        if i == 'D' and Dn < countD:
            Dlens = int(dic1[i].split('-')[Dn])
            for i in range(Dlens):
                readlist.insert(upper, 'None')
            upper += Dlens
            Dn +=1
    read_dict = {}
    for i,j in enumerate(readlist):
        read_dict[i+start_site] = j
    return read_dict

def callbase(samfile, out):
    '''samfile is converted by 'somtools view', so no header informations.
    '''
    f0 = open(samfile, 'r')
    f1 = open(out, 'w')
    f1.write('ref-name\t' + 'position\t' + 'ROi\t' + 'A\t' + 'T\t' + 'C\t' \
+ 'G\t' + 'other\n')
    for i in snp_site_list:
        k = i.split('-')
        O_refname = k[0] #O means snp site information
        O_position = int(k[1])
        refbase = k[2]
        print O_refname + '-' + str(O_position)
        Acount, Tcount, Ccount, Gcount, othercount = 0, 0, 0, 0, 0
        for j in f0:
            R_refname = j.split()[2]
            R_position = int(j.split()[3])
            if O_refname != R_refname:
                continue
            else:
                if O_position < R_position:
                    continue
                CIGAR = j.split()[5]
                dic1, dic2 = parseCIGAR(CIGAR) # dic1 is str, dic2 is digit
                real_readlen = dic2.get('M') + dic2.get('D', 0) # 'I' needn't to consider
                upper = R_position + real_readlen - 1
                if O_position > upper:
                    continue
                if O_position >= R_position and O_position <= upper:
                    CIGARorder = []
                    for i in CIGAR:
                        if i.isalpha() and i != 'H':
                            CIGARorder.append(i)
                    readlist = []
                    old_read = j.split()[9]
                    readlist.extend(old_read)
                    read_dict = parse_read(dic1, CIGARorder, readlist, R_position)
                    yourbase = read_dict[O_position]
                    if yourbase == 'A':
                        Acount += 1
                    elif yourbase == 'T':
                        Tcount += 1
                    elif yourbase == 'C':
                        Ccount += 1
                    elif yourbase == 'G':
                        Gcount += 1
                    else:
                        othercount += 1
        f1.write(O_refname + '\t' + refbase + '\t' + str(O_position) + '\t'\
 + str(Acount) + '\t' + str(Tcount) + '\t' + str(Ccount) + '\t' + \
 str(Gcount) + '\t' + str(othercount) + '\n')
        f0.seek(0,0)
    f0.close()
    f1.close()

if __name__ == '__main__':
    import sys
    s = options.samname
    v = options.vcfname
    o = options.output
    SNP_sites(v)
    callbase(s, o)
