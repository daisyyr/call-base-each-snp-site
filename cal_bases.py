#!/usr/lib/python
#-*- coding:utf8 -*-
import re

def parseCIGAR(cigar):
    '''parse the cigar field and calculate the read truth length.
    return the length in int.'''
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
    '''eg:dic1={'I': '1-1', 'M': '31-20-48'},
    CIGARorder=[M,I,M,I,M],
    read=['A','T','C','T','G','T','C','A','T','C','T','G','T'],
    start_site=2089
    '''
    countM = CIGARorder.count('M') #3
    countI = CIGARorder.count('I') #2
    countD = CIGARorder.count('D') #0
    Mn = 0
    In = 0
    Dn = 0
    upper = 0
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

def positions(vcffile):
    '''get snp sites information from a vcf file, return a list contain those
    informations'''
    f0 = open(vcffile, 'r')
    global snp_site_list
    snp_site_list = []
    for i in f0:
        if i.startswith('#'):
            pass
        else:
            j = '-'.join(i.split()[0:2])
            snp_site_list.append(j)
    f0.close()
    return snp_site_list

def findbase(samfile):
    '''samfile is converted by somtools view, so no header informations
    '''
    f0 = open(samfile, 'r')
    f1 = open('result.txt', 'w')
    f1.write('ref-name\t'+'position\t'+'A\t'+'T\t'+'C\t'+'G\t'+'other\n')
    for i in snp_site_list:
        O_refname = i.split('-')[0] #O means snp site information
        O_position = int(i.split('-')[1])
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
                        print yourbase
                        othercount += 1
        f1.write(O_refname+'\t'+str(O_position)+'\t'+str(Acount)+'\t'+str(Tcount)+'\t'+str(Ccount)+'\t'+str(Gcount)+'\t'+str(othercount)+'\n')
        f0.seek(0,0)
    f0.close()
    f1.close()

if __name__ == '__main__':
    import sys
    positions(str(sys.argv[1]))
    findbase(str(sys.argv[2]))


