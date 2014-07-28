#!/usr/lib/python

from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import getRAlist

msg_usage = 'usage: %prog [-R] real_file [-F] fb_file [-G] gatk_file [-S]
sb_file'
descr ='''Using the lists returned from getRAlist.py to construct plot which
compare the accuracy among the three Callers.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-R', '--real', dest = 'realfile',
                     help = 'Input the common file recording common snp sites
                     info.')
optparser.add_option('-F', '--fb', dest = 'fbfile',
                     help = 'Input the vcf file generated by freebays.')
optparser.add_option('-G', '--gatk', dest = 'gatkfile',
                     help = 'Input the vcf file generated by GATK.')
optparser.add_option('-S', '--sb', dest = 'sbfile',
                     help = 'Input the vcf file generated by samtools.')
options, args = optparser.parse_args()


def doline_plot(realfile, fbfile, gatkfile, sbfile):
    positions, real_R, real_A, fb_R, fb_A, gatk_R, gatk_A, sb_R, sb_A = \
getRAlist.listforplot(realfile, fbfile, gatkfile, sbfile)
    pos = np.array(positions)
    pos_max = max(pos)
    pos_min = min(pos)
    realR = np.array(real_R)
    realA = np.array(real_A)
    fbR = np.array(fb_R)
    fbA = np.array(fb_A)
    gatkR = np.array(gatk_R)
    gatkA = np.array(gatk_A)
    sbR = np.array(sb_R)
    sbA = np.array(sb_A)

    real_fb_R = realR - fbR
    real_gatk_R = realR - gatkR
    real_sb_R = realR - sbR
    real_fb_A = realA - fbA
    real_gatk_A = realA - gatkA
    real_sb_A = realA - sbA

    real_fb_R_max = max(real_fb_R)
    real_gatk_R_max = max(real_gatk_R)
    real_sb_R_max = max(real_sb_R)
    R_max = max(real_fb_R_max, real_gatk_R_max, real_sb_R_max)

    real_fb_A_max = max(real_fb_A)
    real_gatk_A_max = max(real_gatk_A)
    real_sb_A_max = max(real_sb_A)
    A_max = max(real_fb_A_max, real_gatk_A_max, real_sb_A_max)

    real_fb_R_min = min(real_fb_R)
    real_gatk_R_min = min(real_gatk_R)
    real_sb_R_min = min(real_sb_R)
    R_min = min(real_fb_R_min, real_gatk_R_min, real_sb_R_min)


    real_fb_A_min = min(real_fb_A)
    real_gatk_A_min = min(real_gatk_A)
    real_sb_A_min = min(real_sb_A)
    A_min = min(real_fb_A_min, real_gatk_A_min, real_sb_A_min)

    print 'pos_min: %s.'%pos_min  #195013
    print 'pos_max: %s.'%pos_max  #42354375
    print 'R_min: %s.'%R_min      #-1
    print 'R_max: %s.'%R_max      #1472
    print 'A_min: %s.'%A_min      #-10
    print 'A_max: %s.'%A_max      #7356

    plt.figure(1)
    plt.subplot(1,1,1)
    plt.title('R accuracy')
    plt.xlabel('SNP sites')
    plt.ylabel('realR - callersR')

    fblineR, = plt.plot(pos, real_fb_R, 'b-', linewidth=0.5)
    gatklineR, = plt.plot(pos, real_gatk_R, 'r-', linewidth=0.5),
    sblineR, = plt.plot(pos, real_sb_R, 'y-', linewidth=0.5)
    plt.legend([fblineR, gatklineR, sblineR], ['freebayes', 'GATK', 'samtools'])
    plt.axis([195000, 42400000, -1, 700])
    plt.savefig('R_accuracy.pdf')

    plt.figure(2)
    plt.subplot(1,1,1)
    plt.title('A accuracy')
    plt.xlabel('SNP sites')
    plt.ylabel('realA - callersA')

    fblineA, = plt.plot(pos, real_fb_A, 'b-', linewidth=0.5)
    gatklineA, = plt.plot(pos, real_gatk_A, 'r-', linewidth=0.5),
    sblineA, = plt.plot(pos, real_sb_A, 'y-', linewidth=0.5)
    plt.legend([fblineA, gatklineA, sblineA], ['freebayes', 'GATK', 'samtools'])
    plt.axis([195000, 42400000, -10, 2500])
    plt.savefig('A_accuracy.pdf')

if __name__ == '__main__':
    r = options.realfile
    f = options.fbfile
    g = options.gatkfile
    s = options.sbfile
    do_plot(r, f, g, s)

