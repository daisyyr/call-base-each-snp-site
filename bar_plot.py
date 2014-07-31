#!/usr/lib/python
import genotyperange
from numpy import *
import itertools
import matplotlib.pyplot as plt

def forbarplot(x, y, z, n=10):
    RdA, fbgt, gatkgt, sbgt = genotyperange.get_RA_gp(x, y, z)
    RdAlist = []        # reduce the accuracy, only retain one dot
    for i in RdA:
        RdAlist.append('%.1f'%i)
    minv, maxv = min(RdA),max(RdA)
    fb_RA_00 = {}
    fb_RA_01 = {}
    fb_RA_11 = {}
    gatk_RA_00 = {}
    gatk_RA_01 = {}
    gatk_RA_11 = {}
    sb_RA_00 = {}
    sb_RA_01 = {}
    sb_RA_11 = {}
    for RA, fb, gatk, sb in zip(RdAlist, fbgt, gatkgt, sbgt):
        fb_RA_00[RA] = 0
        fb_RA_01[RA] = 0
        fb_RA_11[RA] = 0
        gatk_RA_00[RA] = 0
        gatk_RA_01[RA] = 0
        gatk_RA_11[RA] = 0
        sb_RA_00[RA] = 0
        sb_RA_01[RA] = 0
        sb_RA_11[RA] = 0
    for RA, fb, gatk, sb in zip(RdAlist, fbgt, gatkgt, sbgt):
        if fb == '0/0':
            fb_RA_00[RA] += 1
        if fb == '0/1':
            fb_RA_01[RA] += 1
        if fb == '1/1':
            fb_RA_11[RA] += 1
        if gatk == '0/0':
            gatk_RA_00[RA] += 1
        if gatk == '0/1':
            gatk_RA_01[RA] += 1
        if gatk == '1/1':
            gatk_RA_11[RA] += 1
        if sb == '0/0':
            sb_RA_00[RA] += 1
        if sb == '0/1':
            sb_RA_01[RA] += 1
        if sb == '1/1':
            sb_RA_11[RA] += 1
# get the r/a:count dictionary of each callers' each genotype
    print 'fb_RA_00: %s'%fb_RA_00
    print 'fb_RA_01: %s'%fb_RA_01
    print 'fb_RA_11: %s'%fb_RA_11
    print 'gatk_RA_00: %s'%gatk_RA_00
    print 'gatk_RA_01: %s'%gatk_RA_01
    print 'gatk_RA_11: %s'%gatk_RA_11
    print 'sb_RA_00: %s'%sb_RA_00
    print 'sb_RA_01: %s'%sb_RA_01
    print 'sb_RA_11: %s'%sb_RA_11

    x = []
    y_fb_11, y_fb_01, y_fb_00 = [], [], []
    y_gatk_11, y_gatk_01, y_gatk_00 = [], [], []
    y_sb_11, y_sb_01, y_sb_00 = [], [], []
    for i in sorted(fb_RA_11):
        x.append(i)
        y_fb_11.append(fb_RA_11[i])
    for i in sorted(fb_RA_01):
        y_fb_01.append(fb_RA_01[i])
    for i in sorted(fb_RA_00):
        y_fb_00.append(fb_RA_00[i])
    for i in sorted(gatk_RA_11):
        y_gatk_11.append(gatk_RA_11[i])
    for i in sorted(gatk_RA_01):
        y_gatk_01.append(gatk_RA_01[i])
    for i in sorted(gatk_RA_00):
        y_gatk_00.append(gatk_RA_00[i])
    for i in sorted(sb_RA_11):
        y_sb_11.append(sb_RA_11[i])
    for i in sorted(sb_RA_01):
        y_sb_01.append(sb_RA_01[i])
    for i in sorted(sb_RA_00):
        y_sb_00.append(sb_RA_00[i])
#get list of R/A, each callers' each genotype.
    print 'R/A ratio list: '
    print x, 'the length of R/A ratio list: %s'%len(x)
    print 'freebayse 1/1 in defferent R/A ratio:'
    print y_fb_11, 'the length of freebayse 1/1: %s'%len(y_fb_11)
    print 'freebayse 1/1 in defferent R/A ratio:'
    print y_fb_01, 'the length of freebayse 0/1: %s'%len(y_fb_01)
    print 'freebayse 1/1 in defferent R/A ratio:'
    print y_fb_00, 'the length of freebayse 0/0: %s'%len(y_fb_00)
    print 'gatk 1/1 in defferent R/A ratio:'
    print y_gatk_11, 'the length of gatk 1/1: %s'%len(y_gatk_11)
    print 'gatk 1/1 in defferent R/A ratio:'
    print y_gatk_01, 'the length of gatk 0/1: %s'%len(y_gatk_01)
    print 'gatk 1/1 in defferent R/A ratio:'
    print y_gatk_00, 'the length of gatk 0/0: %s'%len(y_gatk_00)
    print 'samtools 1/1 in defferent R/A ratio:'
    print y_sb_11, 'the length of samtools 1/1: %s'%len(y_sb_11)
    print 'samtools 1/1 in defferent R/A ratio:'
    print y_sb_01, 'the length of samtools 0/1: %s'%len(y_sb_01)
    print 'samtools 1/1 in defferent R/A ratio:'
    print y_sb_00, 'the length of samtools 0/0: %s'%len(y_sb_00)

    return x, y_fb_11, y_fb_01, y_fb_00, y_gatk_11, y_gatk_01, y_gatk_00, \
y_sb_11, y_sb_01, y_sb_00

def bar_plot(x, y, z):
    x, y_fb_11, y_fb_01, y_fb_00, y_gatk_11, y_gatk_01, y_gatk_00, \
y_sb_11, y_sb_01, y_sb_00 = forbarplot(x, y, z)
    fig = plt.figure(1, figsize=(28, 12))
    ax = fig.add_subplot(111)
    bar_width = 0.2
    fb_left = (arange(0.1, len(x)))
    ax.bar(fb_left, array(y_fb_11), bar_width, color='coral')
    ax.bar(fb_left, array(y_fb_01), bar_width, color='yellow', \
bottom=array(y_fb_11))
    ax.bar(fb_left, array(y_fb_00), bar_width, color='orange', \
bottom=array(y_fb_11)+array(y_fb_01))

    ax.bar(fb_left+bar_width, array(y_gatk_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width, array(y_gatk_01), bar_width, color='yellow', \
bottom=array(y_gatk_11))
    ax.bar(fb_left+bar_width, array(y_gatk_00), bar_width, color='orange', \
bottom=array(y_gatk_11)+array(y_gatk_01))

    ax.bar(fb_left+bar_width+bar_width, array(y_sb_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width+bar_width, array(y_sb_01), bar_width, color='yellow', \
bottom=array(y_sb_11))
    ax.bar(fb_left+bar_width+bar_width, array(y_sb_00), bar_width, color='orange', \
bottom=array(y_sb_11)+array(y_sb_01))
    ax.set_xticks(arange(0.1, len(x), 5)+(bar_width+bar_width+bar_width)/2)

    ax.set_xticklabels([i for i in itertools.islice(x, None, None, 5)])
    ax.set_xlabel('R/A')
    ax.set_ylabel('Number of SNPs')
    ax.set_title("Genotypes of three callers in defferent R/A")
    plt.savefig('test.pdf')

#in odrer to show more details, replot by compress date which R/A>4
    print '#'*80
    rex, rey_fb_11, rey_fb_01, rey_fb_00 = [], [], [], []
    rey_gatk_11, rey_gatk_01, rey_gatk_00 = [], [], []
    rey_sb_11, rey_sb_01, rey_sb_00 = [], [], []

    for ra, fb11, fb01, fb00, ga11, ga01, ga00, sb11, sb01, sb00 in \
zip(x, y_fb_11, y_fb_01, y_fb_00, y_gatk_11, y_gatk_01, y_gatk_00, \
y_sb_11, y_sb_01, y_sb_00):
        if float(ra) <= 4:
            rex.append(ra)
            rey_fb_11.append(fb11)
            rey_fb_01.append(fb01)
            rey_fb_00.append(fb00)
            rey_gatk_11.append(ga11)
            rey_gatk_01.append(ga01)
            rey_gatk_00.append(ga00)
            rey_sb_11.append(sb11)
            rey_sb_01.append(sb01)
            rey_sb_00.append(sb00)
        if float(ra) > 4:
            rey_fb_11[-1] += fb11
            rey_fb_01[-1] += fb01
            rey_fb_00[-1] += fb00
            rey_gatk_11[-1] += ga11
            rey_gatk_01[-1] += ga01
            rey_gatk_00[-1] += ga00
            rey_sb_11[-1] += sb11
            rey_sb_01[-1] += sb01
            rey_sb_00[-1] += sb00
    rex[-1] = rex[-1]+'+'
    print 'rex: %s length of rex: %s'%(rex, len(rex))
    print 'rey_fb_11: %s length of rey_fb_11: %s'%(rey_fb_11, len(rey_fb_11))
    print 'rey_fb_01: %s length of rey_fb_01: %s'%(rey_fb_01, len(rey_fb_01))
    print 'rey_fb_00: %s length of rey_fb_00: %s'%(rey_fb_00, len(rey_fb_00))
    print 'rey_gatk_11: %s length of rey_fb_11: %s'%(rey_gatk_11,\
len(rey_gatk_11))
    print 'rey_gatk_01: %s length of rey_fb_01: %s'%(rey_gatk_01,\
len(rey_gatk_01))
    print 'rey_gatk_00: %s length of rey_fb_00: %s'%(rey_gatk_00,\
len(rey_gatk_00))
    print 'rey_sb_11: %s length of rey_fb_11: %s'%(rey_sb_11, len(rey_sb_11))
    print 'rey_sb_01: %s length of rey_fb_01: %s'%(rey_sb_01, len(rey_sb_01))
    print 'rey_sb_00: %s length of rey_fb_00: %s'%(rey_sb_00, len(rey_sb_00))
    fig2 = plt.figure(2, figsize=(28, 12))
    ax = fig2.add_subplot(111)
    bar_width = 0.25
    fb_left = (arange(0.9, len(rex)))
    ax.bar(fb_left, array(rey_fb_11), bar_width, color='coral', \
label='1/1')
    ax.bar(fb_left, array(rey_fb_01), bar_width, color='yellow', \
bottom=array(rey_fb_11), label='0/1')
    ax.bar(fb_left, array(rey_fb_00), bar_width, color='red', \
bottom=array(rey_fb_11)+array(rey_fb_01), label='0/0')
    ax.legend()

    ax.bar(fb_left+bar_width, array(rey_gatk_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width, array(rey_gatk_01), bar_width, color='yellow', \
bottom=array(rey_gatk_11))
    ax.bar(fb_left+bar_width, array(rey_gatk_00), bar_width, color='orange', \
bottom=array(rey_gatk_11)+array(rey_gatk_01))

    ax.bar(fb_left+bar_width+bar_width, array(rey_sb_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width+bar_width, array(rey_sb_01), bar_width, color='yellow', \
bottom=array(rey_sb_11))
    ax.bar(fb_left+bar_width+bar_width, array(rey_sb_00), bar_width, color='orange', \
bottom=array(rey_sb_11)+array(rey_sb_01))
    ax.set_xticks(arange(0.9, len(rex), 5)+(bar_width+bar_width+bar_width)/2)

    ax.set_xticklabels([i for i in itertools.islice(rex, None, None, 5)])
    ax.set_xlabel('R/A', fontsize=30)
    ax.set_ylabel('Number of SNPs', fontsize=30)
    ax.set_title("Genotypes of three callers in defferent R/A", fontsize=40)
    plt.savefig('replot.pdf')

if __name__ == '__main__':
    import sys
    bar_plot(sys.argv[1], sys.argv[2], sys.argv[3])
