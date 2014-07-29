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
    range = linspace(minv, maxv, n)
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

#    print fb_RA_00
#    print fb_RA_01
#    print fb_RA_11
#    print gatk_RA_00
#    print gatk_RA_01
#    print gatk_RA_11
#    print sb_RA_00
#    print sb_RA_01
#    print sb_RA_11
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

    print x, len(x)
#    print y_fb_11, len(y_fb_11)
#    print y_fb_01, len(y_fb_01)
#    print y_fb_00, len(y_fb_00)
#    print y_gatk_11, len(y_gatk_11)
#    print y_gatk_01, len(y_gatk_01)
#    print y_gatk_00, len(y_gatk_00)
#    print y_sb_11, len(y_sb_11)
#    print y_sb_01, len(y_sb_01)
#    print y_sb_00, len(y_sb_00)
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
#    ax.bar(fb_left, array(y_fb_00), bar_width, color='orange', \
#bottom=array(y_fb_11)+array(y_fb_01))

    ax.bar(fb_left+bar_width, array(y_gatk_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width, array(y_gatk_01), bar_width, color='yellow', \
bottom=array(y_gatk_11))
#    ax.bar(fb_left+bar_width, array(y_gatk_00), bar_width, color='orange', \
#bottom=array(y_gatk_11)+array(y_gatk_01))

    ax.bar(fb_left+bar_width+bar_width, array(y_sb_11), bar_width, color='coral')
    ax.bar(fb_left+bar_width+bar_width, array(y_sb_01), bar_width, color='yellow', \
bottom=array(y_sb_11))
#    ax.bar(fb_left+bar_width+bar_width, array(y_sb_00), bar_width, color='orange', \
#bottom=array(y_sb_11)+array(y_sb_01))
    ax.set_xticks(arange(0.1, len(x), 5)+(bar_width+bar_width+bar_width)/2)

    ax.set_xticklabels([i for i in itertools.islice(x, None, None, 5)])
    ax.set_xlabel('R/A')
    ax.set_ylabel('Number of SNPs')
    ax.set_title("Genotypes of three callers in defferent R/A")
    plt.savefig('test.pdf')





if __name__ == '__main__':
    import sys
    bar_plot(sys.argv[1], sys.argv[2], sys.argv[3])
