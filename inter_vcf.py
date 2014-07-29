#/usr/lib/python
from HTSeq import VCF_Reader
from optparse import OptionParser

msg_usage = 'usage: %prog [-F] fb_file [-G] gatk_file [-S] sb_file'
descr ='''get the common snp site from three different softwares.
record those common sites with freebayse describe info.
'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-F', '--fb', dest = 'fbfile',
                     help = 'Input the vcf file generated by freebays.')
optparser.add_option('-G', '--gatk', dest = 'gatkfile',
                     help = 'Input the vcf file generated by GATK.')
optparser.add_option('-S', '--sb', dest = 'sbfile',
                     help = 'Input the vcf file generated by samtools.')
options, args = optparser.parse_args()

def common_vcf(x, y, z):
    '''find the common site in those three files.
    write the common sites to a new file with
    freebayse information.
    '''
    vf1 = VCF_Reader(x)
    vf2 = VCF_Reader(y)
    vf3 = VCF_Reader(z)
    vf1_set = set()
    vf2_set = set()
    vf3_set = set()
    for i in vf1:
        vf1_set.add(i.chrom+'-'+str(i.pos.pos))
    for i in vf2:
        vf2_set.add(i.chrom+'-'+str(i.pos.pos))
    for i in vf3:
        vf3_set.add(i.chrom+'-'+str(i.pos.pos))
    common_set = vf1_set & vf2_set & vf3_set
    rf = open('common_sites.vcf', 'w')
    for i in open(x, 'r'):
        if not x.startswith('#')
            obj = i.split()[0] + '-' + i.split()[1]
            if obj in common_set:
                rf.write(i)
    rf.close()

if __name__ == '__main__':
    fb = options.fbfile
    gatk = options.gatkfile
    sb = options.sbfile
    common_vcf(fb, gatk, sb)

