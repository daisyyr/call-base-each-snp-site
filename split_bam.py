#!/usr/lib/python
import subprocess
'''
split bam files by @SN
'''
def parseheader(bamfiles):
    cmd = 'samtools view -H %s'%bamfiles
    headerobj = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    header = headerobj.communicate()[0]
    SNs = []
    for i in header.split():
        if i.startswith('SN'):
            SNs.append(i.split(':')[1])
    return SNs

def splitbam(bamfiles):
    newfileprefix = '.'.join(bamfiles.split('.')[0:-1])
    for SN in parseheader(bamfiles):
        newfile = newfileprefix +'.'+SN+'.bam'
        cmd = 'nohup samtools view -b %s %s > %s &'%(bamfiles, SN, newfile)
        subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    import sys
    splitbam(sys.argv[1])

