'''
Created on 2014-11-25

@author: liurui
'''
from optparse import OptionParser
import os


parser = OptionParser()
parser.add_option("-1", "--fastq1", dest="fastq1",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="fastq file reads 1")
parser.add_option("-2", "--fastq2", dest="fastq2",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="fastq file reads 2")

parser.add_option("-p","--partnumber",dest="partnumber")
parser.add_option("-o","--outputpath",dest="outputpath")
# (options, args) = parser.parse_args()

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
fq1=open(options.fastq1,'r')
fq2=open(options.fastq2,'r')
partnumber=int(options.partnumber)
outputpath=options.outputpath
writefile1_list=[]
writefile2_list=[]
if __name__ == '__main__':
    for i in range(partnumber):
        writefile1_list.append(open(os.path.join(outputpath, '%s_part%d' % (options.fastq1,i)), 'w'))
        writefile2_list.append(open(os.path.join(outputpath, '%s_part%d' % (options.fastq2,i)), 'w'))
    n=0
    for headerline1 in fq1:
        f_current1=writefile1_list[n%partnumber]
        f_current2=writefile2_list[n%partnumber]
        headerline2=fq2.readline()
        f_current1.write(headerline1)
        f_current2.write(headerline2)
        seqline1=fq1.readline()
        seqline2=fq2.readline()
        plusline1=fq1.readline()
        plusline2=fq2.readline()
        qualitline1=fq1.readline()
        qualitline2=fq2.readline()
        f_current1.write('%s%s%s'%(seqline1,plusline1,qualitline1))
        f_current2.write('%s%s%s'%(seqline2,plusline2,qualitline2))
        n=n+1
    for i in range(partnumber):
        writefile1_list[i].close()
        writefile2_list[i].close()
    fq1.close()
    fq2.close()
        
        
        