# -*- coding: UTF-8 -*-
from optparse import OptionParser
import re
'''
Created on 2017年9月30日

@author: zhangxiaowei
'''
parser = OptionParser()

parser.add_option("-v", "--vcffile", dest="vcffile",help="used only when -d 0")# action="callback",type="string",callback=useoptionvalue_previous1,
"""
0: other format -> vcf
1: vcf -> other format
"""

(options, args) = parser.parse_args()
if __name__ == '__main__':
    inf=open(options.vcffile,'r')
    of=open(options.vcffile+"rewrite","w")
    for line in inf:
        if line[0]=="#":
            if line[:6]!="#CHROM":
                linelist=re.split(r'\s+',line.strip())
            else:
                print(line.strip,file=of)
        linelist=re.split(r'\s+',line.strip())
        if linelist[0]