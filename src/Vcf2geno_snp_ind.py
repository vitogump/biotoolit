# -*- coding: UTF-8 -*-
import re,sys
from optparse import OptionParser
from NGS.BasicUtil import Util, VCFutil
'''
Created on 2015-4-22

@author: liurui
'''
parser = OptionParser()
parser.add_option("-v", "--vcffile", dest="vcffilename",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--configure", dest="configure")
parser.add_option("-s","--software",dest="software",help="GATK or samtools ")
parser.add_option("-1", "--ld-window-kb", dest="ldwinkb")
parser.add_option("-2", "--ld-window", dest="ldwin")
parser.add_option("-d","--dilute",dest="dilute",default="1")
parser.add_option("-o","--outputpre",dest="outputpre")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

configure = open(options.configure, 'r')
if len(sys.argv)!=5:
    print("python Vcf2geno_snp_ind.py [vcf1] [outfilePrefix] [SOFTWARE] [perfix]")
    exit(-1)
VcfFileName = sys.argv[1]
OutFilePrefix = sys.argv[2]

if __name__ == '__main__':
    vcfdata=VCFutil.VCF_Data(options.vcffilename.strip())