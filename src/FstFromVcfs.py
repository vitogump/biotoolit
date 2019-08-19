# -*- coding: UTF-8 -*-
import re, numpy, sys
import biofunc
from itertools import combinations

if len(sys.argv) < 3:
    print("python FstFromVcf.py [vcf1] [vcf2] [vcf3]....")
    exit(-1)
# FstPopfile = open(sys.argv[1], 'r')
FstsMap={}
ComparePopfiles = []

def FstBetween2Pop(Pop1Vcf, Pop2Vcf):
    pop1 = open(Pop1Vcf, 'r')
    pop2 = open(Pop2Vcf, 'r')
    FstChromMap={}
    pop1line = pop1.readline()
    while re.search(r"^##", pop1line) != None:
        pop1line = pop1.readline()
    if re.search(r'^#', pop1line) != None:
        pop1lines = pop1.readline()
    else:
        print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
        exit(-1)
    
    pop2line = pop2.readline()
    while re.search(r"^##", pop2line) != None:
        pop2line = pop2.readline()
    if re.search(r'^#',pop2line)!=None:
        pop2lines=pop2.readline()
    else:
        print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
        exit(-1)
        
               
    pop1.close()
    pop2.close()
    return FstChromMap

if len(sys.argv)==3:
    FstBetween2Pop(sys.argv[1],sys.argv[2])
else:
    pass
# Fst=list(combinations(l2,2))
# for vcfName in sys.argv[2:]:
#     ComparePopfiles.append(open(vcfName, 'r'))
