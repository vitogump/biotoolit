# -*- coding: UTF-8 -*-
import re,os,sys

if len(sys.argv)!=5:
    print("python Vcf2Ped.py [vcf1] [outfilePrefix] [SOFTWARE] [perfix]")
    exit(-1)

VcfFileName = sys.argv[1]
OutFilePrefix = sys.argv[2]
vcffile = open(VcfFileName, "r")
mapfile = open(OutFilePrefix+".map", "w")
pedfile = open(OutFilePrefix+".ped", "w")
SOFTWARE=sys.argv[3].upper()
positionlist=[]
pedmap={}

line = vcffile.readline()
while re.search(r"^##", line) != None:
    line = vcffile.readline()
title=re.split(r"\s+",line.strip())
total_individ= len(title) -9
print(title,len(title),total_individ)
for outName in title[len(title)-total_individ:]:
    pedmap[outName]=[]

currentChromSome=None
for line in vcffile:
    linelist = re.split(r"\s+",line)
    if linelist[3].strip().upper()=='N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip())>1:#when ref is N ,or INDEL ,or multiple allels 
        continue
    
    positionlist.append((linelist[0].replace("scaffold",""),linelist[0]+"_"+linelist[1],0,linelist[1]))
    if SOFTWARE=="GATK":
        GT_idx=(re.split(":",linelist[8])).index("GT")#gatk GT:AD:DP:GQ:PL
        PL_idx=(re.split(":",linelist[8])).index("PL")
        for i in range(total_individ):
            sample=linelist[i+9]
            if len(re.split(":",sample))==1 or re.split(":",sample)[GT_idx]=="./." or  len(re.split(r",",re.split(":",sample)[PL_idx]))!=3:# ./.
                pl="0,0,0"
            else:
                pl=re.split(":",sample)[PL_idx]
                genotype = re.split(":",sample)[GT_idx]                
            if pl!="0,0,0":
                a1=int(re.search(r"(\d)/(\d)",genotype).group(1))
                a2=int(re.search(r"(\d)/(\d)",genotype).group(2))
                alle1=linelist[a1+3].strip()
                alle2=linelist[a2+3].strip()
                pedmap[title[9+i]]+=[alle1,alle2]
            else:
                pedmap[title[9+i]]+=['0','0']
    
    elif SOFTWARE=="SAMTOOLS":
        pass
        for i in range(total_individ):
            Sample=linelist[i+9]
            genotype = re.search(r"([^:]+):([^:]+):([^:]+)",Sample.strip()).group(1)
            pl = re.search(r"([^:]+):([^:]+):([^:]+)",Sample.strip()).group(2)
            if pl!="0,0,0":
                a1=int(re.search(r"(\d)/(\d)",genotype).group(1))
                a2=int(re.search(r"(\d)/(\d)",genotype).group(2))
                alle1=linelist[a1+3].strip()
                alle2=linelist[a2+3].strip()
                pedmap[title[9+i]]+=[alle1,alle2]
            else:
                pedmap[title[9+i]]+=['0','0']
            
for elem in positionlist:
    print(elem[0],elem[1],elem[2],elem[3],sep='\t',file=mapfile)
i=1
for name in pedmap.keys():
    print(i,name,"0","0","1","1","\t".join(pedmap[name]),sep='\t',file=pedfile)
    i+=1       
    
    
    
    