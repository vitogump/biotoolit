# -*- coding: UTF-8 -*-
import re,sys,random,copy


if len(sys.argv)<4:
    print("python Vcf2Ped.py [vcf1] [outfilePrefix] [SOFTWARE] [dilutetodensity]")
    exit(-1)

VcfFileName = sys.argv[1]
OutFilePrefix = sys.argv[2]
vcffile = open(VcfFileName, "r")
mapfile = open(OutFilePrefix+".map", "w")
pedfile = open(OutFilePrefix+".ped", "w")
SOFTWARE=sys.argv[3].upper()


snpPositionlist=[]
pedmap={}
temppedmap={}

line = vcffile.readline()
totaltitles=1
while re.search(r"^##", line) != None:
    line = vcffile.readline()
    totaltitles+=1
if len(sys.argv)==5:
    dilutetodensity=float(sys.argv[4])


else:
    dilutetodensity=None

title=re.split(r"\s+",line.strip())
total_individ= len(title) -9
print(title,len(title),total_individ)
for outName in title[len(title)-total_individ:]:
    pedmap[outName]=[]
    temppedmap[outName]=[]
    tempsnpPositionlist=[]
currentChromSome=None
print("exclude those sites which ref is not N ,INDEL ,or multiple alleles")
totalrecsbeforelastchrom=None
startposoflastchrom=None
firstChrom=False

for line in vcffile:
    linelist = re.split(r"\s+",line.strip())
    
    if dilutetodensity !=None and currentChromSome!=linelist[0]:
        print(linelist[0])
        if   currentChromSome!=None and not firstChrom:
            #this situation,only the lastchrom is undiluted,all other is diluted
            tempsnpPositionlist=copy.deepcopy(snpPositionlist[:totalrecsbeforelastchrom])
            for i in range(total_individ):
                temppedmap[title[9+i]]=copy.deepcopy(pedmap[title[9+i]][:totalrecsbeforelastchrom*2])
            totallocsofcurchr=len(snpPositionlist)-totalrecsbeforelastchrom
            
            
            endposoflastchrom=int(snpPositionlist[-1][-1])
            if totallocsofcurchr!=0:
                dilute=dilutetodensity*(endposoflastchrom-startposoflastchrom)/(1000*totallocsofcurchr)
            else:
                dilute=0

            print(currentChromSome,linelist[0],dilutetodensity,"totallocsofcurchr",totallocsofcurchr,"startposoflastchrom:",startposoflastchrom,"endposoflastchrom",endposoflastchrom,"dilute:",dilute)
            startposoflastchrom=int(linelist[1])
            if dilute<1 and dilute>0:
                VcfRecRandomSelectIdxlist=random.sample([j for j in range(totallocsofcurchr)],int(dilute*totallocsofcurchr)+1)
                VcfRecRandomSelectIdxlist.sort()
                
                for selectedidx in VcfRecRandomSelectIdxlist:
                    print(snpPositionlist[totalrecsbeforelastchrom+selectedidx])
                    tempsnpPositionlist.append(copy.deepcopy(snpPositionlist[totalrecsbeforelastchrom+selectedidx]))
                    for i in range(total_individ):
#                         print(title[9+i],"totalrecsbeforelastchrom",totalrecsbeforelastchrom,selectedidx,"len(snpPositionlist)",len(snpPositionlist),"len(pedmap[title[9+i]])",len(pedmap[title[9+i]]),pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2],pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2+1])
                        temppedmap[title[9+i]]+=[pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2],pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2+1]]
            pedmap=copy.deepcopy(temppedmap)
            snpPositionlist=copy.deepcopy(tempsnpPositionlist)
            temppedmap={};tempsnpPositionlist=[]
            totalrecsbeforelastchrom=len(snpPositionlist)
#             mapfile = open(OutFilePrefix+".map", "w")
#             pedfile = open(OutFilePrefix+".ped", "w")
#             for elem in snpPositionlist:
#                 print(elem[0],elem[1],elem[2],elem[3],sep='\t',file=mapfile)
#             i=1
#             mapfile.close()
#             for name in sorted(pedmap.keys()):
#                 print(i,name,"0","0","1","1","\t".join(pedmap[name]),sep='\t',file=pedfile)
#                 i+=1
#             pedfile.close()
            currentChromSome=linelist[0]
        elif  currentChromSome==None:
            #first time
            currentChromSome=linelist[0]
            firstChrom=True
        elif  currentChromSome!=None and firstChrom:
            firstChrom=False
            tempsnpPositionlist=[]
            for i in range(total_individ):
                temppedmap[title[9+i]]=[]
            endposoflastchrom=int(snpPositionlist[-1][-1]);startposoflastchrom=int(snpPositionlist[0][-1])
            totallocsofcurchr=len(snpPositionlist)
            dilute=dilutetodensity*(endposoflastchrom-startposoflastchrom)/(1000*totallocsofcurchr)
            startposoflastchrom=int(linelist[1])
            print("totallocsofcurchr",totallocsofcurchr,"len(pedmap)",len(pedmap[title[9+i]]))
            if dilute<1 and dilute>0:
                VcfRecRandomSelectIdxlist=random.sample([j for j in range(totallocsofcurchr)],int(dilute*totallocsofcurchr)+1)
                VcfRecRandomSelectIdxlist.sort()
                for selectedidx in VcfRecRandomSelectIdxlist:
                    tempsnpPositionlist.append(copy.deepcopy(snpPositionlist[selectedidx]))
                    for i in range(total_individ):
                        temppedmap[title[9+i]]+=[pedmap[title[9+i]][selectedidx*2],pedmap[title[9+i]][selectedidx*2+1]]
                pedmap=copy.deepcopy(temppedmap)
                temppedmap={}
                snpPositionlist=copy.deepcopy(tempsnpPositionlist)
                tempsnpPositionlist=[]
            totalrecsbeforelastchrom=len(snpPositionlist)
            print("totalrecsbeforelastchrom",totalrecsbeforelastchrom,"should equal to half len(pedmap)",len(pedmap[title[9+i]]))
            currentChromSome=linelist[0]
            
    if linelist[3].strip().upper()=='N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip())>1:#when ref is N ,or INDEL ,or multiple allels 
        continue
    
    snpPositionlist.append((re.sub(r"[\D\W]+","",linelist[0]),re.sub(r"\..*$","",linelist[0])+"_"+linelist[1],0,linelist[1]))
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
else:
    if dilutetodensity !=None:
        #this situation,only the lastchrom is undiluted,all other is diluted
        tempsnpPositionlist=copy.deepcopy(snpPositionlist[:totalrecsbeforelastchrom])
        for i in range(total_individ):
            temppedmap[title[9+i]]=copy.deepcopy(pedmap[title[9+i]][:totalrecsbeforelastchrom*2])
        print("snpPositionlist",len(snpPositionlist),"totalrecsbeforelastchrom",totalrecsbeforelastchrom)
        totallocsofcurchr=len(snpPositionlist)-totalrecsbeforelastchrom
        
        endposoflastchrom=int(snpPositionlist[-1][-1])
        dilute=dilutetodensity*(endposoflastchrom-startposoflastchrom)/(1000*totallocsofcurchr)
        
        
        print("last",dilutetodensity,"totallocsofcurchr",totallocsofcurchr,"startposoflastchrom:",startposoflastchrom,"endposoflastchrom",endposoflastchrom,"dilute:",dilute)
        startposoflastchrom=int(linelist[1])
        if dilute<1 and dilute>0:
            VcfRecRandomSelectIdxlist=random.sample([j for j in range(totallocsofcurchr)],int(dilute*totallocsofcurchr)+1)
            VcfRecRandomSelectIdxlist.sort()
            print(VcfRecRandomSelectIdxlist)
            for selectedidx in VcfRecRandomSelectIdxlist:
                print("last",snpPositionlist[totalrecsbeforelastchrom+selectedidx])
                tempsnpPositionlist.append(copy.deepcopy(snpPositionlist[totalrecsbeforelastchrom+selectedidx]))
                for i in range(total_individ):
                    temppedmap[title[9+i]]+=[pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2],pedmap[title[9+i]][totalrecsbeforelastchrom*2+selectedidx*2+1]]
        pedmap=copy.deepcopy(temppedmap)
        snpPositionlist=copy.deepcopy(tempsnpPositionlist)
        temppedmap={};tempsnpPositionlist=[]
        totalrecsbeforelastchrom=len(snpPositionlist)
        print(len(pedmap[title[9]]),"should equal doulbe ",totalrecsbeforelastchrom)
        currentChromSome=linelist[0]
for elem in snpPositionlist:
    print(elem[0],elem[1],elem[2],elem[3],sep='\t',file=mapfile)
i=1
for name in sorted(pedmap.keys()):
    print(i,name,"0","0","1","1","\t".join(pedmap[name]),sep='\t',file=pedfile)
    i+=1       
mapfile.close()
pedfile.close()    
vcffile.close()
print("finish")
    
    