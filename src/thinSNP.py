'''
Created on 2015-6-3

@author: liurui
'''
import re, sys, random, copy


if len(sys.argv)!=4:
    print("python thinmapped.py  [filePrefix]  [dilutetodensity] [1(map/ped) 2(geno/snp) 3(randomextractline)]")
    exit(-1)
dilutetodensity=float(sys.argv[2])
print(dilutetodensity)

if __name__ == '__main__':
    if sys.argv[3]==1:
        mapfile=open(sys.argv[1].strip()+".map",'r')
        pedfile=open(sys.argv[1].strip()+".ped",'r')
        outmapfile=open(sys.argv[1].strip()+"dilutetodensity"+sys.argv[2]+".map",'w')
        outpedfile=open(sys.argv[1].strip()+"dilutetodensity"+sys.argv[2]+".ped",'w')
        pedlist=[]
        for line in pedfile:
            linelist=re.split(r"\s+",line.strip())
            pedlist.append(linelist)
        curchrom=None;startpos=None
        for line in mapfile:
            maplinelist=re.split(r"\s+",line.strip())
    elif sys.argv[3]=="2":#geno snp
        genofile=open(sys.argv[1].strip()+".geno",'r')
        snpfile=open(sys.argv[1].strip()+".snp",'r')
        outgenofile=open(sys.argv[1].strip()+"dilutetodensity"+sys.argv[2]+".geno",'w')
        outsnpfile=open(sys.argv[1].strip()+"dilutetodensity"+sys.argv[2]+".snp",'w')
        snplinesMapbyChrom={}
        genolinesMapbyChrom={}
        chromorder=[]
        genolines=[]
        snplines=[]
        curchrom=None;startpos=None
        for snpline in snpfile:
            genoline=genofile.readline()
            snplinelist=re.split(r"\s+",snpline.strip())
            if curchrom!=snplinelist[1]:
                if curchrom!=None:
                    snplinesMapbyChrom[curchrom]=[]
                    genolinesMapbyChrom[curchrom]=[]
                    endpos=int(re.split(r"\s+",snplines[-1].strip())[3])
                    
                    dilute = dilutetodensity * (endpos - startpos) / (1000 * len(snplines))
                    print(dilutetodensity,endpos,startpos,dilute)
                    VcfRecRandomSelectIdxlist = random.sample([j for j in range(len(snplines))], int(dilute * len(snplines)))
                    VcfRecRandomSelectIdxlist.sort()
                    print("VcfRecRandomSelectIdxlist",VcfRecRandomSelectIdxlist)
                    if len(VcfRecRandomSelectIdxlist)<=1:
                        VcfRecRandomSelectIdxlist=[int(len(snplines)/2)]
                    for idx in VcfRecRandomSelectIdxlist:
                        snplinesMapbyChrom[curchrom].append(copy.deepcopy(snplines[idx]))
                        genolinesMapbyChrom[curchrom].append(copy.deepcopy(genolines[idx]))
                    startpos=int(snplinelist[3])
                elif curchrom==None and startpos==None:
                    #only the first snpline
                    startpos=int(snplinelist[3])
                genolines=[genoline]
                snplines=[snpline]
                
                curchrom=snplinelist[1]
                chromorder.append(curchrom)
                print(curchrom)
            else:
                genolines.append(genoline)
                snplines.append(snpline)
        else:
            print("final",curchrom,snplines[-1],startpos)
            endpos=int(re.split(r"\s+",snplines[-1].strip())[3])
            snplinesMapbyChrom[curchrom]=[]
            genolinesMapbyChrom[curchrom]=[]
            dilute = dilutetodensity * (endpos - startpos) / (1000 * len(snplines))
            VcfRecRandomSelectIdxlist = random.sample([j for j in range(len(snplines))], int(dilute * len(snplines)))
            VcfRecRandomSelectIdxlist.sort()
            for idx in VcfRecRandomSelectIdxlist:
                print(snplines[idx])
                print(genolines[idx])
                snplinesMapbyChrom[curchrom].append(copy.deepcopy(snplines[idx]))
                genolinesMapbyChrom[curchrom].append(copy.deepcopy(genolines[idx]))
        ###### write out file#################
        for chrom in chromorder:
            for rec in snplinesMapbyChrom[chrom]:
                print(rec,end="",file=outsnpfile)
            for rec in genolinesMapbyChrom[chrom]:
                print(rec,end="",file=outgenofile)
        outgenofile.close()
        genofile.close()
        snpfile.close()
        outsnpfile.close()
    elif sys.argv[3]=="3":#geno snp
        outfile=open(sys.argv[1].strip()+"dilutetodensity"+sys.argv[2],'w')
        f=open(sys.argv[1].strip(),'r')
        fcontext=f.readlines()
        VcfRecRandomSelectIdxlist = random.sample([j for j in range(len(fcontext))], int(dilutetodensity * len(fcontext)))
        VcfRecRandomSelectIdxlist.sort()
        VcfRecRandomSelectIdxlist[0]=0;VcfRecRandomSelectIdxlist[1]=1
        for idx in VcfRecRandomSelectIdxlist:
            print(fcontext[idx],end="",file=outfile)
        outfile.close()
        f.close()
            