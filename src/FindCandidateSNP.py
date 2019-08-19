# -*- coding: UTF-8 -*-
import pickle, re, sys
import biofunc
import numpy as np
from scipy.stats import chisquare as chisq


if len(sys.argv) != 6:
    print(len(sys.argv))
    print("python FindCandidateSNP [pooling data(p/P) or individuals data(I/i)] [altnatehomVCF(aa)] [refhomVCF(rr)] [heterozygoteVCF(ra)(optional)] [outfilename]")
    exit(-1)
else:
    print("\t".join(sys.argv))
pool_individuals = sys.argv[1].lower()
ContainCandidatSNPFileName = sys.argv[2]  # mutation type aa
NoCandidatSNPFileName = sys.argv[3]  # wild type rr
HalfCandidatSNPFileName = None
if len(sys.argv) == 6:
    HalfCandidatSNPFileName = sys.argv[4]

outfilename=sys.argv[5]
indelfile = open(outfilename+"indel", 'w')

def compareSNP(vcfContainCandidateSNP, vcfNoCandidateSNP, halfInVCF, outfileName=None):
    if outfileName != None:
        outputfile = open(outfileName, 'w')
    else:
        outputfile = open("test.txt", 'w')
    try:
        noIndex = pickle.load(open(vcfNoCandidateSNP + ".myindex", 'rb'))
        halfIndex = pickle.load(open(halfInVCF + ".myindex", 'rb'))
    except IOError:
        biofunc.indexVCF(vcfNoCandidateSNP, vcfNoCandidateSNP + ".myindex")
        biofunc.indexVCF(halfInVCF, halfInVCF + ".myindex")
        noIndex = pickle.load(open(vcfNoCandidateSNP + ".myindex", 'rb'))
        halfIndex = pickle.load(open(halfInVCF + ".myindex", 'rb'))
    print("index finished")
    numOfrr_indvds = len(noIndex["title"][9:])
    numOfra_indvds = len(halfIndex["title"][9:])
        
    containF = open(vcfContainCandidateSNP, 'r')
    noF = open(vcfNoCandidateSNP, 'r')
    halfF = open(halfInVCF, 'r')
    
    
# 这个方法超级笨啊！
    currentChrom = None
    lineC = containF.readline()
    
    while lineC:
        if re.search(r'^##', lineC) != None:
            lineC = containF.readline()
            continue
        elif re.search(r'^#', lineC) != None:
            lineC = containF.readline()
            lineClist = re.split(r'\s+', lineC)
            numOfaa_indvds = len(lineClist[9:])
            continue
        lineClist = re.split(r'\s+', lineC.strip())
 ####################    read the records of currentchrom in the  NoCandidatSNPFileName,HalfCandidatSNPFileName into a map,the value of the map is the list which record the line
        if currentChrom != lineClist[0]:
            currentChrom = lineClist[0]
            lineNmap = {}
            lineHmap = {}

            try:            
                noF.seek(noIndex[currentChrom])
                lineN = noF.readline()
            except KeyError:
                print(currentChrom + "didn't find in " + NoCandidatSNPFileName)
                lineN = None  # skip the while loop
            while lineN and (re.split(r'\s+', lineN.strip()))[0] == currentChrom:
                lineNlist = re.split(r'\s+', lineN.strip())
                lineNmap[lineNlist[1]] = lineNlist
                lineN = noF.readline()
                
            try:   
                halfF.seek(halfIndex[currentChrom])    
                lineH = halfF.readline()
            except KeyError:
                print(currentChrom + "didn't find in " + HalfCandidatSNPFileName)
                lineH = None  # skip the while loop
            while lineH and (re.split(r'\s+', lineH.strip()))[0] == currentChrom:
                lineHlist = re.split(r'\s+', lineH.strip())
                lineHmap[lineHlist[1]] = lineHlist
                lineH = halfF.readline()
 ##################################       
        snp_position = lineClist[1]
        
        if snp_position not in lineNmap:
            if pool_individuals == 'p':
                if (snp_position in lineHmap) and lineClist[4] !=lineHmap[snp_position][4]:
                    break
                refalle_count = int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineClist[7]).group(1)) + int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineClist[7]).group(2))
                altalle_count = int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineClist[7]).group(3)) + int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineClist[7]).group(4))
                ob = np.array([refalle_count, altalle_count])
                ex = np.array([0.5, 0.5]) * np.sum(ob)
                
                if chisq(ob, ex, 1)[0] > 3.84 and refalle_count < altalle_count :            
                    if snp_position in lineHmap:
                        # chi-test
                        refalle_count = int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineHmap[snp_position][7]).group(1)) + int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineHmap[snp_position][7]).group(2))
                        altalle_count = int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineHmap[snp_position][7]).group(3)) + int(re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', lineHmap[snp_position][7]).group(4))
                        ob = np.array([refalle_count, altalle_count])
                        ex = np.array([0.5, 0.5]) * np.sum(ob)
                        value = chisq(ob, ex, 1)[0]
                        if value < 3.84:  # if degree of freedom is 1 the chisquare shuold be corrected.i don't whither the scipy do it or not
                            if lineClist[3].strip().upper() == 'N' or len(lineClist[3].strip()) > 1 or len(lineClist[4].strip()) > 1:
                                print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], snp_position, lineHmap[snp_position][3], lineHmap[snp_position][4], lineHmap[snp_position][7], value, sep='\t', file=indelfile)
                            else:
                                print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], snp_position, lineHmap[snp_position][3], lineHmap[snp_position][4], lineHmap[snp_position][7], value, sep='\t', file=outputfile)
                    else:
                        if lineClist[3].strip().upper() == 'N' or len(lineClist[3].strip()) > 1 or len(lineClist[4].strip()) > 1:
                            print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], 0, sep='\t', file=indelfile)
                        else:
                            print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], 0, sep='\t', file=outputfile)
            elif pool_individuals == 'i':
                if (snp_position in lineHmap) and lineClist[4] !=lineHmap[snp_position][4]:
                    break
                for sample in lineClist[9:]:
#                     print(lineClist,sample)
                    genotype = re.search(r"([^:]+):([^:]+):([^:]+)", sample.strip()).group(1)
                    pl = re.search(r"([^:]+):([^:]+):([^:]+)", sample.strip()).group(2)
                    if pl != "0,0,0":
                        a1 = int(re.search(r"(\d)/(\d)", genotype).group(1))
                        a2 = int(re.search(r"(\d)/(\d)", genotype).group(2))
                        if a1 != 1 or a2 != 1:
                            break  # if thus,then abandon this snp.go to the next line immediately without print this snp into outputfile
                    else:
                        continue
                else:  #
                    if snp_position in lineHmap:
                        for sample in lineHmap[snp_position][9:]:
                            genotype = re.search(r"([^:]+):([^:]+):([^:]+)", sample.strip()).group(1)
                            pl = re.search(r"([^:]+):([^:]+):([^:]+)", sample.strip()).group(2)
                            if pl != "0,0,0":
                                a1 = int(re.search(r"(\d)/(\d)", genotype).group(1))
                                a2 = int(re.search(r"(\d)/(\d)", genotype).group(2))
                                if a1 == a2:
                                    print()
                                    break  # if thus,then abandon this snp.go to the next line immediately without print this snp into outputfile
                            else:
                                continue
                        else:
                            if lineClist[3].strip().upper() == 'N' or len(lineClist[3].strip()) > 1 or len(lineClist[4].strip()) > 1:
                                print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], snp_position, lineHmap[snp_position][3], lineHmap[snp_position][4], lineHmap[snp_position][7], value, sep='\t', file=indelfile)
                            else:
                                print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], snp_position, lineHmap[snp_position][3], lineHmap[snp_position][4], lineHmap[snp_position][7], value, sep='\t', file=outputfile)
                    else:
                        if lineClist[3].strip().upper() == 'N' or len(lineClist[3].strip()) > 1 or len(lineClist[4].strip()) > 1:
                            print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], 0, sep='\t', file=indelfile)
                        else:
                            print(currentChrom, lineClist[1], lineClist[2], lineClist[3], lineClist[4], lineClist[5], lineClist[6], lineClist[7], 0, sep='\t', file=outputfile)
        lineC = containF.readline()
        

            
            
        
#         lineN=noF.readline()
#         while lineN:
#             lineNlist=re.split(r'\s+',lineN)
#             
#             if snp_position == lineNlist[1]:
#                 break
#             elif currentChrom != lineNlist[0]:#find snp exist in vcfContainCandidateSNP file but not in vcfNoCandidateSNP.now looking for halfInVCF file
#                print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],0,sep='\t',file=open("test.txt",'a')) 
# #                 halfF.seek(halfIndex[currentChrom])
# #                 for lineH in halfF:
# #                     lineHlist=re.split(r'\s+',lineH)
# #                     if snp_position == lineHlist[1]:
# #                         print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],lineHlist[3],lineHlist[4],lineHlist[7],sep='\t',file=open("test.txt",'a'))
# #                         break# also found in halfInVCF
# #                     elif currentChrom != lineHlist[0]:
# #                         print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],0,sep='\t',file=open("test.txt",'a'))
# #                         break#candidateSNP has been found but didn't found in halfInVCF
#             lineN=noF.readline()
#         else:#find snp exist in vcfContainCandidateSNP file but not in vcfNoCandidateSNP.now looking for halfInVCF file
#             print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],0,sep='\t',file=open("test.txt",'a'))
# #             halfF.seek(halfIndex[currentChrom])
# #             for lineH in halfF:
# #                 lineHlist=re.split(r'\s+',lineH)
# #                 if snp_position == lineHlist[1]:
# #                     print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],lineHlist[3],lineHlist[4],lineHlist[7],sep='\t',file=open("test.txt",'a'))
# #                     break#also found in halfInVCF
# #                 elif currentChrom != lineHlist[0]:
# #                     print(currentChrom,lineNlist[1],lineNlist[2],lineNlist[3],lineNlist[4],lineNlist[5],lineNlist[6],lineNlist[7],0,sep='\t',file=open("test.txt",'a'))
# #                     break#candidateSNP has been found but didn't found in halfInVCF
#         lineC=containF.readline()
    outputfile.close()
    containF.close()
    noF.close()
    halfF.close()
    indelfile.close() 
        

# indexVCF(CandidatSNPExistFileName,CandidatSNPExistFileName+".myindex")
def gethomozygoteSNP(inputvcfFileName, outputvcfFileName):
    vcfin = open(inputvcfFileName, 'r')
    vcfout = open(outputvcfFileName, 'w')
    line = vcfin.readline()
    
    while re.search(r'^##', line) != None:
        print(line, file=vcfout)
        line = vcfin.readline()
    print(line, end="", file=vcfout)
    
    for line in vcfin:
        linelist = re.split(r"\s+", line.strip())
        if linelist[3].strip().upper() == 'N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip()) > 1:
            continue
        genotype = re.search(r"([^:]+):([^:]+):([^:]+)", linelist[9].strip()).group(1)
        if genotype == "1/1" or genotype == "0/0":
            print(line, end="", file=vcfout)
    vcfin.close()
    vcfout.close()

def getheterozygoteSNP(inputvcfFileName, outputvcfFileName):
    vcfin = open(inputvcfFileName, 'r')
    vcfout = open(outputvcfFileName, 'w')
    line = vcfin.readline()
    
    while re.search(r'^##', line) != None:
        print(line, file=vcfout)
        line = vcfin.readline()
    print(line, end="", file=vcfout)
    
    for line in vcfin:
        linelist = re.split(r"\s+", line.strip())
        if linelist[3].strip().upper() == 'N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip()) > 1:
            continue
        genotype = re.search(r"([^:]+):([^:]+):([^:]+)", linelist[9].strip()).group(1)
        if genotype == "0/1" or genotype == "1/0":
            print(line, end="", file=vcfout)
    vcfin.close()
    vcfout.close()

# gethomozygoteSNP(ContainCandidatSNPFileName,"aa.vcf")
# gethomozygoteSNP(NoCandidatSNPFileName,"rr.vcf")
# getheterozygoteSNP(HalfCandidatSNPFileName,"ra.vcf")
         
# compareSNP("aa.vcf","rr.vcf","ra.vcf","afterchitest.myvcf")

compareSNP(ContainCandidatSNPFileName, NoCandidatSNPFileName, HalfCandidatSNPFileName, outfilename)    
    
    
    
    
    
