# -*- coding: UTF-8 -*-
import pickle, re, sys
import biofunc
import numpy as np
'''
Created on 2013-10-15

@author: liurui
'''
if len(sys.argv) != 4:
    print(len(sys.argv))
    print("python FindCandidateSNP [contation] [donotcontain] [outfilename]")
    exit(-1)
else:
    print("\t".join(sys.argv))
ContainCandidatSNPFileName = sys.argv[1]  # mutation type aa
NoCandidatSNPFileName = sys.argv[2]  # wild type rr

outfilename=sys.argv[3]

def compareSNP(vcfContainCandidateSNP, vcfNoCandidateSNP, outfileName=None):
    if outfileName != None:
        outputfile = open(outfileName, 'w')
    else:
        outputfile = open("test.txt", 'w')
    try:
        noIndex = pickle.load(open(vcfNoCandidateSNP + ".myindex", 'rb'))

    except IOError:
        biofunc.indexVCF(vcfNoCandidateSNP, vcfNoCandidateSNP + ".myindex")

        noIndex = pickle.load(open(vcfNoCandidateSNP + ".myindex", 'rb'))

    print("index finished")
    numOfrr_indvds = len(noIndex["title"][9:])

        
    containF = open(vcfContainCandidateSNP, 'r')
    noF = open(vcfNoCandidateSNP, 'r')

    
    
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
 ####################    read the records of currentchrom in the  NoCandidatSNPFileName into a map,the value of the map is the list which record the line
        if currentChrom != lineClist[0]:
            currentChrom = lineClist[0]
            lineNmap = {}


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
                

 ##################################       
        snp_position = lineClist[1]
        
        if snp_position not in lineNmap:
            print(lineC.strip(), file=outputfile)

        lineC = containF.readline()
        
    outputfile.close()
    containF.close()
    noF.close()
        

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

if __name__ == '__main__':
    compareSNP(ContainCandidatSNPFileName, NoCandidatSNPFileName, outfilename)