# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle

def indexVCF(VCFName, indexFileName):
    vcffile = open(VCFName, 'r')
    vcfChromIndex = {}
    line = vcffile.readline()
    
    while re.search(r'^##', line) != None:
        line = vcffile.readline()
    
    if re.search(r'^#', line) != None:
        vcfChromIndex["title"] = re.split(r'\s+', line)
    else:
        print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
        exit(-1)        
    currentChrom = None
    lastPosition = vcffile.tell()
    print(line)
    line = vcffile.readline()
    print(line)
    while line:      
        linelist = re.split(r"\s+", line)
        if currentChrom != linelist[0]:
            currentChrom = linelist[0]
            vcfChromIndex[currentChrom] = lastPosition
            print("test", lastPosition)
        lastPosition = vcffile.tell()

        line = vcffile.readline()
    pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
    vcffile.close()


def getVcfMap(vcfFileName):
    """
    this func is from bio\test\posAroundGene\func.py ,and did some improvement,that is add  INFO = collist[7],and add INFO into a 
    read the vcffile into a map which keys are chrom,values are a list of tuple
    {chrNo:[(pos,REF,ALT),(pos,REF,ALT),,,,,],chrNo:[],,,,,,},the order of the tuples in the list,is according pos,
    you we can search a record by  binary chop search
    """
    vcfMap = {}
    vcfFile = open(vcfFileName, 'r')
    
    line = vcfFile.readline()
    while re.search(r'^##', line) != None:
#        print(line)
        line = vcfFile.readline()
    if re.search(r'^#', line) != None:
        lineslist = vcfFile.readlines()
    else:
        print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'\n" + line)
        exit(-1)
#    print("pass")
    currentLine = 0
    totalRecs = len(lineslist)
    while currentLine != totalRecs:
#        print(currentLine)
        collist = re.split(r'\s+', lineslist[currentLine])
        chrom = collist[0]
        pos = int(collist[1])
        REF = collist[3]
        ALT = collist[4]
        INFO = collist[7]
        if chrom in vcfMap:
            vcfMap[chrom].append((pos, REF, ALT, INFO))
        else:
            vcfMap[chrom] = [(pos, REF, ALT, INFO)]
        currentLine += 1
    vcfFile.close()
    return vcfMap


