# -*- coding: UTF-8 -*-
import re
import sys
import numpy

VcfFileName = sys.argv[1]
FileOutName = sys.argv[2]
vcffile = open(VcfFileName, "r")
fileout = open(FileOutName, "w")
#testfile = open("testfile.txt","w")
windowWidth = 40000
skipover = 20000
windowNo = 1

CNMI = 0;CNMA = 0;HETEROZY = 0
line = vcffile.readline()



while re.search(r"^##", line) != None:
    line = vcffile.readline()
#skip over the annotate lines in the begining of the file

HeterozyRec = {}

winStartPos = 0
#
#while int(collist[1]) - winStartPos>windowWidth:
#    HeterozyRec[chrom].append(0)
#    windowNo+=1
##                        print(pos,currentLine,chrom,len(HeterozyRec[chrom]),"============================================")
#    winStartPos += skipover


#HeterozyRec[chrom] = []

#HeterozyRec[chrom].append(0)
if re.search(r'^#', line) != None:
    linelist = vcffile.readlines()
else:
    exit(-1)



nextLine = -1#currentLineNum == nextStartLineNum 表示还没有找到下一个其实位置
#nextStartPosOfFirstSnp=int(re.split(r"\s+",linelist[0])[1])
isDetecting = True
currentLine = 0
totalRecs = len(linelist)

collist = re.split(r"\s+", linelist[currentLine])
chrom = collist[0]

#
HeterozyRec[chrom] = []
#
start_snp_pos = 0
firstComeInWin = True   
while currentLine != totalRecs:
    line = linelist[currentLine]
    collist = re.split(r"\s+", line)
    
    info = collist[7]
    pos = int(collist[1])
    if re.search(r"INDEL", info) == None:
        #last_snp_pos = pos
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", info)
        #print(dp4.group(1),dp4.group(2),dp4.group(3),dp4.group(4),collist)#,file=testfile
        if chrom == collist[0]:
            if pos - winStartPos <= windowWidth:
                if firstComeInWin:
                    start_snp_pos = pos
                    firstComeInWin = False
                last_snp_pos = pos#for the this window
                refdep = int(dp4.group(1)) + int(dp4.group(2))
                altalleledep = int(dp4.group(3)) + int(dp4.group(4))
                if refdep < altalleledep:
                    CNMI += refdep
                    CNMA += altalleledep
                else:
                    CNMA += refdep
                    CNMI += altalleledep
#                print(pos,winStartPos)
#                if nextLine == currentLine:
#                    start_snp_pos = pos#for this window                
                if pos - winStartPos > skipover and isDetecting:

#                    print(pos,"-",winStartPos,"=",pos - winStartPos,"find find find find**************************************")
                    nextLine = currentLine
                    #start_snp_pos = pos#for the next window
                    isDetecting = False
            else:
                # caculate the HETEROZY of the last window,and init variables go to the next window
                HETEROZY = CNMA * CNMI * 2 / ((CNMA + CNMI) ** 2)
                CNMI = 0;CNMA = 0
                HeterozyRec[chrom].append((start_snp_pos, last_snp_pos, HETEROZY));windowNo += 1
                firstComeInWin = True
                winStartPos += skipover
                if nextLine != -1:
#                    HeterozyRec[chrom].append(0)
#                    print(pos,currentLine,chrom,len(HeterozyRec[chrom]))
                    currentLine = nextLine
                    nextLine = -1
                    isDetecting = True               
                    continue#见到continue说明 下次循环 还在指定位置开始
                else:#the next snp position go beyond more than one windows 
                    while pos - winStartPos > windowWidth:
                        HeterozyRec[chrom].append((0, 0, 0))
                        windowNo += 1
#                        print(pos,currentLine,chrom,len(HeterozyRec[chrom]),"============================================")
                        winStartPos += skipover
                    else:
                        isDetecting = True
                        continue#keep the currentLine stay at the same value #见到continue说明 下次循环 还在指定位置开始

        else:
            print(chrom, collist[0])
            #deal with the aftermost window in the last chrom
            try:
                HETEROZY = CNMA * CNMI * 2 / ((CNMA + CNMI) ** 2)
            except ZeroDivisionError:
                print(currentLine, pos, dp4, collist)
                break
            HeterozyRec[chrom].append((start_snp_pos, last_snp_pos, HETEROZY));windowNo += 1
            firstComeInWin = True
            start_snp_pos = 0;
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))
            if refdep < altalleledep:
                CNMI = refdep
                CNMA = altalleledep
            else:
                CNMA = refdep
                CNMI = altalleledep
#            print(pos,currentLine,chrom,len(HeterozyRec[chrom]))
            chrom = collist[0]
            HeterozyRec[chrom] = []
            winStartPos = 0
            nextLine = -1
            isDetecting = True
#            currentLine = 0
#    else:
#        if chrom != collist[0]:
#            HETEROZY = CNMA * CNMI * 2 / ((CNMA + CNMI) ** 2)
#            HeterozyRec[chrom].append(HETEROZY);windowNo+=1
#            nextLine = -1
#            winStartPos = 0
#            chrom = collist[0]
#            HeterozyRec[chrom]=[]
#            isDetecting = True
#        elif pos - winStartPos > windowWidth:
#            HETEROZY = CNMA * CNMI * 2 / ((CNMA + CNMI) ** 2)
#            HeterozyRec[chrom].append(HETEROZY);windowNo+=1
#            if nextLine != -1:
##                    HeterozyRec[chrom].append(0)
##                    print(pos,currentLine,chrom,len(HeterozyRec[chrom]))
#                currentLine = nextLine
#                nextLine = -1
#                isDetecting = True               
#                continue#见到continue说明 下次循环 还在指定位置开始
#            else:#the next snp position go beyond more than one windows 
#                while pos - winStartPos>windowWidth:
#                    HeterozyRec[chrom].append(0)
#                    windowNo+=1
##                        print(pos,currentLine,chrom,len(HeterozyRec[chrom]),"============================================")
#                    winStartPos += skipover
#                else:
#                    isDetecting = True
#                    continue#keep the currentLine stay at the same value #见到continue说明 下次循环 还在指定位置开始
#            print("kkkkkkkkkkkkkkkkkkk")
#            break

    currentLine += 1
else:
    print("caculate last window")
    windowNo += 1



allwindow = []
for chrom in sorted(HeterozyRec.keys()):
    for i in range(len(HeterozyRec[chrom])):
        #print(chrom+"\t"+str(i)+"\t"+str(HeterozyRec[chrom][i]),file=fileout)
        if HeterozyRec[chrom][i][2] != 0:
            allwindow.append(HeterozyRec[chrom][i][2])
        else:
            pass

expectation = numpy.mean(allwindow)
std0 = numpy.std(allwindow)#ddof=0
std1 = numpy.std(allwindow, ddof=1)

print(pos, windowNo, len(allwindow), expectation, std0, std1)

del allwindow
ZHeterozyRec = {}

for chrom in HeterozyRec.keys():
    ZHeterozyRec[chrom] = []

"""structure of HetrozyRec{chrom:[(first_snp_pos,last_snp_pod,hp),(first_snp_pos,last_snp_pod,hp),.....],chrom:[],.....}
"""
for chrom in sorted(HeterozyRec.keys()):
    for i in range (len(HeterozyRec[chrom])):
        if HeterozyRec[chrom][i][2] != 0:
            ZHeterozyRec[chrom].append((HeterozyRec[chrom][i][2] - expectation) / std1)
            #print((HeterozyRec[chrom][i]-HeterozyRec[chrom][i]*expectation)/(std1*HeterozyRec[chrom][i]),std1*HeterozyRec[chrom][i],(HeterozyRec[chrom][i]-HeterozyRec[chrom][i]*expectation))
        else:
            ZHeterozyRec[chrom].append(0)
        print(chrom + "\t" + str(i) + "\t" + str(HeterozyRec[chrom][i][0]) + "\t" + str(HeterozyRec[chrom][i][1]) + "\t" + str(HeterozyRec[chrom][i][2]) + "\t" + str(ZHeterozyRec[chrom][i]), file=fileout)
#for chrom in sorted(HeterozyRec.keys()):
#    for i in range(len(HeterozyRec[chrom])):
#        print(chrom+"\t"+str(i)+"\t"+str(HeterozyRec[chrom][i])+"\t"+str(ZHeterozyRec[chrom][i]),file=fileout)


vcffile.close()
fileout.close()
#testfile.close()









