'''
Created on 2016-1-27

@author: liurui
'''
from math import floor
from optparse import OptionParser
import re, os


parser = OptionParser()
parser.add_option("-i","--inputwinfile",dest="inputwinfile",help="40k20k file")
# parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
# parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-o","--outfilename",dest="outfilename",help="20k10k file")#

(options, args) = parser.parse_args()
if __name__ == '__main__':
    splitedwinNoseqlist=[]
    infile=open(options.inputwinfile,"r")
    outfile=open(options.outfilename,"w")
    sumtest_txt=open("sumtxt","w")
    infile.readline()
    curchr="none"
    testfile=open("test.txt",'w')
    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
    for line in infile:
        linelist=re.split(r"\t",line.strip())
        if linelist[0].strip()!=curchr:
            
            if splitedwinNoseqlist !=[]:
                subwinstartNo=0;subwinendNo=splitedwinNoseqlist[-1][0][-2]
                print(curchr,"0","1","20000","100",splitedwinNoseqlist[0][1],splitedwinNoseqlist[0][2],sep='\t',file=outfile)
                if splitedwinNoseqlist[0][1]=="NA" or splitedwinNoseqlist[0][2]=="NA":
                    splitedwinNoseqlist[0][1]=0
                    splitedwinNoseqlist[0][2]=0
                if splitedwinNoseqlist[1][1]=="NA" or splitedwinNoseqlist[1][2]=="NA":
                    splitedwinNoseqlist[1][1]=0
                    splitedwinNoseqlist[1][2]=0
                print(curchr,"1","10001","30000","100",(float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[0][2])+float(splitedwinNoseqlist[1][2]))/2,sep='\t',file=outfile)
                print(curchr,"2","20001","40000","100",(float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[0][2])+float(splitedwinNoseqlist[1][2]))/2,sep='\t',file=outfile)
#                 print(curchr,"3","30001","50000","100",(float(splitedwinNoseqlist[2][1])+float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[2][2])+float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/3,sep='\t',file=outfile)
                for subwinNo in range(3,subwinendNo+1):
                    print("\n",curchr,subwinNo,end="\n",file=sumtest_txt)
                    startidx= floor((subwinNo-3)/2)
                    count=0;sumwinvalue=0;sumzvalue=0
                    for tempwinNo in range(startidx,len(splitedwinNoseqlist)):
                        if subwinNo in  splitedwinNoseqlist[tempwinNo][0]:
                            if splitedwinNoseqlist[tempwinNo][1]=="NA" or splitedwinNoseqlist[tempwinNo][2]=="NA":
                                continue
                            count+=1
                            sumwinvalue+=float(splitedwinNoseqlist[tempwinNo][1])
                            sumzvalue+=float(splitedwinNoseqlist[tempwinNo][2])
                            print("+",sumwinvalue,sumzvalue,end="\t",file=sumtest_txt)
                        elif (len(splitedwinNoseqlist)>tempwinNo+2 and subwinNo>splitedwinNoseqlist[tempwinNo+2][0][-1]) or subwinNo>splitedwinNoseqlist[-1][0][-1]:

                            break
                    if count==0:
                        continue
                    print(curchr,str(subwinNo),subwinNo*10000+1,subwinNo*10000+20000,"100",sumwinvalue/count,sumzvalue/count,sep="\t",file=outfile)
            curchr=linelist[0].strip()        
            splitedwinNoseqlist=[[(0,1,2,3),linelist[5],linelist[6]]]
            
        else:
            winNo=int(linelist[1])
            splitinto_subwinstartNo=2*winNo-1
            print(curchr,[(splitinto_subwinstartNo,splitinto_subwinstartNo+1,splitinto_subwinstartNo+2,splitinto_subwinstartNo+3,splitinto_subwinstartNo+4),linelist[5],linelist[6]],file=testfile)
            splitedwinNoseqlist.append([(splitinto_subwinstartNo,splitinto_subwinstartNo+1,splitinto_subwinstartNo+2,splitinto_subwinstartNo+3,splitinto_subwinstartNo+4),linelist[5],linelist[6]])
    else:
        subwinstartNo=0;subwinendNo=splitedwinNoseqlist[-1][0][-2]
        print(curchr,"0","1","20000","100",splitedwinNoseqlist[0][1],splitedwinNoseqlist[0][2],sep='\t',file=outfile)
        if splitedwinNoseqlist[0][1]=="NA" or splitedwinNoseqlist[0][2]=="NA":
            splitedwinNoseqlist[0][1]=0
            splitedwinNoseqlist[0][2]=0
        if splitedwinNoseqlist[1][1]=="NA" or splitedwinNoseqlist[1][2]=="NA":
            splitedwinNoseqlist[1][1]=0
            splitedwinNoseqlist[1][2]=0
#         print(curchr,"0","1","20000","100",float(splitedwinNoseqlist[0][1]),float(splitedwinNoseqlist[0][2]),sep='\t',file=outfile)
        print(curchr,"1","10001","30000","100",(float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[0][2])+float(splitedwinNoseqlist[1][2]))/2,sep='\t',file=outfile)
        print(curchr,"2","20001","40000","100",(float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[0][2])+float(splitedwinNoseqlist[1][2]))/2,sep='\t',file=outfile)
#         print(curchr,"3","30001","50000","100",(float(splitedwinNoseqlist[2][1])+float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/2,(float(splitedwinNoseqlist[2][2])+float(splitedwinNoseqlist[0][1])+float(splitedwinNoseqlist[1][1]))/3,sep='\t',file=outfile)
        for subwinNo in range(3,subwinendNo+1):
            startidx= floor((subwinNo-3)/2)
            count=0;sumwinvalue=0;sumzvalue=0
            for tempwinNo in range(startidx,len(splitedwinNoseqlist)):
                if subwinNo in  splitedwinNoseqlist[tempwinNo][0]:
                    if splitedwinNoseqlist[tempwinNo][1]=="NA" or splitedwinNoseqlist[tempwinNo][2]=="NA":
                        continue
                    count+=1
                    sumwinvalue+=float(splitedwinNoseqlist[tempwinNo][1])
                    sumzvalue+=float(splitedwinNoseqlist[tempwinNo][2])
                elif subwinNo>splitedwinNoseqlist[tempwinNo][0][-1]:
                    break
            if count==0:
                continue
            print(curchr,str(subwinNo),subwinNo*10000+1,subwinNo*10000+20000,"100",sumwinvalue/count,sumzvalue/count,sep="\t",file=outfile)
    testfile.close()
    outfile.close()
    sumtest_txt.close()