'''
Created on 2014-1-1

@author: rui
'''

from optparse import OptionParser
import sys, re,copy,gzip


tempfilename1=""
tempfilename2=""
#the two function below still can't give the -3 and -4 options default value 
def useoptionvalue_previous1(option, opt_str, value, parser):
    tempfilename1=value
#     print(value)
def useoptionvalue_previous2(option, opt_str, value, parser):
    tempfilename2=value
#     print(value)
parser = OptionParser()
parser.add_option("-1", "--reforderfile1", dest="reforderfile1",nargs=2,# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="infile1 afterRetainTheSameLeft")
parser.add_option("-2", "--infile2", dest="infilename2",action="append",nargs=2,# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="infile2 afterRetainTheSameLeft. multple")
parser.add_option("-c","--samecolumns",action="append",dest="columnlist",nargs=2,help="1 2")#ordered by priority level
parser.add_option("-r", "--rmna",action="store_true",default=False,help="rn NA")
parser.add_option("-f", "--fillna",action="store_true",default=False,help="rn NA")
parser.add_option("-w", "--winwidth",dest="winwidth",default=None,help="winwidth")
parser.add_option("-s", "--slideSize",dest="slideSize",default=None,help="slideSize")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

infile1map={}
infile2maplist=[]
priority1th=int(options.columnlist[0][0])-1
priority2nd=int(options.columnlist[0][1])-1
if options.reforderfile1[1].endswith("gz"):
    outfileref=gzip.open(options.reforderfile1[1],'wb')
else:
    outfileref=open(options.reforderfile1[1],'w')
outfile2list=[]

def countRegionlength(selectedWinMap,winwidth,slideSize,winType="winvalue",morethan_lessthan="m",mergeNA=2):
    selectedRegion={}
    for chrom in selectedWinMap:
        selectedWinMap[chrom].sort(key=lambda listRec: int(listRec[1]))
        selectedRegion[chrom]=[]
        mergedRegion=[selectedWinMap[chrom][0]]
        i=1
        while i < len(selectedWinMap[chrom]):
#             print(chrom,selectedWinMap[chrom][i])
#             try:
            extremeValue="NA"
            if int(selectedWinMap[chrom][i-1][1])+1==int(selectedWinMap[chrom][i][1]) or int(selectedWinMap[chrom][i-1][1])*slideSize+winwidth>=int(selectedWinMap[chrom][i][1])*slideSize:#continues win
                mergedRegion.append(selectedWinMap[chrom][i])
            else:#not continues
                #process last region
                Region_start=int(mergedRegion[0][1])*slideSize
                Region_end=int(mergedRegion[-1][1])*slideSize+winwidth
                Nwin=len(mergedRegion)
                extremeValues=[]
                noofsnps=[]
                for e in mergedRegion:
                    if winType=="winvalue":
                        extremeValues.append(float(e[5]))
                    elif winType=="zvalue": 
                        extremeValues.append(float(e[6]))
                    try:
                        noofsnps.append(int(e[4])) 
                    except:
                        noofsnps.append(0)   
                        
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue=min(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=max(extremeValues)
                maxNoSNP=max(noofsnps)
                mixNoSNP=min(noofsnps)  
                selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue,mixNoSNP,maxNoSNP))
                #process this win
                mergedRegion=[selectedWinMap[chrom][i]]
            i+=1
#             except IndexError:
#                 print(i,len(selectedWinMap[chrom]),selectedWinMap[chrom])
#                 exit(-1)
        else:
            Region_start=int(mergedRegion[0][1])*slideSize
            Region_end=int(mergedRegion[-1][1])*slideSize+winwidth
            Nwin=len(mergedRegion)
            extremeValues=[]
            noofsnps=[]
            for e in mergedRegion:
                if winType=="winvalue" and e[5] !="NA":
                    extremeValues.append(float(e[5]))
                elif winType=="zvalue" and e[6] !="NA": 
                    extremeValues.append(float(e[6]))
                try:
                    noofsnps.append(int(e[4]))
                except:
                    noofsnps.append(0)
            if (morethan_lessthan == "m" or morethan_lessthan == "M" ) and len(extremeValues)!=0:
                extremeValue=min(extremeValues)
            elif (morethan_lessthan == "l" or morethan_lessthan == "L" ) and len(extremeValues)!=0:
                extremeValue=max(extremeValues) 
            maxNoSNP=max(noofsnps)
            mixNoSNP=min(noofsnps)                      
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue,mixNoSNP,maxNoSNP))
    if mergeNA!=False and int(mergeNA)>0:
        print("skip")
#         for chrom in selectedRegion:
#             selectedRegion[chrom].sort(key=lambda listRec: int(listRec[1]))
#             i=1
#             idxlist_to_pop=[]
#             while i <len(selectedRegion[chrom]):
#                 winNo_end=str(int(selectedRegion[chrom][i][1]/slideSize))
#                 winNo_start=str(int((selectedRegion[chrom][i-1][2]-winwidth)/slideSize))
# #                 print("select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and  winNo<"+winNo_end)
#                 wincount_to_determine=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
#                 wincount_to_add=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintabletextvalueallwin + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
#                 if len(wincount_to_determine)==0 and len(wincount_to_add)<= int(mergeNA):
#                     if morethan_lessthan == "m" or morethan_lessthan == "M":
#                         extremeValue=min(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
#                     elif morethan_lessthan == "l" or morethan_lessthan == "L":
#                         extremeValue=max(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
#                     maxNoSNP=max(selectedRegion[chrom][i][3],selectedRegion[chrom][i-1][3])
#                     mixNoSNP=min(selectedRegion[chrom][i][3],selectedRegion[chrom][i-1][3])
#                     selectedRegion[chrom][i]=(chrom,selectedRegion[chrom][i-1][1],selectedRegion[chrom][i][2],selectedRegion[chrom][i-1][3]+selectedRegion[chrom][i][3]+len(wincount_to_add),extremeValue,mixNoSNP,maxNoSNP)
#                     idxlist_to_pop.append(i-1)
#                 i+=1
#             else:
#                 idxlist_to_pop.reverse()
#                 for idx_to_pop in idxlist_to_pop:
#                     selectedRegion[chrom].pop(idx_to_pop)
    else:
        for chrom in selectedRegion:
            selectedRegion[chrom].sort(key=lambda listRec: int(listRec[1]))
    totallength=0
    i=0
    for chrom in selectedRegion.keys():
#         print(selectedRegion[chrom])
        for chrom,r_s,r_e,Nwin,extremeValue,mixNoSNP,maxNoSNP in selectedRegion[chrom]:
            i+=1
            totallength+=(r_e-r_s)
    print(i)
    return copy.deepcopy(selectedRegion),totallength
if __name__ == '__main__':
    print(options.columnlist)
    if options.reforderfile1[0].endswith("gz"):
        infileref=gzip.open(options.reforderfile1[0],'rt')
    else:
        infileref=open(options.reforderfile1[0],'r')
    infile2list=[]
    for infilename,outfinename in options.infilename2:
        infile2maplist.append({})
        if infilename.endswith("gz"):
            infile2list.append(gzip.open(infilename,"rt"))
        else:
            infile2list.append(open(infilename,"r"))
        outfile2list.append(open(outfinename,"w"))
#     infile2=open(options.infilename2,'r')
    maporder=[]
#     startbasessMapByChr={}
    for line in infileref:
        linelist=re.split(r"\t",line.strip())
        if options.rmna and len(linelist)>5:
            if linelist[5]=="na" or linelist[5]=="NA" or re.search(r"inf", linelist[5])!=None or linelist[6]=="na" or linelist[6]=="NA" or re.search(r"inf", linelist[5])!=None:
                print("skip",linelist)
                continue
        maporder.append((linelist[priority1th].strip(),linelist[priority2nd].strip()))
        try:
            infile1map[linelist[priority1th].strip()][linelist[priority2nd].strip()]=linelist
        except KeyError:
#             if linelist[2].isdigit():
#                 startbasessMapByChr[linelist[priority1th].strip()]=int(linelist[2])
            infile1map[linelist[priority1th].strip()]={linelist[priority2nd].strip():linelist}
    for idx in range(len(options.infilename2)):
        for line in infile2list[idx]:
            linelist=re.split(r"\t",line.strip())
            if options.rmna and len(linelist)>5:
                if linelist[5]=="na" or linelist[5]=="NA" or re.search(r"inf", linelist[5])!=None or linelist[6]=="na" or linelist[6]=="NA" or re.search(r"inf", linelist[5])!=None:
                    print("skip",linelist)
                    continue
            try:
                infile2maplist[idx][linelist[priority1th].strip()][linelist[priority2nd].strip()]=linelist
            except KeyError:
                infile2maplist[idx][linelist[priority1th].strip()]={linelist[priority2nd].strip():linelist}
    for k1,k2 in maporder:
        Retain=True
        for infile2map in infile2maplist:
            if k1 in infile2map and k2 in infile2map[k1]:
                pass
            elif options.fillna :
                if  k1=="chrNo" and (k2=="winNo" or k2=="firstsnppos"):
                    
                    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile2list[infile2maplist.index(infile2map)])
                    
                    if Retain:
                        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfileref)
                    Retain=(Retain and False)
                    continue
                    print("title")
                if priority2nd==2:
                    try:
                        infile2map[k1][k2]=[str(k1),infile1map[k1][k2][1],int(k2),infile1map[k1][k2][3],"0","NA","NA"]
                    except KeyError:
                        infile2map[k1]={k2:[str(k1),infile1map[k1][k2][1],int(k2),infile1map[k1][k2][3],"0","NA","NA"]}
                elif priority2nd==1:
                    try:
                        infile2map[k1][k2]=[str(k1),str(k2),int(k2)*int(options.slideSize.strip()),int(k2)*int(options.slideSize.strip())+int(options.winwidth.strip()),"0","NA","NA"]
                    except KeyError:
                        infile2map[k1]={k2:[str(k1),str(k2),int(k2)*int(options.slideSize.strip()),int(k2)*int(options.slideSize.strip())+int(options.winwidth.strip()),"0","NA","NA"]}
            else:
                Retain=(Retain and False)
                print("except file1",k1,k2,"not found in file2",options.infilename2[infile2maplist.index(infile2map)])
        if Retain:
            print(*infile1map[k1][k2],sep="\t",file=outfileref)
            for infile2map_idx in range(len(infile2maplist)):
                print(*infile2maplist[infile2map_idx][k1][k2],sep="\t",file=outfile2list[infile2map_idx])

    outfileref.close()
    infileref.close()
    for idx in range(len(infile2maplist)):
        infile2list[idx].close()
        outfile2list[idx].close()
    def printregionlength(filename,winwidth,slideSize):
        infileref=open(filename,"r")
        infileref.readline()
        selectedWinMap={}
        for line in infileref:
            win = re.split(r"\s+",line.strip())
            if win[0] in selectedWinMap:
                selectedWinMap[win[0]].append(win)
            else:
                selectedWinMap[win[0]]=[win]
        if priority2nd  ==1: #this 'if' block is a temp process in case the exception when priority2nd==2      
            regionmap,regionlength=countRegionlength(selectedWinMap,winwidth,slideSize)
            print(filename,regionlength)
        infileref.close()
    if options.winwidth ==None:
        exit()
    printregionlength(options.reforderfile1[0],int(options.winwidth),int(options.slideSize))
    printregionlength(options.reforderfile1[1],int(options.winwidth),int(options.slideSize))


    for infilename,outfinename in options.infilename2:
        printregionlength(infilename,int(options.winwidth),int(options.slideSize))
        printregionlength(outfinename,int(options.winwidth),int(options.slideSize))
    print("only retain the follow fields\n","chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue")
