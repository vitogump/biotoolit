'''
Created on 2015-10-13

@author: liurui
'''
from optparse import OptionParser
import sys, re,copy

parser = OptionParser()
parser.add_option("-1", "--reforderfile1", dest="reforderfile1",nargs=2,# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="infile1 afterRetainTheSameLeft")
parser.add_option("-2", "--infile2", dest="infilename2",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="infile2 afterRetainTheSameLeft. multple")
parser.add_option("-c","--samecolumns",action="append",dest="columnlist",nargs=2,help="1 2")#ordered by priority level
parser.add_option("-r", "--rmna",action="store_true",default=False,help="rn NA")
parser.add_option("-f", "--fillna",action="store_true",default=False,help="rn NA")
parser.add_option("-w", "--winwidth",dest="winwidth",default=20000,help="winwidth")
parser.add_option("-S", "--submark",dest="submark",nargs=2,default=[8,"sexchromosome"],help="submark")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

subfile1map={}
maporder=[]
priority1th=int(options.columnlist[0][0])-1
priority2nd=int(options.columnlist[0][1])-1
if __name__ == '__main__':
    filetosub=open(options.infilename2,'r')
    infileref=open(options.reforderfile1[0],'r')
    outfileref=open(options.reforderfile1[1],'w')
#     outfileref=open(options.infilename2[1],'w')
    for line in filetosub:
        linelist=re.split(r"\s+",line.strip())
        if options.rmna and len(linelist)>5:
            if linelist[5]=="na" or linelist[5]=="NA" or re.search(r"inf", linelist[5])!=None or linelist[6]=="na" or linelist[6]=="NA" or re.search(r"inf", linelist[5])!=None:
                print("skip",linelist)
                continue
        maporder.append((linelist[priority1th].strip(),linelist[priority2nd].strip()))
        try:
            subfile1map[linelist[priority1th].strip()][linelist[priority2nd].strip()]=linelist
        except KeyError:
            subfile1map[linelist[priority1th].strip()]={linelist[priority2nd].strip():linelist}
    for line in infileref:
        linelist=re.split(r"\s+",line.strip())
        if linelist[priority1th].strip() in subfile1map and linelist[priority2nd].strip() in subfile1map[linelist[priority1th].strip()]:
            print(*subfile1map[linelist[0]][linelist[1]],sep="\t",file=outfileref)
        else:
            print(*linelist,sep="\t",file=outfileref)
#         if linelist[int(options.submark[0])-1]==options.submark[1]:
#             print(*subfile1map[linelist[0]][linelist[1]],sep="\t",file=outfileref)
#         else:
#             print(*linelist[:7],sep="\t",file=outfileref)
    outfileref.close()
    infileref.close()
    filetosub.close()
    
    
    """
    f=open(sys.argv[1],"r")
    of=open(sys.argv[2],"w")
    currentmiRNA=None
    seqcountmap={}
    for line in f:
        linelist=re.split(r"\s+",line.strip())
        #terminate last group,and init new
        if currentmiRNA!=linelist[5]:
            maxseq=None;maxseq_count=0
            for seq in seqcountmap.keys():
                if seqcountmap[seq]>maxseq_count:
                    maxseq=seq
                    maxseq_count=seqcountmap[maxseq]
            if maxseq!=None:
                print(">"+currentmiRNA+" "+str(maxseq_count)+"\n"+maxseq,file=of)
            #init new
            currentmiRNA=linelist[5]
            seqcountmap={}
        #collection
        if linelist[3] in seqcountmap:
            seqcountmap[linelist[3]]+=1
        else:
            seqcountmap[linelist[3]]=1
    of.close()
    f.close()
 """       
        