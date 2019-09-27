'''
Created on 2016-2-17

@author: liurui
'''
from optparse import OptionParser
import sys, re,gzip


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
parser.add_option("-2", "--infile2", dest="infilename2",nargs=2,# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="infile2 afterRetainTheSameLeft. multple")
parser.add_option("-c","--samecolumns",action="append",dest="columnlist",nargs=2,help="1 2")#ordered by priority level
parser.add_option("-r", "--rmna",action="store_true",default=False,help="rn NA")
parser.add_option("-f", "--fillna",action="store_true",default=False,help="rn NA")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

priority1th=int(options.columnlist[0][0])-1
priority2nd=int(options.columnlist[0][1])-1


if __name__ == '__main__':

    def printdiff(inreffilename,outreffile,infilename2nd):
        if inreffilename.endswith("gz"):
            infileref=gzip.open(inreffilename,'rt')
        else:
            infileref=open(inreffilename,'r')
#         infileref=open(inreffilename,'r')
        outfileref=open(outreffile,'w')
        if infilename2nd.endswith("gz"):
            infile2list=gzip.open(infilename2nd,"rt")
        else:
            infile2list=open(infilename2nd,"r")        
#         infile2list=open(infilename2nd,"r")
        infile1map={}
        infile2map={}
        maporder=[]
        for line in infileref:
            linelist=re.split(r"\s+",line.strip())
            maporder.append((linelist[priority1th].strip(),linelist[priority2nd].strip()))
            try:
                infile1map[linelist[priority1th].strip()][linelist[priority2nd].strip()]=linelist
            except KeyError:
                infile1map[linelist[priority1th].strip()]={linelist[priority2nd].strip():linelist}
        for line in infile2list:
            linelist=re.split(r"\s+",line.strip())
            try:
                infile2map[linelist[priority1th].strip()][linelist[priority2nd].strip()]=linelist
            except KeyError:
                infile2map[linelist[priority1th].strip()]={linelist[priority2nd].strip():linelist}
        for k1,k2 in maporder:
            Retain=True
            if k1 in infile2map and k2 in infile2map[k1]:
                pass
            else:
                Retain=(Retain and False)
            if not Retain:
                print(*infile1map[k1][k2],sep="\t",file=outfileref)
        infileref.close()
        outfileref.close()
        infile2list.close()
    printdiff(options.reforderfile1[0],options.reforderfile1[1],options.infilename2[0])
    printdiff(options.infilename2[0],options.infilename2[1],options.reforderfile1[0])