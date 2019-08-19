# -*- coding: UTF-8 -*-
import re
from optparse import OptionParser
# first param is a .fasta file need to be split into serverl small .fasta file
# second param is a perfix name
# third param is a specify chromosome,if exist,this programma will extract this chromosome's fasta seq and write in to outfile.if not exist this programme will extract all fasta seq according '>'
parser = OptionParser()
parser.add_option("-r","--reffafilename",dest="reffafilename",help="depth of the folder to output")
parser.add_option("-c","--chrom_extract",dest="chrom_extract",default=None,help="require least chrom length")

parser.add_option("-o","--outputfileprefix",dest="outputfileprefix",default="")

(options, args) = parser.parse_args()
FaFileName = options.reffafilename
OutFilePre = options.outputfileprefix
if options.chrom_extract!=None:
    myChromID = options.chrom_extract
    OutFileName = OutFilePre+myChromID+".fa"
#else:
#    OutFileName =0

fafile = open(FaFileName,"r")
#findex = open(FileIndex,"w")


outfile=None

Found = False
for line in fafile:
    if re.search(r"^[>]",line)!=None:
        if outfile!=None:
            outfile.close()
        outfile=None
        collist = re.split(r"\s+",line)
        chromID = re.search(r"[^>]+",collist[0]).group(0)
        #print(chromID)
        if options.chrom_extract!=None:# this block is code for select a specify chromosome,delete this block,this programma still work for split all ">" in fa file
            if chromID==myChromID:
                Found=True
                print(line,chromID,myChromID)
                i=0##########
                try:
                    print(line,file = outfile,end="")
                except IndexError:
                    outfile=open(OutFileName,"w")
                    print(line,file = outfile,end="")
            else :
                if Found == True:
                    break
                else: 
                    for line in fafile:
                        if re.search(r"^[>]",line)!=None:
                            collist = re.split(r"\s+",line)
                            chromID = re.search(r"[^>]+",collist[0]).group(0)
                            if chromID == myChromID:
                                Found=True
                                outfile=open(OutFileName,"w")
                                print(line,file = outfile,end="")
                                break
        else:
            if outfile!=None:
                print(line,file = outfile,end="")
            else:
                outfile=open(OutFilePre+chromID+".fa","w")       
                print(line,file = outfile,end="")

    else:
#        print(i)
        print(line,file = outfile,end="")

fafile.close()    



