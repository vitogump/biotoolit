'''
Created on 2019年8月19日

@author: liurui
'''
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--fafile", dest="fafile", action="append", help="only one > in each file")
parser.add_option("-l", "--length", dest="length", default=None, help="if None, use the first fafile length of bases per line")

parser.add_option("-o", "--outputprefix", dest="outputpre")

(options, args) = parser.parse_args()

fo = open(options.outputpre, "w")
TOTALLEN=0
if __name__ == '__main__':
    falist = []
    for fn in options.fafile:
        #read header >
        falist.append(open(fn, 'r'))
#     print(falist[0].readline(), file=fo)
    totalseq=">"
    for fa in falist:
        totalseq+=fa.readline().strip()[1:]
    print(totalseq,file=fo)
    #print header done
    totalseq = falist[0].readline().strip()
    lineWidth = len(totalseq)
    rearstr=0
    for fa in falist:
        k = 0
        #complement 
        oline=totalseq[k:k+lineWidth-rearstr]
        TOTALLEN+=len(oline.strip())
        print(oline.strip(),file=fo);k+=(lineWidth-rearstr)
        #new line
        oline=totalseq[k:k+lineWidth-rearstr]
        while len(oline) == lineWidth:
            TOTALLEN+=len(oline.strip())
            print(oline, file=fo);k+=lineWidth
            oline=totalseq[k:k+lineWidth]
        else:
            rearstr=len(oline.strip())
            TOTALLEN+=len(oline.strip())
            print(oline.strip(),end="",file=fo)
        #read 
        totalseq=""
        for inline in fa:
            totalseq+=inline.strip()
        fa.close()
    else:
        k=0
        oline=totalseq[k:k+lineWidth-rearstr]
        TOTALLEN+=len(oline.strip())
        print(oline.strip(),file=fo);k+=(lineWidth-rearstr)
        print("new line")
        oline=totalseq[k:k+lineWidth]
        while len(oline)==lineWidth:
            TOTALLEN+=len(oline.strip())
            print(oline,file=fo);k+=lineWidth
            oline=totalseq[k:k+lineWidth]
        else:
            TOTALLEN+=len(oline.strip())
            print(oline.strip(),file=fo)
    print(TOTALLEN)
    fo.close()
