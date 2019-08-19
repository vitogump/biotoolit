'''
Created on 2016-10-25

@author: liurui
'''
import re,sys
if __name__ == '__main__':
    f=open(sys.argv[1],'r')
    fw=open(sys.argv[2],'w')
    print(f.readline().strip(),file=fw)
    print(f.readline().strip(),file=fw)
    for line in f:
        
        linelist=re.split(r"\s+",line.strip())
        print(len(linelist))
        for e in linelist[:4]:
            print(e,end="\t",file=fw)
        for e in linelist[4:32]:
            print(e)
            ee=re.search(r"^([\d\.][/|][\d\.])",e).group(1)
            print(ee,end="\t",file=fw)
        print(file=fw)
    f.close()
    fw.close()