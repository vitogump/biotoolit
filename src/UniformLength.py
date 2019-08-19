import re,sys,os

if len(sys.argv)!=3:
    print("python UniformLength [fa/fqfile] [minium length]")
    exit(-1)
path  = sys.argv[1]
minlen = int(sys.argv[2])

# fileName=re.subn(r'.*?\\','',inputFileName)
def foreachfile(FileName,minlen):

    if re.search(r'[^.]*$', FileName).group(0).lower()=='fa' or re.search(r'[^.]*$', FileName).group(0).lower()=='fasta':
        infile = open(FileName,'r')
        outfile = open(FileName+".uniform"+sys.argv[2],'a')        
        seq = ''
        preline=''
        for line in infile:
            if re.search(r'^>',line)!= None:
                if len(seq)>=minlen:
                    seq = seq[0:minlen]
                    print(preline,end='',file=outfile)
                    print(seq,file=outfile)
                seq=''
                preline = line
            else:
                seq += line.strip()
                    
    elif re.search(r'[^.]*$', FileName).group(0).lower()=='fq' or re.search(r'[^.]*$', FileName).group(0).lower()=='fastq':
        pass
#print(re.subn(r'.*?\\','',inputfile))
    infile.close()
    outfile.close()

if os.path.isdir(path):
    for FileName in os.listdir(path):
        FileName = os.path.join(path,FileName)
        if not os.path.isdir(FileName):
            foreachfile(FileName,minlen)
    pass
else:
    foreachfile(path,minlen)

