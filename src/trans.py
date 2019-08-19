
import re
f=open("GCF_000298735.2_Oar_v4.0_genomic.gff",'r')
transcripid="";t=False
for eachline in f:
    if "#" in eachline:
        continue
    linelist=re.split(r"\t",eachline.strip())
    if eachline.strip().split()[2]=="gene":
        print(eachline.strip())
        continue
        #print(eachline)
    
    if eachline.strip().split()[2]=="mRNA":
        linelist[2]="transcript"
        tt=re.search(r'transcript_id=(.*)', eachline)
        if tt !=None:
            transcripid=tt.group(1).strip()
            print(*linelist[:8],'transcript_id "'+transcripid+'"; '+linelist[8],sep="\t")
        continue
    print(*linelist[:8],'transcript_id "'+transcripid+'"; '+linelist[8],sep="\t")   




    
