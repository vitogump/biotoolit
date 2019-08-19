'''
Created on 2016-3-16

@author: liurui
'''
from optparse import OptionParser
import re, sys, os


def extractTrscptbyclass_code(gtfFileName,cuffcompare_out,pathtoCNCI,pathtoDNA2protein,pfamlib,pathtogffread,reffa,pfam_scan_pl):
    """cuffcompare_out="cuffcompare_out"
    gtfFileName="transcripts.gtf"
    """
    #extract trscript that contain more than 1 exon from high quality Trscpts
    a=os.system("""less -S """+gtfFileName+""" | cut -d ';' -f1,2 | sed 's/;/\t/g'| sed 's/gene_id//g'| sed 's/transcript_id//g'| sed 's/"//g'| grep -v 'transcript'| awk '{print$10}'|grep 'CUFF'| uniq -c|sed 's/   //g'| sed 's/^  //g'| sed 's/ /\t/g'| awk '{if($1>=2)print$2}' >transcripts_multi_exons_ID """)
    if a!=0:
        
        os.system(""" sed 's/.*gene_id \([^; ]*\);.*/\1/' """+gtfFileName+"""|sort|uniq -c|awk '{gsub(/"/,"");print $2}' > transcripts_multi_exons_ID""")
    #extract transcript longer than 200nt and FPKM greater than 0
    os.system(""" awk '{if($7!=0&&$11>200)print}' """+cuffcompare_out+""".transcripts.gtf.tmap  > """+cuffcompare_out+"""_g.transcripts.gtf.tmap_filtered_FPKM_len_200 """)
    #extract compared transcripts infomation
    os.system(""" grep -wFf transcripts_multi_exons_ID """+cuffcompare_out+"""_g.transcripts.gtf.tmap_filtered_FPKM_len_200 > """+cuffcompare_out+"""_g.transcripts.gtf.tmap_filtered_FPKM_multi_exons_len_200 """)
    #extract class_code
    os.system("""echo -e "j\ni\no\nu\nx"> class_code""")
    os.system(""" grep -wFf class_code """+cuffcompare_out+"_g.transcripts.gtf.tmap_filtered_FPKM_multi_exons_len_200 > candidate_lncRNA1")
    os.system(""" awk '{print $5}' candidate_lncRNA1 > candidate_lncRNA1_ID """)
    os.system(""" grep -wFf candidate_lncRNA1_ID """+gtfFileName +"> candidate_lncRNA1.gtf")
#     return "candidate_lncRNA1.gtf"
# """
    os.system(pathtogffread+" -w candidate_lncRNA.fa -g "+reffa+" candidate_lncRNA1.gtf")
# """
# def rm_n_Nseq(candidate_lncRNA_fa,pathtoCNCI,pathtoDNA2protein,pfamlib,gtfFileName,cuffcompare_out):
    #extract candidate lincRNA
    os.system(""" tr '\n' '\t' candidate_lncRNA.fa  |sed 's/>/\n>/g'| sed 's/\t/#/'| sed 's/\t//g'|sed 's/#/\t/g'| sed 's/gene=//g'|sed '/^$/d'| grep -v 'N'|grep -v 'n'|sed 's/\t/\n/g' >  candidate_lncRNA_2.fa""")
    os.system(""" tr '\n' '\t' candidate_lncRNA.fa  |sed 's/>/\n>/g'| sed 's/\t/#/'| sed 's/\t//g'|sed 's/#/\t/g'| sed 's/gene=//g'|sed '/^$/d'| grep -v 'N'|grep -v 'n' > candidate_lncRNA.hash""")
    curpwd=os.getcwd()
    os.chdir(pathtoCNCI)
    if pathtoCNCI[-1]=="/":
        pathtoCNCI=pathtoCNCI[:-1]

    CNCI_package_name=re.search(r"[^/]*$",pathtoCNCI)
    
    if CNCI_package_name.find("/")!=-1:
        
        os.system("""python ../"""+ CNCI_package_name+"""/CNCI.py -f """+curpwd+""""/candidate_lncRNA_2.fa -o """+curpwd+"""/CNCI_noncoding_out -m ve -p 30""")
    else:
        print("please give the path to updir of CNCI.py")
        exit(0)
    os.chdir(curpwd+"""/CNCI_noncoding_out""")
    os.system("""grep 'noncoding' CNCI.index | cut -f1| sed 's/_/ /g' > CNCI_noncoding_ID """)
    os.system("""grep -wFf CNCI_noncoding_ID ../candidate_lncRNA.hash  | sed 's/\t/\n/g' > ../CNCI_noncoding.fa""")
    os.chdir(curpwd)
    #translate into aa seq
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 1> CNCI_noncoding_AA1.fa""")
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 2 > CNCI_noncoding_AA2.fa""")
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 3 > CNCI_noncoding_AA3.fa""")
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 4 > CNCI_noncoding_AA4.fa""")
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 5 > CNCI_noncoding_AA5.fa""")
    os.system("""perl """+pathtoDNA2protein+""" -s CNCI_noncoding.fa -f 6 > CNCI_noncoding_AA6.fa""")
    #extract aa seq from the DNA aa seq
    os.system("""less -S CNCI_noncoding_AA1.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA1_1.fa """)
    os.system("""less -S CNCI_noncoding_AA2.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA2_1.fa """)
    os.system("""less -S CNCI_noncoding_AA3.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA3_1.fa """)
    os.system("""less -S CNCI_noncoding_AA4.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA4_1.fa """)
    os.system("""less -S CNCI_noncoding_AA5.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA5_1.fa """)
    os.system("""less -S CNCI_noncoding_AA6.fa | tr '\n' '\t'| sed 's/>/\n>/g'|sed '/^$/d'| awk '!(NR%2)'|tr '\t' '\n' | sed '/^$/d'  > CNCI_noncoding_AA6_1.fa """)
    #filter from pfam
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA1_pfam -fasta CNCI_noncoding_AA1_1.fa -dir """+pfamlib )
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA2_pfam -fasta CNCI_noncoding_AA2_1.fa -dir """+pfamlib)
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA3_pfam -fasta CNCI_noncoding_AA3_1.fa -dir  """+pfamlib)
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA4_pfam -fasta CNCI_noncoding_AA4_1.fa -dir """+pfamlib)
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA5_pfam -fasta CNCI_noncoding_AA5_1.fa -dir """+pfamlib)
    os.system(pfam_scan_pl+""" -outfile CNCI_noncoding_AA6_pfam -fasta CNCI_noncoding_AA6_1.fa -dir """+pfamlib)
    #rm trscpt matched with pfam
    os.system("""cat *pfam | grep -v '#'| sed '/^$/d'| cut -d ' ' -f1| cut -d '_' -f1| sort -u > pfam_hits_ID """)
    os.system("""grep -vwFf pfam_hits_ID ./CNCI_noncoding_out/CNCI_noncoding_ID | cut -d ' ' -f1 > removed_pfam_hits_ID """)
    os.system("""less CNCI_noncoding.fa | tr '\n' '\t'| sed 's/>/\n>/g'| sed '/^$/d' > CNCI_noncoding_tab """)
    #extract ID and seq
    os.system("""grep -wFf removed_pfam_hits_ID CNCI_noncoding_tab | tr '\t' '\n'| sed '/^$/d' > lncRNA.fa """)
    os.system("""grep -wFf removed_pfam_hits_ID CNCI_noncoding_tab | cut -f1|sed 's/>//g'| cut -d '.' -f1,2,3| sed 's/CUFF$//g' > lncRNA_ID """)
    #count length ,express value
    os.system("""less lncRNA.fa | tr '\n' '\t'| sed 's/>/\n>/g'| sed '/^$/d' > lncRNA_hash """)
    os.system("""perl length_stas.pl lncRNA_hash > lncRNA_length """)
    os.system("""grep -vwFf class_code"""+ cuffcompare_out+""".transcripts.gtf.tmap | awk '{if($3~/=/)print$13}' > known_transcripts_length""")
    #FPKM of lincRNA
    os.system("""awk 'OFS="\t"{print $5,$7}' """""".transcripts.gtf.tmap>"""+cuffcompare_out+""".transcripts.gtf.tmap_fpkm """  )
    os.system("""grep -wFf lncRNA_ID """+cuffcompare_out+""".transcripts.gtf.tmap_fpkm > lncRNA_fpkm""")
    #FPKM of known transcripts
    os.system("""grep -vwFf class_code  """+cuffcompare_out+""".transcripts.gtf.tmap_fpkm| awk 'OFS="\t"{if($3~/=/)print$1,$7}' > known_transcripts_FPKM """)
    #EXON number statistics
    os.system("""awk '{if($3!~/transcript/)print}' """+gtfFileName+"""| cut -d';' -f1,2| sed 's/;/\t/g'|awk 'OFS="\t"{print $1,$12}'| sed 's/"//g'| uniq -c| sed 's/^ \+//g' >exon_num_ID """)
    os.system("""grep -wFf lncRNA_ID exon_num_ID | awk 'OFS="\t"{print $3,$1}' > lncRNAID_exon_num""")
    os.system("""awk 'OFS="\t"{print $3,$1}' exon_num_ID | grep '^N' > known_transcript_exon_num""")
    
parser = OptionParser()
parser.add_option("-g", "--gtfFileName", dest="gtfFileName",
                  help="gtfFileName require transcript CDS [UTR], transcript_id  \"xxxxx\" ")
# parser.add_option("-t", "--transcripts_multi_exons_ID_filename", dest="transcripts_multi_exons_ID_filename", help="transcripts_multi_exons_ID_filename for output")
parser.add_option("-c", "--cuffcompare_out", dest="cuffcompare_out",default="cuffcompare_out",help="path_to_filename_string come from 'cuffcompare -o cuffcompare_out'")
parser.add_option("-r", "--reffa", dest="reffa",help="path to reffa ,used in gffread")
parser.add_option("-C", "--pathtoCNCI", dest="pathtoCNCI")
parser.add_option("-s","--pfam_scan_pl",dest="pfam_scan_pl",default="pfam_scan.pl")
parser.add_option("-f","--pathtogffread",dest="pathtogffread",default="gffread",help="pathtogffread,default value 'gffread'")
parser.add_option("-p","--pfamlib",dest="pfamlib",default=False)
# parser.add_option("-m","--mode",dest="mode",default=False)
parser.add_option("-2","--pathtoDNA2protein",dest="pathtoDNA2protein",help="path_to_DNA2protein.pl",default="DNA2protein.pl")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
if __name__ == '__main__':
    print(options.pathtoCNCI)
#     extractTrscptbyclass_code(options.gtfFileName,options.transcripts_multi_exons_ID_filename,options.cuffcompare_out)
#     os.system("/pub/tool/cufflinks-2.1.1/gffread -w "+candidate_lncRNA.fa+" -g "+Sus_scrofa.Sscrofa10.2.71.dna_sm.toplevel.fa candidate_lncRNA1.gtf")
    extractTrscptbyclass_code(options.gtfFileName,options.transcripts_multi_exons_ID_filename,options.cuffcompare_out,options.pathtoCNCI,options.pathtoDNA2protein,options.pfamlib,options.pathtogffread,options.reffa,options.pfam_scan_pl)
#     if options.mode.strip()=="1":
#         
#         """ /pub/tool/cufflinks-2.1.1/gffread -w candidate_lncRNA.fa -g Sus_scrofa.Sscrofa10.2.71.dna_sm.toplevel.fa candidate_lncRNA1.gtf """
#     elif options.mode.strip()=="2":
#         rm_n_Nseq(options.paramsfor2[0],options.paramsfor2[1],options.paramsfor2[2],options.paramsfor2[3],options.paramsfor2[4],options.paramsfor2[5])