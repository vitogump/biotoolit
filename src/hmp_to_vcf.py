# -*- coding: UTF-8 -*-
import re,sys,csv,json
from optparse import OptionParser
'''
Created on 2017年9月20日

@author: liurui
'''
parser = OptionParser()
parser.add_option("-c", "--csvdata", dest="csvdata",help="used only when -d 0")
parser.add_option("-v", "--vcfconfig", dest="vcfconfig",help="used only when -d 0")# action="callback",type="string",callback=useoptionvalue_previous1,
parser.add_option("-d","--direction", dest="direction", help="if value is 0 convert -t type -i inputfile into vcf files. or value =1  covert -v assigned vcf files into -t type files")
"""
0: other format -> vcf
1: vcf -> other format
"""
parser.add_option("-t", "--inputtype", dest="inputtype")
parser.add_option("-i","--inputfiles",dest="inputfiles",action="append",help="vcf or hmp or map & ped")

parser.add_option("-o","--outputpre",dest="outputpre")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

infilehanders=[]
"""
    ind_idx,population_idx,mark1_idx,,,,,
"""
useforcsvfield=[3,6,2]#ind_idx,population_idx,mark1_idx


def read_csv(filename,usefultitlelist):
    
    csv_rowsJ=[]
    with open(filename) as csvfile:
        reader=csv.DictReader(csvfile)
        title=reader.fieldnames
        title[useforcsvfield[0]]="indvd"
        title[useforcsvfield[1]]="population"
        population_sets=set()
        popKeyindValue={}
        for row in reader:
            population_sets.add(row[title[useforcsvfield[1]]])# determine how many output vcf files
            csv_rowsJ.extend([{title[i]:row[title[i]] for i in usefultitlelist}])
        csvfile.close()
        print(population_sets,"\n",csv_rowsJ)
        for pop in population_sets:
            popKeyindValue[pop]=[]
            for row in csv_rowsJ:
                if row["population"]==pop:
                    popKeyindValue[pop].append(row["indvd"])
        return csv_rowsJ,popKeyindValue
def write_json(data, json_file,format=None):
    with open(json_file,"w") as f:
        f.write(json.dumps(data,sort_keys=False,indent=4,separators=(",",":"),ensure_ascii=False))
    f.close()
vcfheader="""##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##GATKCommandLine=<ID=UnifiedGenotyper,Version=3.3-0-g37228af,Date="Thu Aug 18 14:50:52 HKT 2016",Epoch=1471503052069,CommandLineOptions="analysis_type=UnifiedGenotyper input_fi
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
"""
if __name__ == '__main__':
    if options.direction.strip()=="0":
        print("inputfiles:",options.inputfiles)
        
        for e in options.inputfiles:
            print(e)
            infilehanders.append(open(e,'r'))
        csv_rowsJ,popKeyindValue=read_csv(options.csvdata, useforcsvfield)
        
        #open out vcf file and make a merged popKeyindValue_temp
        outvcffilelist={}
        popKeyindValue_temp={"allpop":[]}
        for pop in popKeyindValue.keys():
            f=open(pop+".txt",'w')
            outvcffilelist[pop]=open(options.outputpre+pop+".vcf",'w')
            outvcffilelist[pop].write(vcfheader)
            popKeyindValue_temp["allpop"].extend(popKeyindValue[pop])
            for indvd in popKeyindValue[pop]:
                print(indvd,file=f)
            f.close()
        """
 
        """        
        popKeyindValue["allpop"]=popKeyindValue_temp["allpop"]
        outvcffilelist["allpop"]=open(options.outputpre+"allpop"+".vcf",'w')  
        outvcffilelist["allpop"].write(vcfheader)
        
#         write_json(csv_rowsJ,"testtempjson")
#         resp=json.loads(open("testtempjson").read())

#         print(resp[0]["population"])
        if options.inputtype=="hmp":
            titlelist=re.split(r"\s+",infilehanders[0].readline().strip())
            #only for test
            for indvd in titlelist[11:]:
                for pop in popKeyindValue.keys():
                    if indvd in popKeyindValue[pop]:
                        break
                else:
                    print(indvd, "warning : different name between inputfile and csvfile ")
            
            #print vcf title line
            for pop in popKeyindValue.keys():
                
                print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",end="\t",file=outvcffilelist[pop])
                for indvd in titlelist[11:]:
                    if indvd in popKeyindValue[pop]:
                        print(indvd,end="\t",file=outvcffilelist[pop])
                print("",file=outvcffilelist[pop])
            #print vcf recs
            for snp_rec in infilehanders[0]:
                rec_list=re.split(r"\t",snp_rec.strip())
                #filter
                if rec_list[5]!="Build_4.0":
#                     print(rec_list[5],rec_list)
                    continue
                
                for pop in popKeyindValue.keys():# assign to different vcf; here only one 
                    print_listForOnePop=[]
                    ref_allale=re.split(r"/",rec_list[1])[0];alt_allale=re.split(r"/",rec_list[1])[1]# col4 and col5 of outputvcf 
                    print_listForOnePop.extend([rec_list[2],rec_list[3],rec_list[0],ref_allale,alt_allale,"3000","."])#init SNP REC
                    AN=0;AC=0
                    #append indvd info
#                     print(len(titlelist),len(rec_list))
                    for indvd_idx in range(11,len(rec_list)):# ordered as title, travel all indvd
#                         print(,indvd_idx)
#                         print(rec_list)
                        if titlelist[indvd_idx] in popKeyindValue[pop]:
                            if ref_allale==rec_list[indvd_idx][0] and alt_allale==rec_list[indvd_idx][1]:#0/1
                                AN+=2;AC+=1
                                print_listForOnePop.append("0/1:5,5:2:8:98,0,172")#PL0,0,0 may cause some effect in some vcf processing software 
                            elif ref_allale==rec_list[indvd_idx][0] and ref_allale==rec_list[indvd_idx][1]:#0/0
                                AN+=2
                                print_listForOnePop.append("0/0:5,5:2:8:0,30,404")
                            elif alt_allale==rec_list[indvd_idx][0] and alt_allale==rec_list[indvd_idx][1]:#1/1
                                AN+=2;AC+=2
                                print_listForOnePop.append("1/1:5,5:2:8:297,24,0")
                            elif "-" in rec_list[indvd_idx]:
                                print_listForOnePop.append("./.")
                    #append INFO field 
                    try :
                        AF=AC/AN
                    except ZeroDivisionError:
                        continue
                    INFOfield="AC="+str(AC)+";AF="+str(AF)[0:5]+";AN="+str(AN)+""
                    print_listForOnePop.insert(7,INFOfield) 
                    print_listForOnePop.insert(8,"GT:AD:DP:GQ:PL")

                            
                    print("\t".join(print_listForOnePop),sep="\t",file=outvcffilelist[pop])

    elif options.direction.strip()=="1":
        pass
    for pop in outvcffilelist.keys():
        outvcffilelist[pop].close()
    infilehanders[0].close()