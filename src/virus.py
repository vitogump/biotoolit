'''
Created on 2015-6-18

@author: liurui
'''
from optparse import OptionParser

from Bio import Phylo, SeqIO


parser = OptionParser()
parser.add_option("-f", "--fastaseqs", dest="fastaseqs",help="homologous file")

parser.add_option("-j","--conditiontojudge",dest="conditiontojudge",action="append",nargs=3,help="""posfromatg base1 base2
                                                                                                                red   blue""")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
inputfileName=options.fastaseqs
def encode_phyliplines(headers, sequences,maxlen=10):
    """
    This creates the contents of a non-interleaved phylip sequence file.
    @param headers: some header strings
    @param sequences: some sequence strings
    """
    nrows = len(headers)
    ncols = len(sequences[0])
    out_lines = ['%d %d' % (nrows, ncols)]
    for h, seq in zip(headers, sequences):
        out_h = h[:maxlen].ljust(maxlen)
        out_lines.append(out_h + seq)
    return '\n'.join(out_lines)
if __name__ == '__main__':
    muscleout_seqmap={}
    muscleout_seqgenerator=SeqIO.parse(inputfileName,"fasta")
    for seq_rec in muscleout_seqgenerator:
        seq="".join(seq_rec.seq)
        seqtojudge=seq.replace("-","")
        startpos=seqtojudge.find("ATG")

        if startpos==-1 or len(seqtojudge)<=startpos+int(options.conditiontojudge[0][0]):
            print(startpos,"len",seq_rec.id)
            muscleout_seqmap[seq_rec.id]=seq
        elif seqtojudge[startpos+int(options.conditiontojudge[0][0])].strip().upper==options.conditiontojudge[0][1].strip().upper():
#             print(seq_rec.id+"red")
            muscleout_seqmap[seq_rec.id+"red"]=seq
        elif seqtojudge[startpos+int(options.conditiontojudge[0][0])].strip().upper==options.conditiontojudge[0][2].strip().upper():
#             print(seq_rec.id+"blue")
            muscleout_seqmap[seq_rec.id+"blue"]=seq
        else:
            print(seq_rec.id)
            muscleout_seqmap[seq_rec.id]=seq
    #output
    pamlinputcdsheader=[];pamlinputcdsseq=[];maxlenlist=[]
    fastaalnoutfile=open(inputfileName+"marked",'w')
    for seqid in muscleout_seqmap.keys():
        pamlinputcdsheader.append(seqid)
        pamlinputcdsseq.append(muscleout_seqmap[seqid])
        maxlenlist.append(len(seqid))
        print(">"+seqid,file=fastaalnoutfile)
        k = 0
        cdsstrline = muscleout_seqmap[seqid][k:k + 60]
        while len(cdsstrline) == 60:
            print(cdsstrline, file=fastaalnoutfile);k += 60
            cdsstrline = muscleout_seqmap[seqid][k:k + 60]
        else:
            print(cdsstrline, file=fastaalnoutfile)
    fastaalnoutfile.close()
    maxlen=max(maxlenlist)
    try:
        print(encode_phyliplines(pamlinputcdsheader, pamlinputcdsseq,maxlen+2), file=open(inputfileName+"marked.phy",'w'))
    except:
        print("except")