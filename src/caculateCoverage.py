'''
Created on 2014-12-25

@author: liurui
'''

from optparse import OptionParser
parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-c", "--cmdexample", dest="cmdtemplatefile",help="oneline scriptexamplefile")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")


parser.add_option("-s", "--scriptstorepath", dest="scriptstorepath", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",
                  help="1 :means produce cmdline scripts for every terminal folder,the input data should be all the data files under the terminal folder. 2:use all selected data files as the input parameters in the only one cmdline script")
parser.add_option("-I","--Interceptor_depth",dest="Interceptor_depth",default="0",help="depth of the folder to output")
parser.add_option("-l", "--interceptdirs", dest="interceptdirs",action="append", default=[], help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
if __name__ == '__main__':
    pass