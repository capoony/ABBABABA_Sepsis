import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import os

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

for file in os.listdir(options.IN):
    if not file.endswith(".cand"):
        continue
    ID = file.split(".cand")[0]
    printl=d(list)
    for l in load_data(options.IN+"/"+ID+".genes"):
        a=l.rstrip().split()
        printl[a[0]].append(a[1])
    for k,v, in sorted(printl.items()):
        print(ID,"\t".join(k.split(".")),",".join(v),sep="\t")
    