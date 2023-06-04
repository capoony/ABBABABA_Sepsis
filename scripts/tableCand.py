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

DATAcand=d(lambda:d(list))
DATAfst = d(lambda: d(str))
for file in os.listdir(options.IN):
    if "Candidates" not in file:
        continue

    for l in load_data(options.IN+"/"+file):
        if file.endswith(".genes"):
            a=l.rstrip().split()
            DATAcand[file.split(".genes")[0]][a[0]].append(a[1])
        if file.endswith(".txt"):
            a = l.rstrip().split()
            DATAfst[file.split(".txt")[0]][a[0]] = a[1]
    
for k,v in sorted(DATAcand.items()):
    for w,v1 in sorted(v.items()):
        print(k,
              "\t".join(w.split(".")),
              DATAfst[k][w],
              ",".join(v1),
              sep="\t")


