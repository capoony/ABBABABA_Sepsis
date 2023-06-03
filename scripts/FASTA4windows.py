import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import os

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--FASTA", dest="FA", help="Input file")
parser.add_option("--windowsize", dest="WS", help="WindowSize")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


FASTA = d(list)

for l in load_data(options.FA):
    if l.startswith(">"):
        ID = l.split()[0].rstrip()[1:]
        continue
    FASTA[ID].extend(list(l.rstrip()))

for l in load_data(options.IN):
    a = l.rstrip().split(".")
    LN = a[1]
    CH = a[0]
    ST = int(LN)-int(int(options.WS)/2)
    EN = int(LN)+int(int(options.WS)/2)
    print(">"+l.strip())
    print("".join(FASTA[CH][int(ST)-1:int(EN)]))
