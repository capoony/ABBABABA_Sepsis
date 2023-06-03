import sys
import gzip
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--names", dest="NA", help="Input file")
parser.add_option("--pos", dest="PO", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")

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


OUT = d(str)
names = options.NA.split(",")
for i in names:
    OUT[i] = gzip.open(options.OUT+i+".pileup.gz", "wt")
ST = [int(x)*3 for x in options.PO.split(",")]
EN = [int(x)*3+2 for x in options.PO.split(",")]

C = 1
for l in load_data(options.IN):
    a = l.rstrip().split()
    if C % 1000000 == 0:
        print(C, "pos processed")
    for i in range(len(names)):
        PRINT = a[:3]
        PRINT.extend(a[ST[i]:EN[i]])
        # print(PRINT)
        OUT[names[i]].write("\t".join(PRINT))
    C += 1
for i in names:
    OUT[i].close()
