import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--windows", dest="WIN", help="Output file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--windowsize", dest="WS", help="Output file")

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
Pos = d(lambda: d(str))
WS = int(10*int(options.WS)/2)
for l in load_data(options.WIN):
    X = l.rstrip().split()
    OUT[X[0]] = open(options.OUT+"_"+X[0]+".pileup", "wt")
    a = X[0].split(".")
    for i in range(int(a[1])-WS, int(a[1])+WS):
        Pos[a[0]][i] = X[0]

for l in load_data(options.IN):
    a = l.rstrip().split()
    if int(a[1]) in Pos[a[0]]:
        OUT[Pos[a[0]][int(a[1])]].write(l)
