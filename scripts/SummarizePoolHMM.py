import scipy.stats as stats
import sys
import os
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--cand", dest="CA", help="Input file")
parser.add_option("--rand", dest="RA", help="Input file")
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


WS = int(int(options.WS)/2)
print("SweepCand\tSweepRand\tNoSweepCand\tNoSweepRand")
SweepCand,NoSweepCand,SweepRand,NoSweepRand=0,0,0,0

for files in os.listdir(options.CA):
    if not files.endswith(".pred"):
        continue
    
    LEN = int(files.split("/")[-1].split(".")[1])
    Pos = d(str)
    for i in range(LEN-WS, LEN+WS):
        Pos[i]
    NAME = files.split("/")[-1].split(".pred")[0]
    Predlist= []

    for l in load_data(options.CA+"/"+files):
        a = l.rstrip().split()
        if int(a[0]) in Pos:
                Predlist.append(a[-1])
    if Predlist.count("3") > len(Predlist)/4:
        SweepCand += 1
    else:
        NoSweepCand += 1

for files in os.listdir(options.RA):
    if not files.endswith(".pred"):
        continue

    LEN = int(files.split("/")[-1].split(".")[1])
    Pos = d(str)
    for i in range(LEN-WS, LEN+WS):
        Pos[i]
    NAME = files.split("/")[-1].split(".pred")[0]
    RandPredlist = []

    for l in load_data(options.RA+"/"+files):
        a = l.rstrip().split()
        if int(a[0]) in Pos:
                RandPredlist.append(a[-1])
    #print(len(RandPredlist))
    if RandPredlist.count("3")>len(RandPredlist)/4:
        SweepRand +=1
    else:
        NoSweepRand += 1


# importing packages

# creating data
data = [[SweepCand, SweepRand], [NoSweepCand, NoSweepRand]]
print(data)
# performing fishers exact test on the data
odd_ratio, p_value = stats.fisher_exact(data)
print('odd ratio is : ' + str(odd_ratio))
print('p_value is : ' + str(p_value))
