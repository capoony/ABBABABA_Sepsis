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
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

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


HSPs = d(lambda: d(list))
for l in load_data(options.IN):
    query_id, subject_id, scname, subject_len, query_len, pct_identity, aln_length, n_of_mismatches, gap_openings, q_start, q_end, s_start, s_end, e_value, bit_score = l.rstrip().split()
    if int(q_start) > int(q_end):
        HSPs[query_id][int(q_end)] = [int(q_end), int(
            q_start), subject_id, e_value, bit_score, aln_length]
    else:
        HSPs[query_id][int(q_start)] = [int(q_start), int(
            q_end), subject_id, e_value, bit_score, aln_length]

# print(HSPs)


for contig, v in HSPs.items():
    LIST = d(list)
    RA = ""
    for Start, v1 in sorted(v.items()):
        q_start, q_end, subject_id, e_value, bit_score, aln_length = v1
        # print(v1)
        if RA == "":
            RA = [q_start, q_end]
            LIST[bit_score] = v1
            continue
        if (q_start > RA[0]-1000 and q_start < RA[1]+1000) or (q_end < RA[1]-1000 and q_end > RA[0]+1000):
            LIST[bit_score] = v1
            X, Y = RA
            RA = [min(q_start, X), max(q_end, Y)]
        else:
            print(contig+"\t"+"\t".join([str(x)
                  for x in LIST[max(LIST.keys())]]))
            LIST = d(list)
            RA = ""
            RA = [q_start, q_end]
            LIST[bit_score] = v1
        # print(RA)
    print(contig+"\t"+"\t".join([str(x) for x in LIST[max(LIST.keys())]]))
