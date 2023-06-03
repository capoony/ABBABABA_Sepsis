import sys
import math
from optparse import OptionParser, OptionGroup
import collections
# from rpy2.robjects import r
# import rpy2.robjects as robjects
from bisect import bisect

#Author: Martin Kapun
#########################################################   HELP   #########################################################################
usage="python %prog --input input.fst --window-size 200000 --output out  --names Florida,Maine,IiR,IsR,YSR --ylim 0.5 --var --comp IiR,IsR,YSR"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Requirements:

R version 2.15

python module rpy2  

This script calculates average FST in sliding windows based on the output from fst-sliding.pl in PoPoolation2 and FST_weir_cockerham.py. By now, the script only non-overlapping windows, whose size can be determined with --window-size. Additionally, one has to define the names of all populations in the comparsion with the parameter --names. The names have to be in the same order as in the sync file, which was the input for the FST scripts (fst-sliding.pl in PoPoolation2 and FST_weir_cockerham.py). You can restrict the comparsions to only specific populations (which have to be a subset of --names) with the parameter --comp. The script generates two outputs: a text file containing all averaged FSTs and a lineplot for each chromosome separately. The parameter --ylim can be used to reduce the y-axis (e.g. to an maximum FST of 0.8). Additionally, by setting --var you can restrict the analysis to only variable positions and ignore invariable positions, i.e. where FST=0 in all comparisons. --out has to be set to specifiy the output file(s).
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="pi-Data")
parser.add_option("--window-size", dest="bin", help="binsize")
parser.add_option("--output", dest="o", help="output files")
parser.add_option("--names", dest="n", help="")
parser.add_option("--ylim", dest="ylim", help="ylimit in figures",default="1")
parser.add_option("--comp", dest="c", help="pops in the comparison")
parser.add_option("--stat", dest="s", help="Which stat to use: FST-nei:1, FST-wc:2, dxy:3")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
	''' import data either from a gzipped or or uncrompessed file or from STDIN'''
	import gzip         
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y


bin=int(options.bin)
names=options.n.split(",")
comp=options.c.split(",")
ylim=options.ylim
stat=int(options.s)-1

binhash=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))))

comphash={}
for i in comp:
    x,y=i.split(":")
    comphash[str(names.index(x)+1)+":"+str(names.index(y)+1)]=i
    comphash[str(names.index(y)+1)+":"+str(names.index(x)+1)]=i
#print comphash
con=1   
for l in load_data(options.input):
    a=l.split()
    chr,pos=a[:2]
    fstset=a[2:]
    if con%100000==0:
        print con,"processed"
        
    ## loop through all columns defined above and store the counts and sums of FSTs to later calculate the averages
    con+=1
    for pop in fstset:
        GO=pop.split("=")
        if len(GO)!=2:
            continue
        ID,val=GO 
        #print pop
        if ID not in comphash:
            continue
        binhash[chr][int(int(pos)/bin)][comphash[ID]]["count"]
        binhash[chr][int(int(pos)/bin)][comphash[ID]]["value"]
        if "NA" in val:
            continue
        ## use the divisions of positions by bins as index! Clever!! Thanx Robs
        binhash[chr][int(int(pos)/bin)][comphash[ID]]["count"]+=1
        binhash[chr][int(int(pos)/bin)][comphash[ID]]["value"]+=float(val.split(",")[stat])

# r('pdf("'+options.o+'.pdf",width=20,height=15)')
# r('par(mfrow=c(2,1),cex=1.5)')

out=open(options.o+"_binned.fst","w")

out.write("Chrom\tPos\tCount\tComp\tFST\n")
## now loop through the dictionary of dictionaries of dictionaries and print results
for chr,hash1 in sorted(binhash.items()):
    bindict=collections.defaultdict(list)
    binlist=[]
    
    for bins,idhash in sorted(hash1.items()):
        outlist,count=[],0
        binlist.append(bins*bin+float(bin)/2)
        for ID, values in idhash.items():
            count=values["count"]
            if values["count"]==0:
                out.write("\t".join(map(str,[chr,bins*bin+float(bin)/2,count,ID,"NA"]))+"\n")
            else:
                out.write("\t".join(map(str,[chr,bins*bin+float(bin)/2,count,ID,values["value"]/values["count"]]))+"\n")
            #outlist.append(ID+"="+str(values["value"]/values["count"]))
        
            #out.write("\t".join(map(str,[chr,bins*bin+float(bin)/2,count]+outlist))+"\n")
    #r('color<-c("black","red","blue","green")')
    # r.assign("position",robjects.vectors.FloatVector(binlist))
    # test=0
    # i=1
#     for key,val in sorted(bindict.items()):
#         r.assign("FST",robjects.vectors.FloatVector(val))
#         if test==0:
#             r('plot(position,FST,type="l",lwd=3,col=color[1],ylim=c(0,'+ylim+'),main=expression(italic("'+chr+'")))')
#             test=1
#         else:
#             r('points(position,FST,lwd=3,col=color['+str(i)+'],type="l")')
#         i+=1
#     r.assign("names",robjects.vectors.StrVector(sorted(bindict.keys())))
#     r('legend("topleft",names,col=color,lwd=3)')
# r('dev.off()')