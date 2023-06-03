import sys 
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup
from collections import defaultdict as d
import random
import gzip
#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="\npython %prog --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

The purpose of this script is to extract a subsample (e.g. --subsample 0.5) or a defined number of polymorphic SNPs (e.g. --subsample 500) in the Individuals defined with --pop from a consensus file produced with extract_consensus.py (--input). For this SNP-dataset the script will calculate all pairwise r^2 values (according to Hill and Robertson 1968). There will be two output files: The first is a tab delimited file containing a tabular representation of the distance matrix, i.e. the columns consist of Chromosome, Position of SNP1, position of SNP2, r^2, alleles of SNP1 and alleles of SNP2.
The second file is a visual representation of the distance matrix using the LDheatmap R package, which needs to be installed before. This figure will show the physical genomic positions of SNPs used r^2 values highlighted in colors and for chromosomes with cosmopolitan inversion the breakpoints and the name of the inversion. Due to memory constraints, this script can only process one chromosmome at a time, which needs to be defined with (--chromosome).
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="Input consensus file")
parser.add_option("--ind", dest="i", help="indivduals used for the analysis")
parser.add_option("--region", dest="e", help="chromosome used for the analysis",default="NA")
parser.add_option("--N-cutoff", dest="u", help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--min-allele", dest="min", help="minimum allele frequency, e.g. 0.1")
parser.add_option("--window", dest="size", help="window-size, e.g. 200000")
parser.add_option("--step", dest="step", help="stepsize, e.g. 10000")
parser.add_option("--rdist", dest="rdist", help="maximum distance to calculate R-squared")
parser.add_option("--measure", dest="measure", help="D or R")

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

def isolate(x,ind):
    ''' get nucleotides for certain individuals'''
    nuc=""
    for i in ind:
        nuc+=x[i]
    return nuc

def median(x):
    mid =int(len(x)/2)
    sort=sorted(x)
    if len(x)==0:
        return None
    if len(x)%2==0:
        lower=sort[mid-1]
        upper=sort[mid]
        return (float(lower)+float(upper))/2.0
    else:
        return sort[mid]
    
def AfHf(x,y):
    ''' calculate allele and haplotype frequencies for two biallelic loci x and y '''
    a=x.replace("N","")[0]
    b=y.replace("N","")[0]
    # remove N's and only retain individuals without ambiguous allels at both loci
    X,Y=["".join(v) for v in zip(*[w for w in zip(x,y) if not "N" in w])]
    ab=["".join(w) for w in zip(X,Y)]
    pa=X.count(a)/float(len(X))
    pb=Y.count(b)/float(len(Y))
    pab=ab.count(a+b)/float(len(ab))
    return pa,pb,pab

def Rsquared(pa,pb,pab):
    ''' Hill Robertson 1968'''
    D=pab-pa*pb
    if (pa*(1-pa)*pb*(1-pb))==0:
        return "NA"
    r2=(D**2)/(pa*(1-pa)*pb*(1-pb))
    return r2

def Dprime(pA,pB,pAB):
    ''' Lewontin'''
    D=pAB-pA*pB
    pa=1-pA
    pb=1-pB
    if D>=0:
        Dmax = min([pA*pb, pa*pB])
    else:
        Dmax = max([-pA*pB, -pa*pb])
    Dp=D/Dmax
    return Dp

def MeanConf(data, confidence=0.95):
    import numpy as np
    import scipy as sp
    import scipy.stats
    a = 1.0*np.array(data)
    n = len(a)
    m, se, m2 = np.mean(a), scipy.stats.sem(a), np.median(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m-h,m,m2, m+h

def decay(P,R,N):
    r.assign("P",robjects.vectors.IntVector(P))
    r.assign("R",robjects.vectors.FloatVector(R))
    r('distance<-P')
    r('LD.data<-R')
    r('n<-'+N)
    r('HW.st<-c(C=0.00001)')
    r('HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=1000,warnOnly=TRUE))')
    r('tt<-summary(HW.nonlinear)')
    return str(list(r('tt$coefficients[1]'))[0])

individual=map(int,options.i.split(","))
cutoff=float(options.u)
#idhash=dict(zip(individual,positions))

if options.e!="NA":
    region=map(int,options.e.split(":"))

fullsnplist=d(str)
count=1

### read data
for l in load_data(options.input):
    if count%10000000==0:
        print count,"positions processed"
    count+=1
    C,P,D=l.split()
    
    if options.e!="NA":
        if int(P)<region[0] or int(P)>region[1]:
            continue
    
    DI=isolate(D,individual)
    #print l,DI
    
    ## test if less N's than expected
    if DI.count("N")/float(len(DI))>cutoff:
        continue
    
    ## only use positions with more than 1 allele
    if len(set(DI.replace("N","")))==1:
        continue
    
    a1=DI.replace("N","")[0]
    freq=DI.replace("N","").count(a1)/float(len(DI.replace("N","")))
    if freq<float(options.min) or freq>1-float(options.min):
        continue

    fullsnplist[int(P)]=DI

### iterate over data:
window=int(options.size)
step=int(options.step)

BIN=window
end=max(fullsnplist.keys())

while(BIN<=end):
    snpdata=range(BIN-window,BIN)
    midpoint=median(snpdata)
    SNPs = list(set(fullsnplist.keys()).intersection(list(snpdata)))
    BIN+=step
    
    if SNPs==[]:
        continue
    rlist=[]
    plist=[]
    for i in SNPs:
        for j in SNPs:
            if j<=i or abs(j-i)>int(options.rdist):
                continue
            PA,PB,PAB=AfHf(fullsnplist[i],fullsnplist[j])
            if options.measure=="R":
                RSQ=Rsquared(PA,PB,PAB)
            else:
                RSQ=Dprime(PA,PB,PAB)
            #print i,j,PA,PB,PAB,RSQ ## For test only!!
            if RSQ!="NA":
                rlist.append(RSQ)
                plist.append(j-i)
                N=len(fullsnplist[i].replace("N",""))
    if rlist==[]:
        continue
    #print midpoint, SNPs, rlist
    print str(midpoint)+"\t"+"\t".join(map(str,MeanConf(rlist)))+"\t"+decay(plist,rlist,str(N))
        
    
