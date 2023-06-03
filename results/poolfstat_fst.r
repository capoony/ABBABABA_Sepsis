
#install.packages('poolfstat')
library('poolfstat')
library(tidyverse)
library(SuperExactTest)

##### Convert sync file to poolfstat file, and call SNPS

# We first have to give haploid sizes of each pool. Here, I had mostly 40 individuals per pool, and since I am working with diploid species, we multiply that by 2,  to get 80 individuals for most pools.
psizes <- as.numeric(c(50,50,50,50,50,50,50,50,50,50))
#psizes <- as.numeric(c(50,50,50,10))
# Then we give the names of each pool/sample.
pnames <- as.character(c('ZurichCyn','CevennesCyn','EstoniaCyn','PetroiaCyn','SorenbergCyn','ZurichNeo','CevennesNeo','GeschinenNeo','HospentalNeo','SorenbergNeo'))
#pnames <- as.character(c('EstoniaCyn','SorenbergCyn','SorenbergNeo','Sor'))
#pnames <- as.character(c('EstoniaCyn','CevennesCyn','CevennesNeo','Sor'))


SG.pooldata <- popsync2pooldata(sync.file = '/media/inter/mkapun/projects/ABBABABA_Sepsis/data/ABBA_BABA-filtered_4FST.sync.gz', poolsizes = psizes, poolnames = pnames,
                                     min.rc = 4, min.cov.per.pool = 10, max.cov.per.pool = 400,
                                     min.maf = 0.01, noindel = TRUE, nlines.per.readblock = 1e+06)

##### From this file we can compute global and per SNP FSTs
SG.pair.fst <- compute.pairwiseFST(SG.pooldata, method = 'Anova',output.snp.values = T)
SNPs=as.data.frame(SG.pair.fst@PairwiseSnpFST)

Table=data.frame('Chromosome'=SG.pooldata@snp.info,
    'Position'=SG.pooldata@snp.info,
    'SoC:SoN'=SNPs[['SorenbergCyn;SorenbergNeo']],
    'SoC:GeN'=SNPs[['SorenbergCyn;GeschinenNeo']],
    'SoC:PtC'=SNPs[['PetroiaCyn;SorenbergCyn']],
    'SoN:GeN'=SNPs[['GeschinenNeo;SorenbergNeo']],
    'SoN:PtC'=SNPs[['PetroiaCyn;SorenbergNeo']],
    'ZuC:ZuN'=SNPs[['ZurichCyn;ZurichNeo']],
    'ZuC:GeN'=SNPs[['ZurichCyn;GeschinenNeo']],
    'ZuC:PtC'=SNPs[['ZurichCyn;PetroiaCyn']],
    'ZuN:GeN'=SNPs[['ZurichNeo;GeschinenNeo']],
    'ZuN:PtC'=SNPs[['PetroiaCyn;ZurichNeo']])

DATA<-Table %>% 
    ## make windows of 1kb size and avoid scientific numeration
    mutate(bin = cut(Position, 
        breaks=seq(0,10000000,1000),
        dig.lab=10)) %>%
    ## convert bin factors to characters
    extract(bin, c('start', 'end'), '(-?\\d+),(-?\\d+)') %>%
    ## convert start end to numbers
    mutate_at(c('start','end'),as.numeric) %>%
    ## get average windows position
    mutate(window=rowMeans(across(c('start','end')))) %>%
    ## remove start/end
    select(-start,-end) %>%
    ## calculate average FST per Chrom and window
    group_by(Chromosome,window) %>%
    summarise(across(2:last_col(), funs(mean(., na.rm = TRUE))),
        Count=n())%>%
    ## make new window column
    unite(Window,Chromosome,window,sep='.') %>%
    ## only retain windows with more than 10 SNPs
    filter(Count>10) %>%
    select(-Count) %>%
    ## convert from wide to Long format
    gather(Comp,FST,2:last_col())

DATA.windows<-(DATA%>%
    spread(Comp,FST))

    write.table(DATA.windows,
        paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt'),
        row.names=F,
        col.names=F,
        quote = F)

    ## calculate length of unique windows
N=length(unique(DATA))

## get top 10% FST for each comparison 
DATA.top10 <- DATA %>%
     group_by(Comp) %>%
     arrange(Comp, desc(FST)) %>% 
     filter(FST > quantile(FST, .99,na.rm=TRUE))

DATA.top10.wide <- DATA.top10 %>%
    spread(Comp,FST)

## Intraspecific candidates
Intra=c('SoC.PtC_mean','SoN.GeN_mean','ZuC.PtC_mean','ZuN.GeN_mean')

for (i in seq(1,length(Intra),1)){
    DATA.intra <- DATA.top10%>%
        pivot_wider(names_from=Comp,values_from=FST) %>%
        select(Window,Intra[i])

    dir.create(paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',Intra[i]))
    write.table(na.omit(DATA.intra),
        paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',Intra[i],'/',Intra[i],'_Candidates_FST.txt'),
        row.names=F,
        col.names=F,
        quote = F)
}


## Interspecific candidates
Inter1=c('SoC.SoN_mean','SoC.SoN_mean','ZuC.ZuN_mean','ZuC.ZuN_mean')
Inter2=c('SoC.GeN_mean','SoN.GeN_mean','ZuC.PtC_mean','ZuN.GeN_mean')
InterL=c('SoC_Inter','SoN_Inter','ZuC_Inter','ZuN_Inter')

for (i in seq(1,length(Inter1),1)){
    DATA.inter <- na.omit(DATA %>%
        spread(Comp,FST) %>%
        mutate(DeltaFST=.data[[Inter1[i]]] - .data[[Inter2[i]]]) %>%
        select(Window,DeltaFST) %>%
        arrange(Window, desc(DeltaFST))) %>% 
        filter(DeltaFST > quantile(DeltaFST, .99))

    dir.create(paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',InterL[i]))
    write.table(data.frame(DATA.inter,DATA.inter),
        paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',InterL[i],'/',InterL[i],'_Candidates_FST.txt'),,
        row.names=F,
        col.names=F,
        quote = F)
}


