### poolfstat 

##  4     ,5,    6,    7,    8,    9,    10,   11,  12,   13,   14
## "ZuC","PhC","PtC","SoC","MoC","ZuN","MoN","GeN","HoN","SoN","Sor"


mkdir -p /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4

NAMES=(SoC_PtC SoC_PhC SoN_GeN SoN_HoN ZuC_PtC ZuC_PhC ZuN_GeN ZuN_HoN MoC_PtC MoC_PhC MoN_GeN MoN_HoN)
CompFull=("6,7,13,14" "5,7,13,14" "11,13,7,14" "12,13,7,14" "6,4,9,14" "6,5,9,14" "11,9,4,14" "12,9,4,14" "6,8,10,14" "5,8,10,14" "11,10,8,14" "12,10,8,14")
P1full=(PtC PhC GeN HoN PtC PhC GeN HoN PtC PhC GeN HoN PtC PhC GeN HoN)
P2full=(SoC SoC SoN SoN ZuC ZuC ZuN ZuN MoC MoC MoN MoN)
P3full=(SoN SoN SoC SoC ZuN ZuN ZuC ZuC MoN MoN MoC MoC)
P4full=(Sor Sor Sor Sor Sor Sor Sor Sor Sor Sor Sor Sor)


for index in ${!NAMES[@]}; do
     P1=${P1full[index]}
     P2=${P2full[index]}
     P3=${P3full[index]}
     P4=${P4full[index]}
     Comp="-f1-3,"${CompFull[index]}
     Name=${NAMES[index]}

     echo cut ${Comp}

    gunzip -c /media/inter/mkapun/projects/ABBABABA_Sepsis/data/ABBA_BABA-filtered_4poolFST.sync.gz \
        | cut ${Comp} \
        | gzip > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}.sync.gz

    echo """

    #install.packages('poolfstat')
    library('poolfstat')
    library(tidyverse)
    library(SuperExactTest)

    ##### Convert sync file to poolfstat file, and call SNPS

    # We first have to give haploid sizes of each pool. Here, I had mostly 40 individuals per pool, and since I am working with diploid species, we multiply that by 2,  to get 80 individuals for most pools.
    psizes <- as.numeric(c(50,50,50,10))
    # Then we give the names of each pool/sample.
    pnames <- as.character(c('$P1','$P2','$P3','$P4'))

    SG.pooldata <- popsync2pooldata(sync.file = '/media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}.sync.gz', 
        poolsizes = psizes, 
        poolnames = pnames,
        min.rc = 4, min.cov.per.pool = 10, 
        max.cov.per.pool = 400,
        min.maf = 0.01, 
        noindel = TRUE, 
        nlines.per.readblock = 1e+06)


    STAT<-compute.fstats(SG.pooldata,nsnp.per.bjack.block = 1000,computeDstat = TRUE)
    NEW<-head(STAT@Dstat.values,3)
    write.table(NEW,
    file='/media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}_f4.stat',
    quote=F)

    """ > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}.r

    Rscript /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}.r

done

NAMES=(SoC_PhC SoN_GeN SoN_HoN ZuC_PtC ZuC_PhC ZuN_GeN ZuN_HoN MoC_PtC MoC_PhC MoN_GeN MoN_HoN)

awk 'NR<3' /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/SoC_PtC_f4.stat \
    >> /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Dstat.txt

for index in ${!NAMES[@]}; do
     Name=${NAMES[index]}

     awk 'NR==2' /media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4/${Name}_f4.stat \
    >> /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Dstat.txt

done

echo """
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

Table=data.frame('Chromosome'=SG.pooldata@snp.info$Chromosome,
    'Position'=SG.pooldata@snp.info$Position,
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
    spread(Comp,FST))$Window

    write.table(DATA.windows,
        paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt'),
        row.names=F,
        col.names=F,
        quote = F)

    ## calculate length of unique windows
N=length(unique(DATA$Window))

## get top 10% FST for each comparison 
DATA.top10 <- DATA %>%
     group_by(Comp) %>%
     arrange(Comp, desc(FST)) %>% 
     filter(FST > quantile(FST, .99,na.rm=TRUE))

DATA <- DATA%>%
    group_by(Comp) %>%
     arrange(Comp, desc(FST)) %>% 
     mutate(Canidates = ifelse(FST > quantile(FST, .99,na.rm=TRUE),
                       "yes", "no")) %>%
    arrange(Comp, desc(Window))

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
Inter1=c('SoC.SoN_mean'),'SoC.SoN_mean','ZuC.ZuN_mean','ZuC.ZuN_mean')
Inter2=c('SoC.PtC_mean'),'SoN.GeN_mean','ZuC.PtC_mean','ZuN.GeN_mean')
InterL=c('SoC_Inter'),'SoN_Inter','ZuC_Inter','ZuN_Inter')

for (i in seq(1,length(Inter1),1)){
    DATA.inter <- na.omit(DATA %>%
        spread(Comp,FST) %>%
        mutate(DeltaFST=.data[[Inter1[i]]] - .data[[Inter2[i]]]) %>%
        select(Window,DeltaFST) %>%
        arrange(Window, desc(DeltaFST))) %>% 
        filter(DeltaFST > quantile(DeltaFST, .99))

    dir.create(paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',InterL[i]))
    write.table(data.frame(DATA.inter$Window,DATA.inter$DeltaFST),
        paste0('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/',InterL[i],'/',InterL[i],'_Candidates_FST.txt'),,
        row.names=F,
        col.names=F,
        quote = F)
}


### plot 

DATA.So.Inter.SoC <- DATA%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,SoC.SoN_mean,SoC.PtC_mean) %>%
     mutate(FST=SoC.SoN_mean -SoC.PtC_mean,
        Name="Sorenberg_Interspecific:SoC")%>%
     select(-SoC.SoN_mean,-SoC.PtC_mean)%>%
     arrange(desc(FST)) %>% 
     mutate(Canidates = ifelse(FST > quantile(FST, .99,na.rm=TRUE),
                       "yes", "no")) %>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.So.Inter.SoN <- DATA%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,SoC.SoN_mean,SoN.GeN_mean) %>%
     mutate(FST=SoC.SoN_mean -SoN.GeN_mean,
        Name="Sorenberg_Interspecific:SoN")%>%
     select(-SoC.SoN_mean,-SoN.GeN_mean)%>%
        mutate(Canidates = ifelse(FST > quantile(FST, .99,na.rm=TRUE),
                    "yes", "no")) %>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.SoN.Intra <- DATA%>%
    filter(Comp == "SoN.GeN_mean") %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     mutate(FST=SoN.GeN_mean,
        Name="Sorenberg_Intraspecific:SoN")%>%
     select(-SoN.GeN_mean)%>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.SoC.Intra <- DATA%>%
    filter(Comp == "SoC.PtC_mean") %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,SoC.PtC_mean) %>%
     mutate(FST=SoC.PtC_mean,
        Name="Sorenberg_Intraspecific:SoC")%>%
     select(-SoC.PtC_mean)%>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.plot=rbind(DATA.So.Inter.SoC,
    DATA.SoC.Intra,
    DATA.So.Inter.SoN,
    DATA.SoN.Intra
    )

data_cum <- DATA.plot %>% 
  group_by(Chrom) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chrom, bp_add)

DATA.plot <- DATA.plot %>% 
  inner_join(data_cum, by = "Chrom") %>% 
  mutate(bp_cum = Pos + bp_add)

axis_set <- DATA.plot %>% 
  group_by(Chrom) %>% 
  summarize(center = mean(bp_cum))

 PLOT <- ggplot(DATA.plot, aes(x = bp_cum, y = FST, 
    color = as_factor(Chrom), 
    size = FST)) +
    ggtitle("Sörenberg")+
  facet_grid(Name~.)+
  geom_vline(data = . %>% filter(Canidates == "yes"),
    aes(xintercept = bp_cum),
    linewidth=0.2)+
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$Chrom, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$Chrom)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "FST") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  theme(strip.text.y.right = element_text(angle = 0))

PLOT

ggsave("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Sorenberg_FST.png",
     PLOT,
     width=16,
     height=6)

ggsave("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Sorenberg_FST.pdf",
     PLOT,
     width=16,
     height=6)

### plot Zurich

DATA.Zu.Inter.ZuC <- DATA%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,ZuC.ZuN_mean,ZuC.PtC_mean) %>%
     mutate(FST=ZuC.ZuN_mean -ZuC.PtC_mean,
        Name="Zurich_Interspecific:ZuC")%>%
     select(-ZuC.ZuN_mean,-ZuC.PtC_mean)%>%
     arrange(desc(FST)) %>% 
     mutate(Canidates = ifelse(FST > quantile(FST, .99,na.rm=TRUE),
                       "yes", "no")) %>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.Zu.Inter.ZuN <- DATA%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,ZuC.ZuN_mean,ZuN.GeN_mean) %>%
     mutate(FST=ZuC.ZuN_mean -ZuN.GeN_mean,
        Name="Zurich_Interspecific:ZuN")%>%
     select(-ZuC.ZuN_mean,-ZuN.GeN_mean)%>%
        mutate(Canidates = ifelse(FST > quantile(FST, .99,na.rm=TRUE),
                    "yes", "no")) %>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.ZuN.Intra <- DATA%>%
    filter(Comp == "ZuN.GeN_mean") %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     mutate(FST=ZuN.GeN_mean,
        Name="Zurich_Intraspecific:ZuN")%>%
     select(-ZuN.GeN_mean)%>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.ZuC.Intra <- DATA%>%
    filter(Comp == "ZuC.PtC_mean") %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,Canidates,ZuC.PtC_mean) %>%
     mutate(FST=ZuC.PtC_mean,
        Name="Zurich_Intraspecific:ZuC")%>%
     select(-ZuC.PtC_mean)%>%
     separate(Window,
        into=c("Chrom","Pos"),
        sep="\\.") %>%
    mutate_at(c('Pos'),as.numeric)%>%
    arrange(Chrom,Pos)

DATA.plot=rbind(DATA.Zu.Inter.ZuC,
    DATA.ZuC.Intra,
    DATA.Zu.Inter.ZuN,
    DATA.ZuN.Intra
    )

data_cum <- DATA.plot %>% 
  group_by(Chrom) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chrom, bp_add)

DATA.plot <- DATA.plot %>% 
  inner_join(data_cum, by = "Chrom") %>% 
  mutate(bp_cum = Pos + bp_add)

axis_set <- DATA.plot %>% 
  group_by(Chrom) %>% 
  summarize(center = mean(bp_cum))

 PLOT <- ggplot(DATA.plot, aes(x = bp_cum, y = FST, 
    color = as_factor(Chrom), 
    size = FST)) +
    ggtitle("Zürich")+
  facet_grid(Name~.)+
  geom_vline(data = . %>% filter(Canidates == "yes"),
    aes(xintercept = bp_cum),
    linewidth=0.2)+
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$Chrom, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$Chrom)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "FST") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  theme(strip.text.y.right = element_text(angle = 0))

PLOT

ggsave("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Zurich_FST.png",
     PLOT,
     width=16,
     height=6)

ggsave("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Zurich_FST.pdf",
     PLOT,
     width=16,
     height=6)

"""> /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolfstat_fst.r

Rscript /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolfstat_fst.r

### Do PoolHMM

###install Scipy==0.15.1

Tests=(SoC_Inter SoC_Inter) SoC.PtC_mean SoN_Inter SoN_Inter SoN.GeN_mean ZuC_Inter ZuC_Inter ZuC.PtC_mean ZuN_Inter ZuN_Inter ZuN.GeN_mean)
Samples=(SorenbergCyn SorenbergNeo) SorenbergCyn SorenbergNeo SorenbergCyn SorenbergNeo ZurichCyn ABBAZurichNeo ZurichCyn ABBAZurichNeo ZurichCyn ABBAZurichNeo)

for index in ${!Tests[@]}; do
     CO=${Tests[index]}
     CY=${Samples[index]}

    echo """

    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CY}
     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py \
          --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/${CY}.pileup.gz \
          --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CO}_Candidates_FST.txt \
          --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CY}/${CY} \
          --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CY}
     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CY}/${CY}*.pileup; do


          tmp=\${i##*/}
          ID=\${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix \${ID} \
               -n 50 \
               --only-spectrum \
               --theta 0.005 \
               -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix \${ID} \
               -n 50 \
               --pred \
               -k 0.000000000000001 \
               -s \${ID} \
               -e sanger
     done

    """ > /media/inter/mkapun/projects/ABBABABA_Sepsis/shell/${CO}_${CY}_poolhmm.sh

    sh /media/inter/mkapun/projects/ABBABABA_Sepsis/shell/${CO}_${CY}_poolhmm.sh &

done 

Tests=(SoC_Inter SoC_Inter SoC.PtC_mean SoN_Inter SoN_Inter SoN.GeN_mean ZuC_Inter ZuC_Inter ZuC.PtC_mean ZuN_Inter ZuN_Inter ZuN.GeN_mean)
Samples=(SorenbergCyn SorenbergNeo SorenbergCyn SorenbergNeo SorenbergCyn SorenbergNeo ZurichCyn ABBAZurichNeo ZurichCyn ABBAZurichNeo ZurichCyn ABBAZurichNeo)

for index in ${!Tests[@]}; do
     CO=${Tests[index]}
     CY=${Samples[index]}

    echo """

    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/RAND_${CY}

     # now random dataset

     WC=(\$(wc -l /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/${CO}_Candidates_FST.txt))
     N=\${WC[0]}

     sort -R /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt |
        head -\${N}  \
        >/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/Rand_FST.txt

     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py \
          --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/${CY}.pileup.gz \
          --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/Rand_FST.txt \
          --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/RAND_${CY}/RAND_${CY} \
          --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/RAND_${CY}

     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${CO}/RAND_${CY}/RAND_${CY}_*.pileup; do

          tmp=\${i##*/}
          ID=\${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix \${ID} \
               -n 50 \
               --only-spectrum \
               --theta 0.005 \
               -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix \${ID} \
               -n 50 \
               --pred \
               -k 0.000000000000001 \
               -s \${ID} \
               -e sanger
     done

     """ > /media/inter/mkapun/projects/ABBABABA_Sepsis/shell/${CO}_${CY}_rand_poolhmm.sh

    sh /media/inter/mkapun/projects/ABBABABA_Sepsis/shell/${CO}_${CY}_rand_poolhmm.sh &

done

cp /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/*/*_Candidates_FST.txt \
/media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap


## now test if there is significant overlap among windows 
mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap

First=(SoC_Inter SoC_Inter SoN_Inter SoN_Inter ZuC_Inter ZuC_Inter ZuN_Inter ZuN_Inter SoC_Inter SoC_Inter SoN_Inter SoN_Inter)
Second=(SoC.PtC_mean SoN.GeN_mean SoC.PtC_mean SoN.GeN_mean ZuC.PtC_mean ZuN.GeN_mean ZuC.PtC_mean ZuN.GeN_mean ZuC_Inter ZuN_Inter ZuC_Inter ZuN_Inter)
Names=(SoCISoC SoCISoN SoNISoC SoNISoN ZuCiZuC ZuCIZuN ZuNIZuC ZuNIZuN SoCIZuCI SoCIZuNI SoNIZuCI SoNIZuNI)


for index in ${!Names[@]}; do
     NA=${Names[index]}
     FI=${First[index]}
     SE=${Second[index]}

    echo """

    library(SuperExactTest)

    D1=read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${FI}/${FI}_Candidates_FST.txt')
    D2=read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/${SE}/${SE}_Candidates_FST.txt')
    N=nrow(read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt'))

    Comp<- list(
     ${FI}=D1[,1],
     ${SE}=D2[,1])

    res=supertest(Comp,n=N)

    write.table(data.frame('NAME'=summary(res)[1],
            'GENES'=summary(res)[8]),
            file='/media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/${NA}.txt',
            row.names=F,
            quote = F,
            sep='\t')

    """ > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap.r

    Rscript /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap.r

done

### extract candidates
for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/*.txt

do

    ID=${i%.*}
    grep '^11' ${ID}.txt \
        | awk -F "\t" '{print $8}' \
        | sed 's/, /\n/g' > ${ID}.cand

done


### is there an excess of sweeps ins candidate windows

mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergCyn \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/SorenbergCyn \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/SoC_Inter_SorenbergCyn.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergNeo \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/SorenbergNeo \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/SoC_Inter_SorenbergNeo.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC.PtC_mean/RAND_SorenbergCyn \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC.PtC_mean/SorenbergCyn \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/SoC.PtC_mean_SorenbergCyn.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoN.GeN_mean/RAND_SorenbergNeo \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoN.GeN_mean/SorenbergNeo\
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/SoN.GeN_mean_SorenbergNeo.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC_Inter/RAND_ZurichCyn \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC_Inter/ZurichCyn \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/ZuC_Inter_ZurichCyn.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC_Inter/RAND_ABBAZurichNeo \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC_Inter/ABBAZurichNeo \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/ZuC_Inter_ABBAZurichNeo.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC.PtC_mean/RAND_ZurichCyn \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuC.PtC_mean/ZurichCyn \
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/ZuC.PtC_mean_ZurichCyn.txt 

python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/SummarizePoolHMM.py \
    --rand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuN.GeN_mean/RAND_ABBAZurichNeo \
    --cand /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuN.GeN_mean/ABBAZurichNeo\
    --windowsize 1000 \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm_excess/ZuN.GeN_mean_ABBAZurichNeo.txt 


#### Test if significant overlap?

for i in /Volumes/MartinResearch2/Wolf2019/ABBABABA_Sepsis/results/overlap/*_FST.txt; do

     ID=${i%.*}

     # python3 /Volumes/MartinResearch2/Wolf2019/scripts/FASTA4windows.py \
     #      --input ${ID}.txt \
     #      --FASTA /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
     #      --windowsize 1000 \
     #      >${ID}.fasta

     sh /Volumes/MartinResearch2/Wolf2019/shell/doBlast_ABBA.sh \
          ${ID} &
done

for i in /Volumes/MartinResearch2/Wolf2019/ABBABABA_Sepsis/results/overlap/*.cand; do

     ID=${i%.*}

     # python3 /Volumes/MartinResearch2/Wolf2019/scripts/FASTA4windows.py \
     #      --input ${ID}.cand \
     #      --FASTA /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
     #      --windowsize 1000 \
     #      >${ID}.fasta

     sh /Volumes/MartinResearch2/Wolf2019/shell/doBlast_ABBA.sh \
          ${ID} &
done

## summarize for canidates:

python /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/tableCand.py \
    --input /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/CandWindGenesFST.txt

python /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/tableOverlap.py \
    --input /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/OverlapWindGenesFST.txt



### Summarize Overlap

awk 'NR==1' /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/SoCISoC.txt \
    > /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Overlap.txt

for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/*.txt

do

if [[ ${i} != *"FST"* ]];then
    awk 'NR==4' $i >> /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Overlap.txt
fi

done

### Summarize GO Terms

for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/*.go

do

tmp=${i##*/}
ID=${tmp%.*}

echo $ID

awk -v ID=$ID '{print ID"\t"$0}' $i >> /media/inter/mkapun/projects/ABBABABA_Sepsis/results/GO.txt

done

### Summarize Pi

echo """
library(tidyverse)
pnames <- as.character(c('Chrom','Pos','MoC','PhC','PtC','SoC','ZuC','ZuN','MoN','GeN','HoN','SoN'))
Spec <- as.character(c('S.neocynipsea','S.neocynipsea','S.cynipsea','S.neocynipsea','S.cynipsea','S.cynipsea','S.cynipsea','S.neocynipsea','S.cynipsea','S.neocynipsea'))
DATA=read.table("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/ABBABABA_new_10000_10000.pi",header=F)

colnames(DATA)<-pnames

DATA.means <- DATA %>%
    gather(Sample,Value,3:last_col()) %>%
    group_by(Sample)%>%
    summarise(Mean = mean(Value, na.rm=TRUE), SD = sd(Value, na.rm=TRUE), SE=SD/sqrt(n()))

DATA.means$Species <- Spec

DATA.new <- DATA.means%>%
    arrange(Species,Sample)

PLOT<- ggplot(DATA.new,aes(x=Sample,y=Mean,fill=Species))+
    geom_bar(stat="identity")+
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2,
                 position=position_dodge(.9))+
    scale_fill_manual(values=c("blue","red"))+
    ylab("pi")+
    xlab("")+
    theme_bw()

ggsave("/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Pi_plot_new.pdf",
    PLOT,
    width=5,
    height=3)

    ### Make Heatmap from FST

""" >

echo """


#install.packages('poolfstat')
library('poolfstat')
library(tidyverse)
library(SuperExactTest)

##### Convert sync file to poolfstat file, and call SNPS

psizes <- as.numeric(c(50,50,50,50,50,50,50,50,50,50,10))
pnames <- as.character(c('ZuC','MoC','PhC','PtC','SoC','ZuN','MoN','GeN','HoN','SoN','Sor'))

SG.pooldata <- popsync2pooldata(sync.file = '/media/inter/mkapun/projects/ABBABABA_Sepsis/data/ABBA_BABA-filtered_4poolFST_new.sync.gz', poolsizes = psizes, poolnames = pnames,
                                     min.rc = 4, min.cov.per.pool = 10, max.cov.per.pool = 400,
                                     min.maf = 0.01, noindel = TRUE, nlines.per.readblock = 1e+06)

##### From this file we can compute global and per SNP FSTs
FST <- compute.pairwiseFST(SG.pooldata, method = 'Identity')
pdf('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/FST_heatmap.pdf',
    width=6,
    height=6,)

heatmap(FST,col=colorRampPalette(c("yellow", "orange", "red"))(256))
legend(x="topright", 
    legend=c(min(FST@values[,1]), mean(FST@values[,1]), max(FST@values[,1])),fill=c("yellow","orange","red"))
dev.off()
