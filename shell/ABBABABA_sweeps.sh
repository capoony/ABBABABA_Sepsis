# thoracica

mkdir -p /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho

#CevennesCyn,EstoniaCyn,PetroiaCyn,SorenbergCyn,ZurichCyn,ABBAZurichNeo,CevennesNeo,GeschinenNeo,HospentalNeo,SorenbergNeo

gunzip -c /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.sync.gz |
     cut -f 1-3,6,7,12,13,15,17,19,20,21,25 |
     grep -v ".:.:.:.:.:." |
     gzip >/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_4FST.sync.gz

python2 /Volumes/MartinResearch2/Wolf2019/scripts/pipeline4Publication/scripts/TrueWindows.py \
     --badcov /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA_BS.txt.gz \
     --te /Volumes/MartinResearch2/Wolf2019/references/seto_repeat_annotation.gff3 \
     --window 1000 \
     --step 1000 \
     --chromosomes /Volumes/MartinResearch2/Wolf2019/analyses/checkFASTA/coveragedist/ABBABABA_GoodCont.txt \
     --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/1k_windows.txt

cut -f 1-2,5,6,11,12,14,16,18,19,20,24 /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/1k_windows.txt-1000-1000.txt \
     >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/ABBABABA_1k_windows.txt

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/pipeline4Publication/scripts/PopGen-var.py \
     --input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_4FST.sync.gz \
     --pool-size 50,50,50,50,50,50,50,50,50,50 \
     --min-count 2 \
     --window 1000 \
     --step 1000 \
     --sitecount /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/ABBABABA_1k_windows.txt \
     --min-sites-frac 0.75 \
     --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/ABBABABA_1k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/pipeline4Publication/scripts/FST.py \
     --pool-size 50,50,50,50,50,50,50,50,50,50 \
     --input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_4FST.sync.gz \
     --minimum-count 2 \
     --minimum-cov 10 |
     gzip >/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA.fst.gz

CevennesCyn,EstoniaCyn,PetroiaCyn,SorenbergCyn,ZurichCyn,ABBAZurichNeo,CevennesNeo,GeschinenNeo,HospentalNeo,SorenbergNeo

python2.7 /Volumes/MartinResearch1/Inv_sequencing/pooled_karyos/scripts/binning_summary.py \
     --input /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA.fst.gz \
     --window-size 1000 \
     --output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA \
     --names CevennesCyn,EstoniaCyn,PetroiaCyn,SorenbergCyn,ZurichCyn,ABBAZurichNeo,CevennesNeo,GeschinenNeo,HospentalNeo,SorenbergNeo \
     --comp SorenbergCyn:SorenbergNeo,SorenbergCyn:PetroiaCyn,SorenbergCyn:EstoniaCyn,SorenbergNeo:GeschinenNeo,SorenbergNeo:HospentalNeo,PetroiaCyn:GeschinenNeo,ZurichCyn:PetroiaCyn,ZurichCyn:EstoniaCyn,ABBAZurichNeo:HospentalNeo,ZurichCyn:ABBAZurichNeo,ABBAZurichNeo:GeschinenNeo,SorenbergCyn:GeschinenNeo,ZurichCyn:GeschinenNeo \
     --ylim 0.5 \
     --stat 1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCPtC
mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCPhC
mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCPtC
mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCPhC
mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoC_Inter
mkdir /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuC_Inter

echo '''
library(tidyverse)
library(SuperExactTest)

## load data
DATA=na.omit(read.table("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA_binned.fst",
     header=T))

## combine Chrom and Pos to windows
DATA<-DATA %>%
     unite(Window,Chrom,Pos,sep=".") 

## calculate length of unique windows
N=length(unique(DATA$Window))

## get top 10% FST for each comparison 
DATA.top10 <- DATA %>%
     filter(Count>10) %>%
     group_by(Comp) %>%
     arrange(Comp, desc(FST)) %>% 
     filter(FST > quantile(FST, .99))

## SoC:PtC and SoC:PhC
DATA.SoC.intra <- DATA.top10%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,`SorenbergCyn:PetroiaCyn`,`SorenbergCyn:EstoniaCyn`)


## compare windows with Superexacttest
Comp.SoC<- list(
     "SoCPtC"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:PetroiaCyn"]]),]$Window,
     "SoCPhc"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:EstoniaCyn"]]),]$Window)

write.table(data.frame("SoCPtC"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:PetroiaCyn"]]),]$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCPtC/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

write.table(data.frame("SoCPhC"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:EstoniaCyn"]]),]$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCPhC/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

res=supertest(Comp.SoC,n=N)
write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCPtC_vs_SoCPhC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")

## Interspecific

DATA.SoC.inter <- na.omit(DATA %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     mutate(DeltaFST=`SorenbergCyn:SorenbergNeo` -`SorenbergCyn:GeschinenNeo`) %>%
     select(Window,DeltaFST) %>%
     arrange(Window, desc(DeltaFST))) %>% 
     filter(DeltaFST > quantile(DeltaFST, .99))

write.table(data.frame("Allo"=DATA.SoC.inter$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoC_Inter/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

Comp.SoC.Inter.PtC<- list(
     "Sym"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:PetroiaCyn"]]),]$Window,
     "Allo"=DATA.SoC.inter$Window)

res=supertest(Comp.SoC.Inter.PtC,n=N)

write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCInter_vs_SoCPtC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")

Comp.SoC.Inter.PhC<- list(
     "Sym"=DATA.SoC.intra[!is.na(DATA.SoC.intra[["SorenbergCyn:EstoniaCyn"]]),]$Window,
     "Allo"=DATA.SoC.inter$Window)

res=supertest(Comp.SoC.Inter.PhC,n=N)

write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/SoCInter_vs_SoCPhC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")

##ZuC:PtC andZuC:PhC
DATA.ZuC.intra <- DATA.top10%>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     select(Window,`ZurichCyn:PetroiaCyn`,`ZurichCyn:EstoniaCyn`)


## compare windows with Superexacttest
Comp.ZuC<- list(
     "ZuCPtC"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:PetroiaCyn"]]),]$Window,
     "ZuCPhc"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:EstoniaCyn"]]),]$Window)

write.table(data.frame("ZuCPtC"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:PetroiaCyn"]]),]$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCPtC/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

write.table(data.frame("ZuCPhC"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:EstoniaCyn"]]),]$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCPhC/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

res=supertest(Comp.ZuC,n=N)
write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCPtC_vs_ZuCPhC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")

## Interspecific

DATA.ZuC.inter <- na.omit(DATA %>%
     pivot_wider(names_from=Comp,values_from=FST) %>%
     mutate(DeltaFST=`ZurichCyn:ABBAZurichNeo`-`ZurichCyn:GeschinenNeo`) %>%
     select(Window,DeltaFST) %>%
     arrange(Window, desc(DeltaFST))) %>% 
     filter(DeltaFST > quantile(DeltaFST, .99))

write.table(data.frame("Allo"=DATA.ZuC.inter$Window),
     "/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuC_Inter/Candidates_FST.txt",
     row.names=F,
     col.names=F,
     quote = F)

Comp.ZuC.Inter.PtC<- list(
     "Sym"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:PetroiaCyn"]]),]$Window,
     "Allo"=DATA.ZuC.inter$Window)

res=supertest(Comp.ZuC.Inter.PtC,n=N)

write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCInter_vs_ZuCPtC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")

Comp.ZuC.Inter.PhC<- list(
     "Sym"=DATA.ZuC.intra[!is.na(DATA.ZuC.intra[["ZurichCyn:EstoniaCyn"]]),]$Window,
     "Allo"=DATA.ZuC.inter$Window)

res=supertest(Comp.ZuC.Inter.PhC,n=N)

write.table(data.frame("NAME"=summary(res)[1],"GENES"=summary(res)[8]),
     file="/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/ZuCInter_vs_ZuCPhC.set.txt",
     row.names=F,
     quote = F,
     sep="\t")


# DATA2 <- DATA%>%
#      pivot_wider(names_from=Comp,values_from=FST) %>%
#      select(-`SorenbergCyn:SorenbergNeo`,-`PetroiaCyn:GeschinenNeo`) %>%
#      mutate(DeltaFST=abs(`SorenbergCyn:PetroiaCyn` -`SorenbergNeo:GeschinenNeo`),Name="Mean_SymVsAllo_WithinSpecies")%>%
#      select(-`SorenbergCyn:PetroiaCyn`, -`SorenbergNeo:GeschinenNeo`)
# DATA=rbind(DATA1,DATA2)

# data_cum <- DATA %>% 
#   group_by(Chrom) %>% 
#   summarise(max_bp = max(Pos)) %>% 
#   mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
#   select(Chrom, bp_add)

# DATA <- DATA %>% 
#   inner_join(data_cum, by = "Chrom") %>% 
#   mutate(bp_cum = Pos + bp_add)

# axis_set <- DATA %>% 
#   group_by(Chrom) %>% 
#   summarize(center = mean(bp_cum))

#  PLOT <- ggplot(DATA, aes(x = bp_cum, y = DeltaFST, 
#                                   color = as_factor(Chrom), size = DeltaFST)) +
#   facet_grid(Name~.)+
#   geom_vline(data = . %>% filter(DeltaFST > 0.75),
#              aes(xintercept = bp_cum))+
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set$Chrom, breaks = axis_set$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
#   scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$Chrom)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "delta FST") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
#   )

# PLOT

# ggsave("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA_binned.png",
#      PLOT,
#      width=16,
#      height=6)

''' >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/FST_cand.r

Rscript /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/FST_cand.r

### prepare for PoolHMM
mkdir -p /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups

samples=(CevennesCyn EstoniaCyn PetroiaCyn SorenbergCyn ZurichCyn ABBAZurichNeo CevennesNeo GeschinenNeo HospentalNeo SorenbergNeo)
pos=(3 4 9 10 12 14 16 17 18 22)

samples=(SorenbergNeo HospentalNeo)
pos=(22 18)

for index in ${!samples[@]}; do
     Start=$((pos[index] * 3 + 1))
     End=$((pos[index] * 3 + 3))
     ID=${samples[index]}

     echo $Start $End $ID

     gunzip -c /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA.mpileup.gz |
          cut -f1-3,${Start}-${End} |
          gzip >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups/${ID}.pileup.gz
done

### make PoolHMM
comparisons=(ZuC_Inter ZuCPhC ZuCPtC SoC_Inter SoCPhC SoCPtC)
Cyn=(ZurichCyn ZurichCyn ZurichCyn SorenbergCyn SorenbergCyn SorenbergCyn)
Neo=(ABBAZurichNeo ABBAZurichNeo ABBAZurichNeo SorenbergNeo SorenbergNeo SorenbergNeo)
PrefixCyn=(ZuC ZuC ZuC SoC SoC SoC)
PrefixNeo=(ZuN ZuN ZuN SoN SoN SoN)

for index in ${!comparisons[@]}; do
     CO=${comparisons[index]}
     CY=${Cyn[index]}
     NE=${Neo[index]}
     PC=${PrefixCyn[index]}
     PN=${PrefixNeo[index]}

     python3 /Volumes/MartinResearch2/Wolf2019/scripts/Pileup4PoolHMM.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups/${CY}.pileup.gz \
          --windows /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST.txt \
          --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/${PC} \
          --windowsize 1000

     python3 /Volumes/MartinResearch2/Wolf2019/scripts/Pileup4PoolHMM.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups/${NE}.pileup.gz \
          --windows /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST.txt \
          --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/${PN} \
          --windowsize 1000

     for i in /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/*.pileup; do

          cd /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}

          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /Volumes/MartinResearch2/Wolf2019/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix ${ID} \
               -n 50 \
               --only-spectrum \
               --theta 0.005 \
               -r 100

          python2 /Volumes/MartinResearch2/Wolf2019/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix ${ID} \
               -n 50 \
               --pred \
               -k 0.000000000000001 \
               -s ${ID} \
               -e sanger
     done

     # now random dataset

     WC=($(wc -l /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST.txt))
     N=${WC[0]}

     awk 'NR>1' /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/FST/ABBA_BABA_binned.fst |
          sort -R |
          head -${N} |
          awk '{print $1"."int($2)}' \
               >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_FST.txt

     python3 /Volumes/MartinResearch2/Wolf2019/scripts/Pileup4PoolHMM.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups/${CY}.pileup.gz \
          --windows /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_FST.txt \
          --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_${PC} \
          --windowsize 1000

     python3 /Volumes/MartinResearch2/Wolf2019/scripts/Pileup4PoolHMM.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/pileups/${NE}.pileup.gz \
          --windows /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_FST.txt \
          --output /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_${PN} \
          --windowsize 1000

     for i in /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Rand_*.pileup; do

          cd /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}

          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /Volumes/MartinResearch2/Wolf2019/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix ${ID} \
               -n 50 \
               --only-spectrum \
               --theta 0.005 \
               -r 100

          python2 /Volumes/MartinResearch2/Wolf2019/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py \
               --prefix ${ID} \
               -n 50 \
               --pred \
               -k 0.000000000000001 \
               -s ${ID} \
               -e sanger
     done

done

### get FASTA sequences

for index in ${!comparisons[@]}; do
     CO=${comparisons[index]}
     CY=${Cyn[index]}
     NE=${Neo[index]}
     PC=${PrefixCyn[index]}
     PN=${PrefixNeo[index]}

     python /Volumes/MartinResearch2/Wolf2019/scripts/SummarizePoolHMM.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO} \
          --windowsize 1000 \
          >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_sweeps.txt

done

### make PoolHMM

for index in ${!comparisons[@]}; do
     CO=${comparisons[index]}

     python3 /Volumes/MartinResearch2/Wolf2019/scripts/FASTA4windows.py \
          --input /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST.txt \
          --FASTA /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
          --windowsize 1000 \
          >/Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST.fasta

done

### summarize PoolHMM

### BLAST candidate windows
comparisons=(SoC_Inter SoCPhC SoCPtC ZuC_Inter ZuCPhC ZuCPtC)

for index in ${!comparisons[@]}; do
     CO=${comparisons[index]}

     sh /Volumes/MartinResearch2/Wolf2019/shell/doBlast_ABBA.sh \
          /Volumes/MartinResearch2/Wolf2019/analyses/PopGen/tho/sweeps/${CO}/Candidates_FST

done
