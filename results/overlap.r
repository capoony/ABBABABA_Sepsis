

    library(SuperExactTest)

    D1=read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoN_Inter/SoN_Inter_Candidates_FST.txt')
    D2=read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/ZuN_Inter/ZuN_Inter_Candidates_FST.txt')
    N=nrow(read.table('/media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt'))

    Comp<- list(
     SoN_Inter=D1[,1],
     ZuN_Inter=D2[,1])

    res=supertest(Comp,n=N)

    write.table(data.frame('NAME'=summary(res)[1],
            'GENES'=summary(res)[8]),
            file='/media/inter/mkapun/projects/ABBABABA_Sepsis/results/overlap/SoNIZuNI.txt',
            row.names=F,
            quote = F,
            sep='\t')

    
