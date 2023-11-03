

    #install.packages('poolfstat')
    library('poolfstat')
    library(tidyverse)
    library(SuperExactTest)

    ##### Convert sync file to poolfstat file, and call SNPS

    # We first have to give haploid sizes of each pool. Here, I had mostly 40 individuals per pool, and since I am working with diploid species, we multiply that by 2,  to get 80 individuals for most pools.
    psizes <- as.numeric(c(50,50,50,10))
    # Then we give the names of each pool/sample.
    pnames <- as.character(c('GeN','SoN','SoC','Sor'))

    SG.pooldata <- popsync2pooldata(sync.file = '/media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4_new/SoN_GeN.sync.gz', 
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
    file='/media/inter/mkapun/projects/ABBABABA_Sepsis/results/f4_new/SoN_GeN_f4.stat',
    quote=F)

    
