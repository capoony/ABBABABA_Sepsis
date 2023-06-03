

    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/RAND_ABBAZurichNeo

     # now random dataset

     WC=($(wc -l /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/ZuN.GeN_mean_Candidates_FST.txt))
     N=${WC[0]}

     sort -R /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt |
        head -${N}          >/media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/Rand_FST.txt

     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py           --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/ABBAZurichNeo.pileup.gz           --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/Rand_FST.txt           --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/RAND_ABBAZurichNeo/RAND_ABBAZurichNeo           --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/RAND_ABBAZurichNeo

     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN.GeN_mean/RAND_ABBAZurichNeo/RAND_ABBAZurichNeo_*.pileup; do

          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --only-spectrum                --theta 0.005                -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --pred                -k 0.000000000000001                -s ${ID}                -e sanger
     done

     
