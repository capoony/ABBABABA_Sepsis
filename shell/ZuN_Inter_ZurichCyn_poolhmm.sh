

    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN_Inter/ZurichCyn
     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py           --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/ZurichCyn.pileup.gz           --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN_Inter/ZuN_Inter_Candidates_FST.txt           --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN_Inter/ZurichCyn/ZurichCyn           --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN_Inter/ZurichCyn
     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/ZuN_Inter/ZurichCyn/ZurichCyn*.pileup; do


          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --only-spectrum                --theta 0.005                -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --pred                -k 0.000000000000001                -s ${ID}                -e sanger
     done

    
