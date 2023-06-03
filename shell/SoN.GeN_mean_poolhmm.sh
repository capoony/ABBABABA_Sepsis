


    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/SoN.GeN_mean/SorenbergNeo
     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py           --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/SorenbergNeo.pileup.gz           --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/SoN.GeN_mean/SoN.GeN_mean_Candidates_FST.txt           --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/SoN.GeN_mean/SorenbergNeo/SorenbergNeo           --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/SoN.GeN_mean/SorenbergNeo
     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/SoN.GeN_mean/SorenbergNeo/SorenbergNeo*.pileup; do


          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --only-spectrum                --theta 0.005                -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --pred                -k 0.000000000000001                -s ${ID}                -e sanger
     done

    
