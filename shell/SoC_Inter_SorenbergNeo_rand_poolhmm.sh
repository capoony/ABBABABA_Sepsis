

    mkdir /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergNeo

     # now random dataset

     WC=($(wc -l /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/SoC_Inter_Candidates_FST.txt))
     N=${WC[0]}

     sort -R /media/inter/mkapun/projects/ABBABABA_Sepsis/results/Windows.txt |
        head -${N}          >/media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/Rand_FST.txt

     python3 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/Pileup4PoolHMM.py           --input /media/inter/mkapun/projects/ABBABABA_Sepsis/data/SorenbergNeo.pileup.gz           --windows /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/Rand_FST.txt           --output /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergNeo/RAND_SorenbergNeo           --windowsize 1000

    cd /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergNeo

     for i in /media/inter/mkapun/projects/ABBABABA_Sepsis/results/poolhmm/SoC_Inter/RAND_SorenbergNeo/RAND_SorenbergNeo_*.pileup; do

          tmp=${i##*/}
          ID=${tmp%.pileup*}

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --only-spectrum                --theta 0.005                -r 100

          python2 /media/inter/mkapun/projects/ABBABABA_Sepsis/scripts/PoolHMM-1.4.4-custom.1/pool-hmm.py                --prefix ${ID}                -n 50                --pred                -k 0.000000000000001                -s ${ID}                -e sanger
     done

     
