#!/bin/bash 
#SBATCH --job-name=0.9compare
#SBATCH --partition=super
#SBATCH -N 16 #nodes
#SBATCH -n 126 #tasks
#SBATCH -c 1 # 1 core per task
#SBATCH --time=05-00:00:00
#SBATCH --output=0.9compare_compare.%j.out
#SBATCH --error=0.9compare_compare.%j.err

ARGS=(0.9 0.8 0.7 0.6 0.5 0.2)
for file in *fasta
do
    for i in {0..5}
    do
        for split in ${ARGS[i]}
            do srun -N 1 -n 1 compare.sh -t 511145 -f $file -d $split &
        done
    done
done
wait
