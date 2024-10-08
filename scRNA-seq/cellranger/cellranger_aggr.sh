#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64GB
#SBATCH -J AGG_9_to_18
#SBATCH -o AGG_9_to_18.%j.out
#SBATCH -e AGG_9_to_18.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn


cellranger aggr --id=AGG_9_to_18 --csv=AGG_9_to_18.csv