#!/bin/bash
#SBATCH -c 16
#SBATCH -J jobname
#SBATCH -o jobname.%j.out
#SBATCH -e jobname.%j.err
#SBATCH --mem=128G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liuyifang@westlake.edu.cn

# Reset the SECONDS variable
SECONDS=0

# Activate the specific conda environment for R
source /home/liuxiaodongLab/liuyifang/miniconda3/etc/profile.d/conda.sh
conda activate R-4.4.0
which R

# Run the R script (assuming it's in the same folder as this shell script)
Rscript jobname.R

# Calculate the total running time in hours
total_time=$SECONDS
hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total running time: ${hours} hours"