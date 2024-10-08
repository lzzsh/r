#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=128GB
#SBATCH -J GRNboost
#SBATCH -o GRNboost.%j.out
#SBATCH -e GRNboost.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn
#conda activate pyscenic

chmod +x 3.1_arboreto_with_multiprocessing.py
./3.1_arboreto_with_multiprocessing.py \
    exprMat_filtered_TB.loom \
    ../../cisTarget_databases/mm10/allTFs_mm.txt \
    --method grnboost2 \
    --output adj.tsv \
    --num_workers 15 \
    --seed 123

Rscript 3.2_generate_regulon.R