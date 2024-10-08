#!/bin/bash

for i in {9..17}
do
cat <<EOT > E${i}_5.sh
#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=128GB
#SBATCH -J E${i}_5
#SBATCH -o E${i}_5.%j.out
#SBATCH -e E${i}_5.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn


cellranger count --id=sample${i} \\
		 --transcriptome=/storage/liuxiaodongLab/liaozizhuo/data/refdata-gex-mm10-2020-A \\
		 --fastqs=/storage/liuxiaodongLab/liuyifang/Projects/YutingFu/Placenta_Project/rawData/20240102-CQT2023120601-F001-rawdata_fastqs/Sample_E${i}_5/ \\
                 --create-bam=true \\
                 --sample=E${i}_5
EOT
done