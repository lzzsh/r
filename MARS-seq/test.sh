#!/bin/bash

# load the necessary tools and environment
module load /soft/modules/modulefiles/bioinfo/star/2.7.9a
module load /soft/modules/modulefiles/bioinfo/samtools/1.16.1
module load /soft/modules/modulefiles/bioinfo/sratoolkit/2.11.2

# split fastq files from sra
cd /storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/rawdata
for file in $(ls *.sra);do
        prefix=${file%%.*}
        fastq-dump --split-3 ${prefix}.sra
        mv ${prefix}_1.fastq ${prefix}_r1.fq
        mv ${prefix}_2.fastq ${prefix}_r2.fq
done
gzip *.fq

# convert mars-seq into 10x Genomics
cd /storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/convert_10x
python /storage/zhangkaiLab/liaozizhuo5001/project/scrna_TE/velocity.py convert --input ../rawdata --output ./converted --threads 8

# cut adapter
for file in $(ls ./converted/*_r1.fq.gz);do
        prefix=${file%%_*}
        prefix=$(basename $prefix)
        cutadapt -m 20 -A 'T{68}' --pair-filter=both -o ./cutadapt/${prefix}_R1.trimmed.fq.gz -p ./cutadapt/${prefix}_R2.trimmed.fq.gz \
        ./converted/${prefix}_r1.fq.gz ./converted/${prefix}_r2.fq.gz
        sleep 1
done

# starsolo
cd /storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/convert_10x/starsolo
for file in $(ls ../cutadapt/*_R1.trimmed.fq.gz);do
        prefix=${file%%_*}
        prefix=$(basename $prefix)
        STAR --genomeDir /storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/genome/gencode_vM23 --readFilesCommand zcat  \
        --readFilesIn ../cutadapt/${prefix}_R2.trimmed.fq.gz ../cutadapt/${prefix}_R1.trimmed.fq.gz \
        --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 7 --soloUMIstart 8 --soloUMIlen 8 \
        --outSAMtype BAM SortedByCoordinate --soloFeatures Gene Velocyto \
        --soloCBwhitelist ../whitelist/${prefix}.txt --outSAMattributes NH HI AS nM CR CY UR UY \
        --outFileNamePrefix ${prefix}_
        mv ${prefix}_Aligned.sortedByCoord.out.bam ${prefix}.bam
        sleep 1
done

# divide bam files into embryo units
cd /storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/convert_10x/extract_embryo
for file in $(ls ../starsolo/*.bam);do
	prefix=${file%%.bam}
	prefix=$(basename $prefix)
	python extract_embryo.py --input ../starsolo/${prefix} --output .
	sleep 1
done


cd storage/zhangkaiLab/liaozizhuo5001/data/scrna_TE/convert_10x/scTE
for file in $(ls ../starsolo/*.bam);do
	prefix=${file%%.*}
	prefix=$(basename $prefix)
	scTE -i ../starsolo/${prefix}.bam -o ${prefix} -x mm10.exclusive.idx  True -CB CR -UMI UR
	sleep 1
done

