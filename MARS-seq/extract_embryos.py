#!/usr/bin/env python
import pysam
import pandas as pd
import argparse
import logging

def extract_reads(SRR_sample, output_bam):
    # 读取样本信息
    amplication_batch = pd.read_csv("amplication_batch.csv")
    metadata = pd.read_table("GSE267870_cell_metadata_mars_seq.tsv", sep="\t")
    pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)

    # 合并信息
    merged_data = amplication_batch.merge(metadata, left_on='Amp_batch_ID', right_on='Plate')
    selected_columns = ['Amp_batch_ID', 'embryo', 'well_barcode', 'GSM', 'SRR']
    cr_tag_values = merged_data[selected_columns]
    cr_tag_values = cr_tag_values[cr_tag_values['embryo'] != "empty"]
    # print(cr_tag_values)

    # 获取胚胎和barcode对应信息
    embryo_info = cr_tag_values.groupby('Amp_batch_ID')['embryo'].agg(lambda x: list(set(x)))
    embryo_info = embryo_info.reset_index()
    embryo_info = embryo_info.merge(amplication_batch[['Amp_batch_ID', 'GSM', 'SRR']], left_on='Amp_batch_ID', right_on='Amp_batch_ID')

    # 将结果转换为DataFrame
    embryo_info_df = embryo_info.reset_index()
    # print(embryo_info_df)

    # 打开输入BAM文件
    input_bamfile = pysam.AlignmentFile(f"{SRR_sample}.bam", "rb")

    # 在embryo_info DataFrame中查找与SRR样本相关的行
    sample_info = embryo_info[embryo_info['SRR'] == SRR_sample]

    # 如果没有找到匹配项，则返回错误消息
    if sample_info.empty:
        print(f"No information found for SRR sample {SRR_sample}")
        return

    sample_info = sample_info.explode('embryo')

    # 获取胚胎信息
    embryos = sample_info['embryo']

    for i in embryos:
        cr_tag_values = cr_tag_values[cr_tag_values['embryo'] == i]

        # 获取cell barcode列表
        cr_tag_values = cr_tag_values['well_barcode']

        # 创建输出BAM文件
        output_bamfile = pysam.AlignmentFile(f"{output_bam}/{i}.bam", "wb", template=input_bamfile)

        # 遍历BAM文件中的读取记录
        for read in input_bamfile:
            # 检查读取的tags属性
            tags = dict(read.tags)

            # 检查是否包含指定的CR:Z标签
            if "CR:Z" in tags and tags["CR:Z"] in cr_tag_values:
                # 将匹配的读取写入输出BAM文件
                output_bamfile.write(read)

        # 关闭BAM文件
        input_bamfile.close()
        output_bamfile.close()

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s  %(message)s", "%d-%m-%Y %H:%M:%S")

    arg_parser = argparse.ArgumentParser(description="extract reads from different embryos.")

    # embryo extraction
    arg_parser.add_argument("--input", type=str, help="Input folder with bam", required=True)
    arg_parser.add_argument("--output", type=str, help="output folder", required=True)
    arg_parser.add_argument("--version", "-v", action="version", version="embryos 0.1")

    args = arg_parser.parse_args()
    extract_reads(args.input, args.output)

# Example command
# extract_reads --input SRR29084300 --output .