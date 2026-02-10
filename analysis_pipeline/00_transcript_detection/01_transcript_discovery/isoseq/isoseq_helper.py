# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/1 14:27
@Email   : xuechao@szbl.ac.cn
@Desc    :
"""
import re
import sys
import os
import argparse
import logging
import polars as pl
from statistics import median
import numpy as np
import pandas as pd
# import pysam


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def summarize_bam_stats(input_path, output_prefix):
    def safe_median(x):
        return median(x) if len(x) > 0 else 0
    output_summary=f"{output_prefix}.bam_summary.tsv"
    output_length=f"{output_prefix}.mapped_length.tsv"
    output_qscore=f"{output_prefix}.read_Qscore.tsv"
    summary_rows = []
    bam=input_path
    sample = os.path.basename(bam).split('.')[0]
    print(f"Processing {sample} ...")
    bamfile = pysam.AlignmentFile(bam, "rb")
    total_reads = 0
    mapped_reads = 0
    lengths_total = []
    lengths_mapped = []
    qscore_total = []
    qscore_mapped = []
    mapqs_mapped = []
    for read in bamfile.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        total_reads += 1
        length = read.query_length or 0
        lengths_total.append(length)
        if read.query_qualities:
            q = np.mean(read.query_qualities)
        else:
            q = 0
        qscore_total.append(q)
        if not read.is_unmapped:
            mapped_reads += 1
            lengths_mapped.append(length)
            qscore_mapped.append(q)
            mapqs_mapped.append(read.mapping_quality)
    bamfile.close()
    summary_rows.append({
        "sample": sample,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "total_median_length": safe_median(lengths_total),
        "total_min_length": np.nanmin(lengths_total),
        "total_max_length": np.nanmax(lengths_total),
        "mapped_median_length": safe_median(lengths_mapped),
        "mapped_min_length": np.nanmin(lengths_mapped),
        "mapped_max_length": np.nanmax(lengths_mapped),
        "total_median_Qscore": safe_median(qscore_total),
        "mapped_median_Qscore": safe_median(qscore_mapped),
        "mapped_median_MAPQ": safe_median(mapqs_mapped),
        "mapping_rate(%)": round(mapped_reads / total_reads * 100, 4) if total_reads > 0 else 0
    })
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(output_summary, index=False)
    print(f"✅ Summary saved: {output_summary}")
    with open(output_length, "w") as f:
        f.write('\n'.join([f'{x}' for x in lengths_total]))
        print(f"✅ Length saved: {output_length}")
    with open(output_qscore, "w") as f:
        f.write('\n'.join([f'{x}' for x in qscore_total]))
        print(f"✅ Read Qscore saved: {output_qscore}")
    pass

def __count_FLNC(sample_file, mapping_df):
    sample_name = os.path.basename(sample_file).split(".")[0]
    sample_df = pl.read_csv(
        sample_file,
        has_header=False,
        new_columns=["read_id"],
        separator="\t",
        infer_schema_length=1000
    )
    joined = sample_df.join(mapping_df, on="read_id", how="inner")
    grouped = (
        joined.group_by("isoform")
        .agg(pl.len())
        .rename({"len": sample_name})
    )
    return grouped, sample_name

def count_FLNC_per_sample(out_root_dir,chr):
    chr_read_stat_path = f'{out_root_dir}/chr_based/{chr}/collapse/merged.collapsed.read_stat.txt'
    out_path = f'{out_root_dir}/chr_based/{chr}/collapse/merged.collapsed.FLNC.count.txt'
    sample_files = [f'{out_root_dir}/sample_based/{s}/{s}.flnc.read_id'
                    for s in os.listdir(f'{out_root_dir}/sample_based')]
    logger.info(f"Chr: {chr}")
    logger.info("Loading mapping file into Polars...")
    mapping_df = pl.read_csv(
        chr_read_stat_path,
        has_header=False,
        new_columns=["read_id", "isoform"],
        separator="\t",
        infer_schema_length=1000
    )
    logger.info(f"Mapping rows = {mapping_df.height}")

    merged = None
    for sample_file in sample_files:
        logger.info(f'start {os.path.basename(sample_file)}')
        grouped, sample_name = __count_FLNC(sample_file, mapping_df)
        if merged is None:
            merged = grouped
        else:
            logger.info(f"Merging sample {sample_name}...")
            merged = merged.join(grouped, on="isoform", how="full", coalesce=True)

    merged = merged.fill_null(0).sort("isoform")
    logger.info(f"Writing output: {out_path}")
    merged.write_csv(out_path)
    logger.info("Done.")

def filter_isoforms(in_prefix, out_prefix, min_all_read_count=10):
    min_all_read_count=int(min_all_read_count)
    min_n_tissue = 1
    min_read_per_sample = 1
    df = pd.read_csv(f"{in_prefix}.flnc_count.txt", index_col=0)
    raw_n_isoform=df.shape[0]
    stats_df = pd.DataFrame(index=df.index)
    stats_df['n_samples'] = (df>=min_read_per_sample).sum(axis=1)
    stats_df['n_total'] = df.sum(axis=1)
    keep = df[(stats_df['n_total']>=min_all_read_count)&(stats_df['n_samples']>=min_n_tissue)]
    keep_ids = set(keep.index)
    keep.to_csv(f"{out_prefix}.flnc_count.txt", sep='\t')
    logger.info(f'raw isoform: {raw_n_isoform}; filtered isoform: {keep.shape[0]}; '
        f'MIN_READ={min_all_read_count}; MIN_N={min_n_tissue}; MIN_READ_PER_SAMPLE={min_read_per_sample}')
    with open(f"{in_prefix}.gff") as gff_in, open(f"{out_prefix}.gff",'w') as gff_out:
        for line in gff_in:
            # if line.startswith('#') or (m:=re.search(r'transcript_id "([^"]+)"', line)) and m.group(1) in keep_ids:
            if line.startswith('#'):
                gff_out.write(line)
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m and m.group(1) in keep_ids:
                gff_out.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Help programs for isoseq pipeline',
    )
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands'
    ) #        required=True,
    subparsers.required = True

    # QC
    parser_qc = subparsers.add_parser(
        'qc',
        help='FLNC quality control'
    )
    parser_qc.add_argument(
        '-i', '--input',
        required=True,
        type=str,
        help='BAM path'
    )
    parser_qc.add_argument(
        '-o', '--output-prefix',
        required=True,
        type=str,
        help='Output prefix'
    )

    # count FLNC read
    parser_count = subparsers.add_parser(
        'count',
        help='Count supported read of isoforms for each sample'
    )
    parser_count.add_argument(
        '-o', '--out-root-dir',
        required=True,
        type=str,
        help='Output root directory (required)'
    )
    parser_count.add_argument(
        '-c', '--chr',
        required=True,
        type=str,
        help='Chromosome identifier (required)'
    )

    # pre-filtering isoform to run annotation
    parser_filter = subparsers.add_parser(
        'filter_isoform',
        help='pre-filtering collapsed isoform to run annotation'
    )
    parser_filter.add_argument(
        '-i', '--in-prefix',
        required=True,
        type=str,
        help='Output root directory (required)'
    )
    parser_filter.add_argument(
        '-o', '--out-prefix',
        required=True,
        type=str,
        help='Chromosome identifier (required)'
    )
    parser_filter.add_argument(
        '-m', '--min-count',
        required=True,
        type=str,
        help='Chromosome identifier (required)'
    )

    args = parser.parse_args()
    if args.command == 'qc':
        summarize_bam_stats(args.input, args.output_prefix)
    elif args.command == 'count':
        count_FLNC_per_sample(args.out_root_dir, args.chr)
    elif args.command == 'filter_isoforms':
        filter_isoforms(args.in_prefix, args.out_prefix, args.min_count)

    else:
        logger.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()