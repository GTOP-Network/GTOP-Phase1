# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/30 11:07
@Desc    : Run Iso-Seq pipeline in HPC or Single node.
"""
import logging
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import os
import subprocess
import sys

import pandas as pd
import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# config for computer evn
REF_DIR='/lustre/home/cxue/raw_data/GMTiP/ref/LRS'

PRIMER_PATH=f'{REF_DIR}/IsoSeq_v2_primers_12.fasta'
REF_GENOME_INDEX=f'{REF_DIR}/genome_ISOSEQ.mmi'
CAGE_PEAK = f"{REF_DIR}/human.refTSS_v3.1.hg38.bed"
POLYA = f"{REF_DIR}/mouse_and_human.polyA_motif.txt"
REF_GENOME_FASTA = f"{REF_DIR}/genome.fa"
REF_GTF = f"{REF_DIR}/gencode.v47.annotation.gtf"
REF_SORTED_GTF = f"{REF_DIR}/gencode.v47.annotation.sorted.gtf.gz"
REF_GTF_FASTA = f"{REF_DIR}/gencode.v47.annotation.fa"

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
ROOT_OUTPUT_DIR=ARGS[2]
LOG_NAME = ARGS[3]
N_TASK = int(ARGS[4])
NT_PER_TASK = int(ARGS[5])

LOAD_ISOSEQ_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate isoseq'
LOAD_PBINDEX_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate pbtk'
LOAD_SAMTOOLS_ENVS_CMD='module load samtools'
LOAD_SQANTI3_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate sqanti3'
LOAD_isoLASER_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate isoLASER'
LOAD_PICARD_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate picard_env'
LOAD_SALMON_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate salmon_env'
LOAD_POLARS_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate polars_env'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'

helper_py=f'{LOAD_PYSAM_ENV_CMD} && python {CURRENT_DIR}/isoseq_helper.py'

def make_dir(*dirname):
    for sdir in dirname:
        os.makedirs(sdir,exist_ok=True)

def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')

def run_commands_threadpool(cmd_list, main_log, max_workers=4,main_log_prefix='batch'):
    '''

    :param cmd_list:
        [
            [{'cmd':'ls -l', 'log_path': 'ls_log/1.log'},
            {'cmd':'df -h', 'log_path': 'df_log/1.log'}],
        ]
    :param main_log:
    :param max_workers:
    :param main_log_prefix:
    :return:
    '''
    make_dir(os.path.dirname(main_log))
    log(f"[INFO] Batch start. Total tasks: {len(cmd_list)}, Max workers: {max_workers}", main_log)
    def worker(task_cmds):
        for cmd_it in task_cmds:
            cmd=cmd_it['cmd']
            # cmd=' '.join(cmd.split())
            task_log=cmd_it['log_path']
            logs_dir=os.path.dirname(task_log)
            make_dir(logs_dir)
            main_log_path=f'{logs_dir}/{main_log_prefix}.log'
            # if isinstance(cmd, str):
            #     cmd_str = cmd
            #     cmd_exec = cmd.split()
            # else:
            #     cmd_exec = cmd
            #     cmd_str = " ".join(cmd)
            log(f"[INFO] Start task: {' '.join(cmd.split())}", main_log_path)
            try:
                with open(task_log, "w", encoding="utf-8") as lf:
                    process = subprocess.Popen(
                        cmd,
                        shell=True,
                        stdout=lf,
                        stderr=subprocess.STDOUT,
                        cwd=logs_dir
                    )
                    retcode = process.wait()
                if retcode == 0:
                    log(f"[INFO] Finished task: {cmd}", main_log_path)
                else:
                    log(f"[ERROR] Task failed ({retcode}): {cmd}", main_log_path)
            except Exception as e:
                log(f"[EXCEPTION] Task crashed: {cmd}\n{e}", main_log_path)

    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cmds in cmd_list:
            futures.append(executor.submit(worker, cmds))
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log(f"[ERROR] Exception in worker thread: {e}", main_log)
    log(f"[INFO] Batch finished. Logs in {main_log}.")

def lima_cmd(bam_path, sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
        {LOAD_ISOSEQ_ENVS_CMD} && 
        lima
        {bam_path} 
        {PRIMER_PATH} 
        {out_dir}/{sample_id}.fl.bam
        --isoseq 
        --num-threads {nt} 
        --dump-clips 
        --store-unbarcoded 
        --peek-guess 
        --log-file {out_dir}/{sample_id}.lima_inner.log
    '''
    cmd=' '.join(cmd.split())
    return cmd

def refine_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
        {LOAD_ISOSEQ_ENVS_CMD} && 
        bam_path=$(find {out_dir} -maxdepth 1 -type f -name "{sample_id}.fl.*.bam" ! -name "*.unbarcoded.bam" | head -n 1) &&
        isoseq refine 
        "$bam_path" 
        {PRIMER_PATH} 
        {out_dir}/{sample_id}.flnc.bam
        --require-polya 
        --num-threads {nt}" 
    '''
    cmd=' '.join(cmd.split())
    return cmd

def refine_with_ployA_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
        {LOAD_ISOSEQ_ENVS_CMD} && 
        bam_path=$(find {out_dir} -maxdepth 1 -type f -name "{sample_id}.fl.*.bam" ! -name "*.unbarcoded.bam" | head -n 1) &&
        isoseq refine 
        "$bam_path" 
        {PRIMER_PATH} 
        {out_dir}/{sample_id}.flnc_polyA.bam
        --num-threads {nt} && rm -rf "$bam_path" 
    '''
    cmd=' '.join(cmd.split())
    return cmd

def cluster_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
    {LOAD_ISOSEQ_ENVS_CMD} && 
    isoseq cluster2 {out_dir}/{sample_id}.flnc.bam {out_dir}/{sample_id}.clustered.bam 
    --sort-threads {nt} 
    --num-threads {nt}
    '''
    cmd = ' '.join(cmd.split())
    return cmd

def pbmm2_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
    {LOAD_ISOSEQ_ENVS_CMD} && 
    pbmm2 align {REF_GENOME_INDEX} {out_dir}/{sample_id}.clustered.bam {out_dir}/{sample_id}.mapped.bam 
      --preset ISOSEQ --unmapped --sort -j {nt}
      --sample {sample_id} --log-level INFO
    '''
    cmd = ' '.join(cmd.split())
    return cmd

def flnc_pbmm2_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    cmd=f'''
    {LOAD_ISOSEQ_ENVS_CMD} && 
    pbmm2 align {REF_GENOME_INDEX} {out_dir}/{sample_id}.flnc.bam {out_dir}/{sample_id}.flnc_mapped.bam 
      --preset ISOSEQ --unmapped --sort -j {nt}
      --sample {sample_id} --log-level INFO
    '''
    cmd = ' '.join(cmd.split())
    return cmd

def split_chr_samtools_cmd_bak(sample_id, out_dir, nt):
    chr_split_dir=f'{out_dir}/{sample_id}.chr_bam'
    make_dir(chr_split_dir)
    cmd = f"""
    {LOAD_SAMTOOLS_ENVS_CMD} && 
    for chr in $(samtools idxstats {out_dir}/{sample_id}.mapped.bam | cut -f1 | grep -v '*'); do
        echo "[`date '+%F %T'`] processing $chr ..."
        samtools view -@ {nt} -b {out_dir}/{sample_id}.mapped.bam  $chr > {chr_split_dir}/$chr.bam
        samtools index {chr_split_dir}/$chr.bam
    done
    """
    return cmd


def split_chr_samtools_cmd(sample_id, out_dir, nt, bam_prefix='mapped', bam_dir_suffix='chr_bam'):
    chr_split_dir = f"{out_dir}/{sample_id}.{bam_dir_suffix}"
    make_dir(chr_split_dir)
    cmd = f"""
    {LOAD_SAMTOOLS_ENVS_CMD} && \\
    mkdir -p {chr_split_dir} && \\
    echo "[`date '+%F %T'`] Starting chromosome split for {sample_id} ..." && \\

    normal_chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM" && \\
    others_list={chr_split_dir}/others.list && > $others_list && \\

    for chr in $(samtools idxstats {out_dir}/{sample_id}.{bam_prefix}.bam | cut -f1 | grep -v '*'); do \\
        if echo "$normal_chrs" | grep -wq "$chr"; then \\
            echo "[`date '+%F %T'`] processing $chr ..." && \\
            samtools view -@ {nt} -b {out_dir}/{sample_id}.{bam_prefix}.bam $chr > {chr_split_dir}/$chr.bam && \\
            samtools index {chr_split_dir}/$chr.bam; \\
        else \\
            echo "$chr" >> $others_list; \\
        fi; \\
    done && \\

    if [ -s $others_list ]; then \\
        echo "[`date '+%F %T'`] merging non-standard contigs into others.bam ..." && \\
        samtools view -@ {nt} -b {out_dir}/{sample_id}.{bam_prefix}.bam $(cat $others_list) > {chr_split_dir}/others.bam && \\
        samtools index {chr_split_dir}/others.bam; \\
    fi && \\

    echo "[`date '+%F %T'`] done."
    """
    return cmd


def hifi_to_mapped_bam(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    '''
    Sample based multiple task.
    From HiFi read (.bam file) to chr-seperated mapped genome (.bam file).
    :return:
    '''
    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/sample_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    for i in df.index:
        series_cmds=[]
        bam_path=df.loc[i,'bam_path']
        sample_id=df.loc[i,'sample_id']
        # limma
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=lima_cmd(bam_path, sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_lima.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # isoseq refine
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=refine_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_refine.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # isoseq refine, retain polyA seq for IsoQuant (an alternative tools).
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=refine_with_ployA_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_refine_polyA.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # extract read id from FLNC.bam
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd = (f"{LOAD_SAMTOOLS_ENVS_CMD} && samtools view -@ {nt_per_task} {out_dir}/{sample_id}.flnc.bam |"
               f" awk '{{print $1}}' > {out_dir}/{sample_id}.flnc.read_id")
        log_path=f'{out_dir}/{sample_id}.run_extract_flnc_read.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # map FLNC read to genome
        cmd=flnc_pbmm2_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_pbmm2.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # FLNC QC
        bam_path=f'{out_dir}/{sample_id}.flnc_mapped.bam'
        out_prefix=f'{out_dir}/{sample_id}.flnc_mapped.QC_stat'
        cmd = f"{helper_py} qc -i {bam_path} -o {out_prefix}"
        log_path = f'{out_dir}/{sample_id}.run_flnc_QC.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # isoseq cluster
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=cluster_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_cluster.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # pbmm2 mapping
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=pbmm2_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_pbmm2.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # splitting chr by samtools
        out_dir=f'{OUTPUT_DIR}/{sample_id}'
        cmd=split_chr_samtools_cmd(sample_id, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_split_chr.log'
        # series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)
    pass


def collapse_and_annotation(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    '''
    CHR based multiple task.
    From chr-seperated mapped bam to gtf;
    Annotate and filter chr-seperated gtf.
    :return:
    '''
    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/chr_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    sample_ids=[]
    excluded_samples=['AGTEX-CA241-5032-LN-YGE8']
    for f in os.listdir(f'{ROOT_OUTPUT_DIR}/sample_based'):
        if os.path.isdir(f'{ROOT_OUTPUT_DIR}/sample_based/{f}') and f not in excluded_samples:
            sample_ids.append(f)
    logger.info(f"load {len(sample_ids)} samples.")
    for i in df.index:
        series_cmds=[]
        chr_id=df.loc[i,'chr']
        root_out_dir=f'{OUTPUT_DIR}/{chr_id}'
        # merge mapped bam by chr.
        out_dir=f'{root_out_dir}/merged_bam'
        make_dir(out_dir)
        log_path=f'{out_dir}/merge_bam.log'
        merged_bam = f'{out_dir}/merge.mapped.bam'
        bam_files=[f'{ROOT_OUTPUT_DIR}/sample_based/{id}/{id}.chr_bam/{chr_id}.bam' for id in sample_ids]
        bam_files_path=f'{out_dir}/merge_bam_path.txt'
        with open(bam_files_path,'w') as bw:
            bw.write('\n'.join(bam_files))
        cmd=f'''
        {LOAD_SAMTOOLS_ENVS_CMD} && 
        samtools merge -f -@ {nt_per_task} -b {bam_files_path} {merged_bam} && 
        samtools index -@ {nt_per_task} {merged_bam}
        '''
        cmd = ' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # collapse by chr-based bam.
        out_dir=f'{root_out_dir}/collapse'
        make_dir(out_dir)
        log_path=f'{out_dir}/run_collapse.log'
        inner_log_path=f'{out_dir}/inner_collapse.log'
        out_gff=f'{out_dir}/merged.collapsed.gff'
        cmd=f'''
        {LOAD_ISOSEQ_ENVS_CMD} && 
        isoseq collapse {merged_bam} {out_gff}
        --num-threads {nt_per_task}  
        --log-file {inner_log_path}
        '''
        cmd = ' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def __merge_read_count_per_sample():
    read_count_path=f'{ROOT_OUTPUT_DIR}/flnc_count_merged/raw_isoform_count_matrix.csv'
    chr_based_dir=f'{ROOT_OUTPUT_DIR}/chr_based'
    out = read_count_path
    os.makedirs(os.path.dirname(out), exist_ok=True)
    chrs = [chr for chr in
            "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM others".split()]
    dfs=[]
    def add_chr2isoform(x:str,chr):
        arr=x.split('.')
        return '.'.join([arr[0],f'{chr}']+arr[1:])

    for chr in chrs:
        chr_idx = chrs.index(chr)+1
        log(f"[INFO] Chr: {chr}")
        flnc = f'{chr_based_dir}/{chr}/collapse/merged.collapsed.FLNC.count.txt'
        df=pl.read_csv(flnc)
        df = df.with_columns(
            pl.col("isoform").map_elements(
                lambda x: add_chr2isoform(x,chr_idx),
                return_dtype=pl.String
            )
        )
        dfs.append(df)
    log("[INFO] start combine")
    fdf = pl.concat(dfs, how="vertical")
    fdf.write_csv(out)
    log("[INFO] Done.")
    pass


def __merge_sqanti3_annotation():
    chr_dir=f'{ROOT_OUTPUT_DIR}/chr_based'
    out_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    make_dir(os.path.dirname(out_prefix))

    def insert_chr_into_pb(orig_id, chr_tag):
        if orig_id.startswith('PB'):
            return f'PB.{chr_tag}.{orig_id[3:]}'
        if orig_id.startswith('ENSG'):
            return orig_id
        else:
            return f'{orig_id}.{chr_tag}'

    def extract_chr_tag(chr_tag):
        chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM others".split()
        chr_id=chrs.index(chr_tag)+1
        return f'{chr_id}'

    def gtf_replace_pb_attrs(line, n):
        line = line.strip()
        attr_tx_re = re.compile(r'(transcript_id\s+")([^"]+)(";)')
        attr_gene_re = re.compile(r'(gene_id\s+")([^"]+)(";)')
        line = attr_tx_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), n)}{m.group(3)}',
                              line)
        line = attr_gene_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), n)}{m.group(3)}',
                                line)
        return line

    # Process per-chromosome directories
    dfs=[]
    chr_dir_tag={}
    for d in os.listdir(chr_dir):
        chr_tag = extract_chr_tag(d)
        chr_dir_tag[d]=chr_tag
    sorted_chrs=sorted(chr_dir_tag.keys(),key=lambda x: chr_dir_tag[x])
    # merge annotation
    for d in sorted_chrs:
        chr_tag = chr_dir_tag[d]
        print(f"Processing class: {d} with chr_tag='{chr_tag}'")
        # locate files
        class_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3_RulesFilter_result_classification.txt'
        if not os.path.exists(class_path):
            print(f'bad {d}')
            continue
        df = pd.read_csv(class_path, sep="\t", dtype=str, low_memory=False)
        tx_col = 'isoform'
        gene_col = 'associated_gene'
        # keep original columns for traceability
        df[f"orig_{tx_col}"] = df[tx_col].astype(str)
        df[f"orig_{gene_col}"] = df[gene_col].astype(str)
        df["source_chr"] = chr_tag
        df[tx_col] = df[tx_col].apply(lambda x: insert_chr_into_pb(x, chr_tag))
        df[gene_col] = df[gene_col].apply(lambda x: insert_chr_into_pb(x, chr_tag))
        dfs.append(df)
    mdf=pd.concat(dfs, ignore_index=True)
    mdf.to_csv(f'{out_prefix}.RulesFilter_result_classification.txt', sep="\t", index=False)
    log(f'save merged annotation.')

    # merge gtf
    gtf_output_lines=[]
    for d in sorted_chrs:
        chr_tag = extract_chr_tag(d)
        print(f"Processing gtf: {d} with chr_tag='{chr_tag}'")
        gtf_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3.filtered.gtf'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                # replace transcript_id if present and matches PB
                attrs = gtf_replace_pb_attrs(attrs, chr_tag)
                parts[8] = attrs
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.gtf', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged gtf.')

    # merge cds
    gtf_output_lines = []
    for d in sorted_chrs:
        chr_tag = extract_chr_tag(d)
        print(f"Processing gtf: {d} with chr_tag='{chr_tag}'")
        gtf_path = f'{chr_dir}/{d}/sqanti3/sqanti3_corrected.cds.gff3'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                # replace transcript_id if present and matches PB
                attrs = gtf_replace_pb_attrs(attrs, chr_tag)
                parts[8] = attrs
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.full.cds.gff3', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged CDS.')


def __custom_sqanti3_filter():
    suffix='.custom'
    min_read_a_sample = 5
    min_support_samples=2
    sqanti3_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    flnc_count_path=f'{ROOT_OUTPUT_DIR}/flnc_count_merged/raw_isoform_count_matrix.csv'
    df=pd.read_csv(flnc_count_path, index_col=0)
    sum_read=df.sum(axis=1)
    support_samples = (df >= min_read_a_sample).sum(axis=1)
    # 1. add supported all reads
    sq_df=pd.read_csv(f'{sqanti3_prefix}.RulesFilter_result_classification.txt',sep='\t',index_col=0)
    sq_df['support_read']=sum_read[sq_df.index]
    # 2. add supported samples (with read >= min_read_a_sample)
    sq_df['support_samples']=support_samples[sq_df.index]
    # save
    advanced_filter_path=f'{sqanti3_prefix}.advanced_filter.txt'
    sq_df.to_csv(advanced_filter_path,sep='\t')

    # update gtf and annotation file
    out_gtf=f'{sqanti3_prefix}{suffix}.gtf'
    out_cds_gtf = f'{sqanti3_prefix}{suffix}.cds.gff3'
    out_annot = f'{sqanti3_prefix}{suffix}.RulesFilter_result_classification.txt'
    df=pd.read_csv(f'{sqanti3_prefix}.advanced_filter.txt',sep='\t')
    df.loc[((df['structural_category']!='full-splice_match') & (df["polyA_motif_found"] == False)),'filter_result']='Artifact'
    df.loc[((df['structural_category']!='full-splice_match') & (df['support_samples']<min_support_samples)),'filter_result']='Artifact'
    isoform=df.loc[df['filter_result']=='Isoform',:].shape[0]
    log(f'isoform: {isoform}')
    # update annotation file
    df.to_csv(out_annot,sep='\t',index=False)
    log(f'save updated annotation file')
    transcripts=set(df.loc[df['filter_result']=='Isoform','isoform'].unique().tolist())
    log(f'isoform: {len(transcripts)}')
    # update gtf
    gtf_path=f'{sqanti3_prefix}.gtf'
    with open(gtf_path, 'r') as in_f, open(out_gtf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            tid_match = re.search(r'transcript_id\s+"([^"]+)"', parts[8])
            if tid_match and tid_match.group(1) in transcripts:
                out_f.write(line)
    log(f'save updated gtf')
    # update cds
    cds_gtf_path=f'{sqanti3_prefix}.full.cds.gff3'
    with open(cds_gtf_path, 'r') as in_f, open(out_cds_gtf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            tid_match = re.search(r'transcript_id\s+"([^"]+)"', parts[8])
            if tid_match and tid_match.group(1) in transcripts:
                out_f.write(line)
    log(f'save updated cds gff')

def __merge_collapsed_gtf():
    chr_dir=f'{ROOT_OUTPUT_DIR}/chr_based'
    out_prefix=f'{ROOT_OUTPUT_DIR}/collapsed_merged/isoseq'
    make_dir(os.path.dirname(out_prefix))

    def insert_chr_into_pb(orig_id, chr_tag):
        if orig_id.startswith('PB'):
            return f'PB.{chr_tag}.{orig_id[3:]}'
        if orig_id.startswith('ENSG'):
            return orig_id
        else:
            return f'{orig_id}.{chr_tag}'

    def extract_chr_tag(chr_tag):
        chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM others".split()
        chr_id=chrs.index(chr_tag)+1
        return f'{chr_id}'

    def gtf_replace_pb_attrs(line, n):
        line = line.strip()
        attr_tx_re = re.compile(r'(transcript_id\s+")([^"]+)(";)')
        attr_gene_re = re.compile(r'(gene_id\s+")([^"]+)(";)')
        line = attr_tx_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), n)}{m.group(3)}',
                              line)
        line = attr_gene_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), n)}{m.group(3)}',
                                line)
        return line

    # Process per-chromosome directories
    dfs=[]
    chr_dir_tag={}
    for d in os.listdir(chr_dir):
        chr_tag = extract_chr_tag(d)
        chr_dir_tag[d]=chr_tag
    sorted_chrs=sorted(chr_dir_tag.keys(),key=lambda x: chr_dir_tag[x])

    # merge gtf
    gtf_output_lines=[]
    for d in sorted_chrs:
        chr_tag = extract_chr_tag(d)
        print(f"Processing gtf: {d} with chr_tag='{chr_tag}'")
        gtf_path = f'{chr_dir}/{d}/collapse/merged.collapsed.gff'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                # replace transcript_id if present and matches PB
                attrs = gtf_replace_pb_attrs(attrs, chr_tag)
                parts[8] = attrs
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.gtf', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged gtf.')

def __merge_collapsed_gtf_fast():
    import os, re

    chr_dir = f'{ROOT_OUTPUT_DIR}/chr_based'
    out_prefix = f'{ROOT_OUTPUT_DIR}/collapsed_merged/isoseq'
    make_dir(os.path.dirname(out_prefix))
    chrs = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM others".split()
    chr2tag = {c: str(i + 1) for i, c in enumerate(chrs)}
    def insert_chr_into_pb(orig_id, chr_tag):
        if orig_id.startswith("PB"):
            return f"PB.{chr_tag}.{orig_id[3:]}"
        if orig_id.startswith("ENSG"):
            return orig_id
        return f"{orig_id}.{chr_tag}"
    tx_re = re.compile(r'(transcript_id\s+")([^"]+)(";)')
    gene_re = re.compile(r'(gene_id\s+")([^"]+)(";)')
    def replace_attrs(attr, chr_tag):
        attr = tx_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), chr_tag)}{m.group(3)}', attr)
        attr = gene_re.sub(lambda m: f'{m.group(1)}{insert_chr_into_pb(m.group(2), chr_tag)}{m.group(3)}', attr)
        return attr
    chr_dirs = sorted(
        os.listdir(chr_dir),
        key=lambda d: int(chr2tag.get(d, 999))
    )
    out_gtf = f"{out_prefix}.gtf"
    with open(out_gtf, "w") as fout:
        for d in chr_dirs:
            chr_tag = chr2tag.get(d)
            if chr_tag is None:
                continue
            print(f"Processing gtf: {d} -> chr_tag={chr_tag}")
            gtf_path = f"{chr_dir}/{d}/collapse/merged.collapsed.gff"

            with open(gtf_path) as fh:
                for ln in fh:
                    if not ln or ln[0] == "#":
                        continue

                    parts = ln.rstrip("\n").split("\t", 8)
                    if len(parts) < 9:
                        continue

                    parts[8] = replace_attrs(parts[8], chr_tag)
                    fout.write("\t".join(parts) + "\n")
    log("save merged gtf.")


def merge_result():
    # merge isoseq collapsed gff
    log(f'merging collapsed gtf')
    __merge_collapsed_gtf_fast()


if __name__ == '__main__':
    # sample based task
    if STEP == 'hifi_to_mapped_bam':
        hifi_to_mapped_bam(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    # chr based task
    if STEP == 'bam_to_isoform':
        collapse_and_annotation(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    # single task
    if STEP == 'merge':
        merge_result()



