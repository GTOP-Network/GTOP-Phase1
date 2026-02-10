# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/24 09:16
@Email   : xuechao@szbl.ac.cn
@Desc    :
"""

# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/30 11:07
@Desc    : Run Iso-Seq pipeline in HPC or Single node.
"""
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import os
import subprocess
import sys
from statistics import median

import numpy as np
import pandas as pd
import polars as pl

# config for computer evn

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))



REF_DIR='/lustre/home/cxue/raw_data/GMTiP-LR-RNA-Seq/ref'
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
ROOT_OUTPUT_DIR = ARGS[2]
LOG_NAME = ARGS[3]
N_TASK = int(ARGS[4])
NT_PER_TASK = int(ARGS[5])

LOAD_ISOSEQ_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate isoseq'
LOAD_PBINDEX_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate pbtk'
LOAD_SAMTOOLS_ENVS_CMD='module load samtools'

LOAD_SQANTI3_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate sqanti3'
SQANTI3_DIR = f"/lustre/home/cxue/software/sqanti3/release_sqanti3"

# LOAD_SQANTI3_ENVS_CMD='source ~/.bashrc && mamba activate SQANTI3.env'
# SQANTI3_DIR = '/lustre/home/cxue/software/sqanti3_5.2.1'

LOAD_isoLASER_ENVS_CMD='module load anaconda && source ~/.bashrc && conda activate isoLASER'
LOAD_PICARD_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate picard_env'
LOAD_SALMON_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate salmon_env'
LOAD_POLARS_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate polars_env'

helper_py=f'{LOAD_POLARS_ENVS_CMD} && python {CURRENT_DIR}/alter_tool_helper.py'


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



def sqanti3_annot_filter(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    '''
    :return:
    '''

    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/chr_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    SQANTI3_QC = f"{SQANTI3_DIR}/sqanti3_qc.py"
    SQANTI3_FILTER = f"{SQANTI3_DIR}/sqanti3_filter.py"
    for i in df.index:
        series_cmds=[]
        chr_id=df.loc[i,'chr']
        root_out_dir=f'{OUTPUT_DIR}/{chr_id}'
        gtf_path=f'{root_out_dir}/gtf/merged.gtf'
        # sqanti3 qc and filter.
        log(f'start sqanti3')
        out_dir = f'{root_out_dir}/sqanti3'
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)
        filter_out_dir = f'{out_dir}/filter'
        make_dir(out_dir,filter_out_dir)
        log_path = f'{out_dir}/run_sqanti3.log'
        out_prefix = 'sqanti3'

        chunk_size=8
        if chr_id in ['chrM','other']:
            chunk_size=1
        cmd = f'''
                {LOAD_SQANTI3_ENVS_CMD} && 
                python {SQANTI3_QC}
                --isoforms {gtf_path}
                --refGTF {REF_GTF}
                --refFasta {REF_GENOME_FASTA}
                --CAGE_peak {CAGE_PEAK}
                --polyA_motif_list {POLYA}
                --saturation
                --report skip
                --output {out_prefix}
                --dir {out_dir}
                --cpus {nt_per_task}
                --chunks {chunk_size}
                &&
                python {SQANTI3_FILTER} rules 
                --sqanti_class {out_dir}/{out_prefix}_classification.txt 
                --filter_gtf {out_dir}/{out_prefix}_corrected.gtf
                --filter_faa {out_dir}/{out_prefix}_corrected.faa
                --skip_report
                --filter_mono_exonic
                --output {out_prefix} 
                --dir {filter_out_dir}
                --cpus {nt_per_task}
                '''
        cmd=' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def merge_sqanti3_annotation():
    chr_dir=f'{ROOT_OUTPUT_DIR}/chr_based'
    out_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    make_dir(os.path.dirname(out_prefix))
    # Process per-chromosome directories
    dfs=[]
    sorted_chrs=sorted([x for x in os.listdir(chr_dir) if x.startswith('chr')])

    # merge annotation
    for d in sorted_chrs:
        print(f"Processing class: {d}")
        # locate files
        class_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3_RulesFilter_result_classification.txt'
        if not os.path.exists(class_path):
            print(f'bad {d}')
            continue
        df = pd.read_csv(class_path, sep="\t", dtype=str, low_memory=False)
        dfs.append(df)
    mdf=pd.concat(dfs, ignore_index=True)
    mdf.to_csv(f'{out_prefix}.RulesFilter_result_classification.txt', sep="\t", index=False)
    log(f'save merged annotation.')

    # merge gtf
    gtf_output_lines=[]
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3.filtered.gtf'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.gtf', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged gtf.')

    # merge cds
    gtf_output_lines = []
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/sqanti3_corrected.cds.gff3'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.full.cds.gff3', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged CDS.')

    # merge predicted protein sequencing
    gtf_output_lines = []
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/sqanti3_corrected.faa'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                gtf_output_lines.append(ln.strip()+'\n')
    with open(f'{out_prefix}.faa', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged protein sequence.')


def custom_sqanti3_filter():
    min_read_a_sample = 5
    min_support_samples = 2
    sqanti3_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    suffix='.custom'
    flnc_count_path=f'{ROOT_OUTPUT_DIR}/flair_quant/combined/raw_gtf.transcript.count.flair.tsv'

    # 1. add supported all reads
    sq_df=pd.read_csv(f'{sqanti3_prefix}.RulesFilter_result_classification.txt',sep='\t',index_col=0)
    df = pd.read_csv(flnc_count_path, index_col=0, sep='\t')
    df = df.reindex(df.index.union(sq_df.index), fill_value=0)
    sum_read=df.sum(axis=1)
    print(df)
    print(sum_read)
    print(sq_df)
    support_samples = (df >= min_read_a_sample).sum(axis=1)
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
    # all isoform with >= 10 supported reads
    df.loc[(df['support_read']<10),'filter_result']='Artifact'
    # non-FSM
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


if __name__ == '__main__':
    if STEP == 'run_sqanti3':
        sqanti3_annot_filter(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'merge':
        merge_sqanti3_annotation()
    if STEP == 'custom_filter':
        custom_sqanti3_filter()