# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/29 21:19
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import os.path
import re
import sys
import polars as pl

import pandas as pd

from util import run_commands_threadpool

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
LOG_NAME = ARGS[2]
N_TASK = int(ARGS[3])
NT_PER_TASK = int(ARGS[4])

LOAD_FLAIR_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate flair'




def flair_quant(conf_path, n_task, log_name, nt_per_task):
    '''
    Flair quant.
    :return:
    '''
    parallel_cmds=[]
    df=pd.read_csv(conf_path)
    for i,row in df.iterrows():
        fa_path=row['fa_path']
        fastq_path=row['fastq_path']
        out_dir=row['out_dir']
        log_path=f'{out_dir}/run_flair_quant.log'
        series_cmds=[]
        os.makedirs(out_dir, exist_ok=True)
        cmd=f'''
            {LOAD_FLAIR_ENVS_CMD} && 
            flair quantify -r {fastq_path} -i {fa_path} --threads {nt_per_task} --sample_id_only 
            --output {out_dir}/flair_quant 
        '''
        cmd=' '.join(re.split(r'\s+', cmd))
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def combine_flair():
    main_dir='/lustre/home/cxue/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'
    quant_dir=f'{main_dir}/flair_quant'
    ref_name='raw_gtf'
    flair_task_out_dir=f'{quant_dir}/output/{ref_name}'
    flair_combined_out_dir=f'{quant_dir}/combined'
    os.makedirs(flair_combined_out_dir, exist_ok=True)

    merged_df = None
    for task in os.listdir(flair_task_out_dir):
        print(f'start {task}')
        count_path=f'{flair_task_out_dir}/{task}/flair_quant.counts.tsv'
        if not os.path.exists(count_path):
            continue
        df = pl.read_csv(count_path,separator='\t',has_header=True)
        gene_col = df.columns[0]
        df = df.with_columns(pl.col(gene_col).cast(pl.String))
        df = df.unique(subset=[gene_col], keep="first")
        if merged_df is None:
            merged_df = df
        else:
            merged_gene_col = merged_df.columns[0]
            df_renamed = df.rename({gene_col: merged_gene_col})
            merged_df = merged_df.join(
                df_renamed,
                on=merged_gene_col,
                how="outer",
                coalesce=True
            )
    if merged_df is not None:
        merged_df = merged_df.fill_null(0)
        merged_df = merged_df.sort(merged_df.columns[0])
    merged_df.write_csv(f'{flair_combined_out_dir}/{ref_name}.transcript.count.flair.tsv',separator='\t')
    print(f'save {ref_name}: {merged_df.shape}')



def __count_to_tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    lengths = lengths.loc[counts.index]
    rpk = counts.div(lengths / 1000, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return tpm

def __isoform_to_gene_sum(df: pd.DataFrame, iso2gene: pd.Series) -> pd.DataFrame:
    df = df.copy()
    df["gene_id"] = iso2gene.loc[df.index].values
    sdf=df.groupby("gene_id").sum()
    return sdf


if __name__ == '__main__':
    if STEP == 'flair_quant':
        flair_quant(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'combine':
        combine_flair()