# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/11/23 21:43
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import os
import sys

import numpy as np
import pandas as pd

from util import get_LRS_correct_sample_id

RUN_LOG_DIR=f'/lustre/home/cxue/project/GMTiP-RNA/20260131/run_log'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
py_name=os.path.basename(__file__)[:-3]

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'
LOAD_POLARS_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate polars_env'


def _load_flair_quant_conf_df():
    main_dir='/lustre/home/cxue/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'
    quant_dir=f'{main_dir}/flair_quant'
    raw_gtf=f'{main_dir}/raw_isoform/tool2_filter.fa'
    fq_dir=f'/lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/flair/output/flair'
    excluded_samples=['AGTEX-CA241-5032-LN-YGE8']
    conf_dir=f'{quant_dir}/conf'
    output_dir=f'{quant_dir}/output'
    os.makedirs(conf_dir, exist_ok=True)
    n=0
    for f in os.listdir(fq_dir):
        if f.startswith('AGTEX') and os.path.isdir(os.path.join(fq_dir, f)) and f not in excluded_samples:
            fa=f'{fq_dir}/{f}/{f}.flnc.fastq'
            if not os.path.exists(fa):
                raise Exception(f'{fa} does not exist')
            sample_id=get_LRS_correct_sample_id(f)
            data=[[sample_id,'cond1','batch1',fa]]
            df=pd.DataFrame(data)
            df.to_csv(f'{conf_dir}/{sample_id}.flair_quant_fa.conf.txt', index=False, sep='\t', header=False)
            n+=1
    print(f'find {n} samples fq')
    # all task
    ref_fas={'raw_gtf':raw_gtf}
    data = []
    for ref_name,ref_fa in ref_fas.items():
        for f in os.listdir(conf_dir):
            sid=f.split('.')[0]
            data.append([f'{conf_dir}/{f}',ref_fa,f'{output_dir}/{ref_name}/{sid}'])
    df=pd.DataFrame(data,columns=['fastq_path','fa_path','out_dir'])
    return df

def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        py_name='flair_quant_pipeline.py',
        conf_sep=',',has_header=True
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    # assign task to nodes, build a sample list file as para for pipeline script per node.
    if df is not None:
        sub_dfs = np.array_split(df, NODES)
    else:
        sub_dfs=[None,]
    cmds = []
    for i, sdf in enumerate(sub_dfs, 1):
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{py_name}_{STEP}/node{i}'
        os.makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.txt'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False, sep= conf_sep, header=has_header)
        # build hpc job submit list
        shell_lines = [
            f'#!/bin/bash',
            f'#SBATCH -J {STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            LOAD_POLARS_ENVS_CMD,
            f'python {CURRENT_DIR}/{py_name} {STEP} {conf_path} {node_log_path} '
            f'{N_TASK_PER_NODE} {NT_PER_TASK}'
        ]
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmd = f'''
            sbatch {hpc_job_shell}
        '''
        cmd = ' '.join(cmd.split())
        cmds.append(cmd)
    for i, cmd in enumerate(cmds, 1):
        print(cmd)
        os.system(cmd)



if __name__ == '__main__':
    COMPUTER='HPC'
    ARGS=sys.argv[1:]
    STEP=ARGS[0]
    ## for step 1: sample-based multiple task.
    if STEP == 'flair_quant':
        df=_load_flair_quant_conf_df()
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 176
        CPU_PER_NODE = 10
        NT_PER_TASK = 12
        N_TASK_PER_NODE = 1
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    ## for step 2: combine
    if STEP == 'combine':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 1
        CPU_PER_NODE = 10
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)
