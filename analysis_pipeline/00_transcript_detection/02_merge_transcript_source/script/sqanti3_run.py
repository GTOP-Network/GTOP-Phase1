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


RUN_LOG_DIR=f'/lustre/home/cxue/project/GMTiP-RNA/20260131/run_log'
OUTPUT_DIR='/lustre/home/cxue/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'
LOAD_POLARS_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate polars_env'


def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        py_name='sqanti3_pipeline.py',
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
            f'python {CURRENT_DIR}/{py_name} {STEP} {conf_path} {OUTPUT_DIR} {node_log_path} '
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
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]
    ## for step 1: chr-based multiple task.
    if STEP == 'run_sqanti3':
        chr_df = pd.DataFrame([f'chr{x}' for x in [i for i in range(1,23)]+['X','Y','M']],
                              columns=['chr'])
        PARTITION='cpuPartition,cu-1'
        NODES = 26
        CPU_PER_NODE = 32
        NT_PER_TASK = 32
        N_TASK_PER_NODE = 1
        submit_job(chr_df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    if STEP == 'merge':
        PARTITION='cpuPartition,cu-1'
        NODES = 1
        CPU_PER_NODE = 1
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    if STEP == 'custom_filter':
        PARTITION='cpuPartition,cu-1'
        NODES = 1
        CPU_PER_NODE = 32
        NT_PER_TASK = 32
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)