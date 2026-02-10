# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/5 09:24
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
import subprocess
import sys
from os import makedirs

import numpy as np
import pandas as pd

# config for computer evn


RUN_LOG_DIR=f'/lustre/home/cxue/project/GMTiP-RNA/20260131/run_log'
OUTPUT_DIR=f'/lustre/home/cxue/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'

def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        pipeline_path = 'enhanced_gtf_pipeline.py',
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
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{STEP}/node{i}'
        makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.csv'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False)
        # build hpc job submit list
        shell_lines = [
            f'#!/bin/bash',
            f'#SBATCH -J {STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            f'python {CURRENT_DIR}/{pipeline_path} {STEP} {conf_path} {OUTPUT_DIR} {node_log_path} '
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
    COMPUTER=ARGS[0]
    STEP=ARGS[1]
    ## single task
    if STEP == 'enhanced_gtf':
        PARTITION='cu-debug,cpuPartition,cu-1'
        NODES = 1
        CPU_PER_NODE = 1
        NT_PER_TASK = 20
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)


