# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/5 22:30
@Desc    : Description
"""
import os
import threading
import subprocess
from datetime import datetime
from tqdm import tqdm
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    os.makedirs(os.path.dirname(main_log),exist_ok=True)
    log(f"[INFO] Batch start. Total tasks: {len(cmd_list)}, Max workers: {max_workers}", main_log)
    def worker(task_cmds):
        for cmd_it in task_cmds:
            cmd=cmd_it['cmd']
            cmd = re.sub(r'\s+', ' ', cmd)
            task_log=cmd_it['log_path']
            logs_dir=os.path.dirname(task_log)
            os.makedirs(logs_dir, exist_ok=True)
            main_log_path=f'{logs_dir}/{main_log_prefix}.log'
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


def get_LRS_correct_sample_id(sample_id):
    ind_codes = {'CC331': 'CC311'}
    tissue_merge = {'0260': '0262','FZDM':'4151','QZDM':'4150'}
    arr=sample_id.split('-')
    ind_id=arr[1]
    if ind_id in ind_codes:
        ind_id=ind_codes[ind_id]
        arr[1]=ind_id
    tissue_id=arr[2]
    if tissue_id in tissue_merge:
        tissue_id=tissue_merge[tissue_id]
        arr[2]=tissue_id
    return '-'.join(arr)

