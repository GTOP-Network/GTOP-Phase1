#!/bin/bash
#SBATCH -J merge_gtf_job
#SBATCH -o /lustre/home/cxue/project/GMTiP-RNA/20260131/merge_gtf.o
#SBATCH -e /lustre/home/cxue/project/GMTiP-RNA/20260131/merge_gtf.e
#SBATCH -p cpuPartition
#SBATCH -N 1
#SBATCH -n 10

source ~/.bashrc && mamba activate polars_env

CURRENT_DIR="/lustre/home/cxue/src/python3/GTOP-RNA/public/analysis_pipeline/LR_RNA/merge_transcript_source/script"
python "$CURRENT_DIR/merge_pipeline.py"
