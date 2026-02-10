#!/bin/bash
# flair-correct(s) -> concat -> flair-collapse

#SBATCH --job-name=flr
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --array=1-25
#SBATCH --output=tmp/flr/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/flr/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/flair"

chrs="${workdir}/input/chrs.list"
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chrs")
refg="${workdir}/ref2/chrs/${chr}.fa"
refgtf="${workdir}/ref2/chrs/${chr}.gtf"
indir="/lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/flair/output/flair"
outdate=$(date "+%Y-%m-%d")
outdir="${workdir}/output/flair/concat/${outdate}/${chr}"
mkdir -p $outdir

# merge the bed files after running FLAIR correct
# If using PacBio reads, be careful doing this, as the reads may not have unique IDs.
# You may need to label each read with its sample ID to keep the read IDs unique.
function exec_concat() {
    echo "#### start merge-beds at $(date "+%Y-%m-%d %H:%M:%S") ####"
    cat ${indir}/*/chrs/${chr}.bed > ${outdir}/${chr}_all_corrected.bed && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# collapse by flair-collapse with recommended parameters for human.
function exec_collapse_human() {
    echo "#### start flair-collapse(human) at $(date "+%Y-%m-%d %H:%M:%S") ####"
    flair collapse -g ${refg} -r ${indir}/*/chrs/AGTEX-*.${chr}.flnc.fastq \
    -q ${outdir}/${chr}_all_corrected.bed \
    -o ${outdir}/${chr} -t ${THREADS} --gtf ${refgtf} \
    --stringent --check_splice --generate_map --annotation_reliant generate && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}

# workflow after merging bed files
function run_collapse_after_merging_beds() {
    echo "#### start collapse after merging beds at $(date "+%Y-%m-%d %H:%M:%S")"
    echo '#### merge bed files'
    exec_concat && \
    echo '#### flair collapse' && \
    exec_collapse_human && \
    echo "#### end workflow at $(date "+%Y-%m-%d %H:%M:%S")"
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate longrw_flair
    run_collapse_after_merging_beds
}
main
