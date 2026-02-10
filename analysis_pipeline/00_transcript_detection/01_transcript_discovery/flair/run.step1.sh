#!/bin/bash
# flnc(.bam) -> fastq -> flair-align -> flair-correct -> flair-collapse

#SBATCH --job-name=flair
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --array=1-176
#SBATCH --output=tmp/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/flair"
refg="${workdir}/ref/genome.fa"
refgtf="${workdir}/ref/gencode.v47.annotation.gtf"

sampleids="${workdir}/input/sampleids.txt"
sampleid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sampleids")
indir="/lustre/home/cxue/project/GMTiP-RNA/20251031/long_read/HPC/output/sample_based/${sampleid}"
outdir="${workdir}/output/flair/${sampleid}"
chrdir="${indir}/chrs"
mkdir -p $chrdir

# generate fastq
function exec_bam2fq() {
    # iso-seq flnc bam -> fastq
    echo "#### start flnc.bam-to-fastq at $(date "+%Y-%m-%d %H:%M:%S") ####"
    samtools fastq --threads ${THREADS} ${indir}/${sampleid}.flnc.bam > ${outdir}/${sampleid}.flnc.fastq && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# align with flair-align
function exec_flair_align() {
    echo "#### start flair-align at $(date "+%Y-%m-%d %H:%M:%S") ####"
    flair align -g ${refg} -r ${outdir}/${sampleid}.flnc.fastq -o ${outdir}/${sampleid} -t ${THREADS} && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# correct with flair-correct
function exec_correct() {
    echo "#### start flair-correct at $(date "+%Y-%m-%d %H:%M:%S") ####"
    flair correct -q ${outdir}/${sampleid}.bed -g ${refg} -f ${refgtf} \
    -o ${outdir}/${sampleid} -t ${THREADS} --ss_window 5 && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# split bam by chr then convert bam to fastq.
function split_bam_by_normal_chr() {
    echo "#### start bam-spliting at $(date "+%Y-%m-%d %H:%M:%S") ####"
    local bam="${outdir}/${sampleid}.bam"
    local chrs=($(seq 1 22 | sed 's/^/chr/') chrX chrY chrM) && \
    for chr in "${chrs[@]}"; do
        samtools view -b -@ ${THREADS} "${bam}" "${chr}" | \
        samtools fastq -@ ${THREADS} - > "${chrdir}/${sampleid}.${chr}.flnc.fastq"
    done
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# split bed by chr
function split_bed_by_chr() {
    echo "#### start bed-spliting at $(date "+%Y-%m-%d %H:%M:%S") ####"
    local bed="${outdir}/${sampleid}_all_corrected.bed"
    cd $chrdir && \
    awk '{print > $1".bed"}' $bed && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}

# workflow
function run_flair() {
    echo "#### start flnc-to-flair-workflow at $(date "+%Y-%m-%d %H:%M:%S")"
    echo '#### convert bam to fastq' && \
    exec_bam2fq && \
    echo '#### flair align' && \
    exec_flair_align && \
    echo '#### flair correct' && \
    exec_correct && \
    echo '#### bam split' && \
    split_bam_by_normal_chr && \
    echo '#### bed split' && \
    split_bed_by_chr && \
    echo "#### end workflow at $(date "+%Y-%m-%d %H:%M:%S")"
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate longrw_flair
    run_flair
}
main
