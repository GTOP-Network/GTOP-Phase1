#!/bin/bash

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/flair/output/flair/concat"

prefix="flair"
outdir="${workdir}/all"
mkdir -p $outdir

function merge_gtfs() {
    local chrs=($(seq 1 22 | sed 's/^/chr/') chrX chrY chrM)
    touch ${outdir}/${prefix}.isoforms.gtf
    for chr in "${chrs[@]}"; do
        sed -n '1,$p' ${workdir}/${chr}/${chr}.isoforms.gtf >> ${outdir}/${prefix}.isoforms.gtf
    done
}

# main func
function main() {
    merge_gtfs
}
main