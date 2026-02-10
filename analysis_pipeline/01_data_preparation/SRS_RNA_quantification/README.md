# GTOP SRS RNA Expression Quantification
***
This directory contains code to replicate the gene expression quantification procedure done as part of the GTOP publication. This procedure is broadly split into three sections:
1. Quality control with FastQC.
2. Read mapping with STAR.
3. Gene-level quantification with RNA-SeQC.

### Quality control with FastQC
***
1. rawdir: Raw sequencing data directory (.fq.gz).
2. baseDir: Project base directory.
3. outDir: Quality control output directory.
4. sample_list: a file listing sample names.

### Read mapping with STAR
***
The `run_mapping` function performs alignment of paired-end RNA-seq reads to the reference genome using STAR (v2.7.3a) and generates sorted BAM files for downstream analysis.
1. dir: Project base directory.
2. start_index: STAR genome Index directory.
3. ref_gtf: Gene annotation file.
4. BAM_dir: bam output directory.
5. `task_use.txt`: a text file listing sample names.

### Gene-level quantification with RNA-SeQC.
***
The `run_rnaseqc` function performs RNA-SeQC quantification for aligned RNA-seq BAM files.
1. dir: tools and Gene annotation file directory.
2. bamDir: bam output directory.
3. outDir: quantification output directory.
4. curDir: Project base directory.
5. `task_use.txt`: a text file listing sample names.