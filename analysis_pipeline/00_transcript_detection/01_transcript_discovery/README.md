# Transcript detection using three complementary approaches

To comprehensively identify transcript models from long-read RNA sequencing data, we implemented three complementary transcript discovery pipelines. Two reference-guided approaches (Bambu and FLAIR) were used to leverage existing genome annotations, while a reference-free Iso-Seq pipeline enabled unbiased transcript reconstruction. All tools were executed with recommended parameters unless otherwise specified.

---

## Transcript detection with Bambu (v3.12.0)

Bambu was applied in a reference-guided manner to identify novel transcript structures while preserving annotated isoforms. Intermediate annotation and read-class objects were reused to improve computational efficiency and reproducibility.

```bash
cd bambu
```

### Step 1: Preparation of annotation objects

```bash
qsub run.step1.sh
```

Genome annotations were preprocessed using the `prepareAnnotations()` function to generate `bambuAnnotations` objects. These objects were cached and reused across multiple runs to avoid repeated preprocessing.

### Step 2: Generation and reuse of read class files (rcFiles)

```bash
qsub run.step2.sh
```

Read classes were constructed and stored as read class files (`rcFiles`). These files were reused in subsequent Bambu runs to reduce computational overhead when performing transcript discovery under different configurations or sample combinations.

### Step 3: Transcript discovery without quantification

```bash
qsub run.step3.sh
```

Transcript discovery was performed using aligned reads (BAM files), reference annotations, and the reference genome sequence. Expression quantification was disabled (`quant = FALSE`) to focus exclusively on transcript structure identification. Bambu was run with `NDR = 1`, and downstream filtering was applied to obtain high-confidence transcript models. The output is a `GRangesList` object containing both annotated and novel transcripts.

---

## Transcript detection with FLAIR (v2.2.0)

FLAIR was used as an independent reference-guided method with explicit splice-site correction based on genome annotations.

```bash
cd flair
```

### Step 1: Read alignment and splice-site correction

```bash
qsub run.step1.sh
```

Raw reads in FASTA or FASTQ format were aligned to the reference genome. SAM files were converted to BED12 format, followed by splice-site correction using genome annotations. Corrected reads were split by chromosome for parallel processing.

### Step 2: Isoform assembly

```bash
qsub run.step2.sh
```

Corrected reads were collapsed into high-confidence isoforms using recommended FLAIR parameters optimized for both known and novel transcript detection.

### Step 3: Merging chromosome-level annotations

```bash
qsub run.step3.sh
```

Chromosome-specific GTF files were merged to generate a unified transcript annotation.

---

## Transcript detection with Iso-Seq (v4.3.0)

A reference-free Iso-Seq pipeline was implemented to enable unbiased transcript discovery. Processing was parallelized at both the sample and chromosome levels to improve scalability.

```bash
cd isoseq
```

### Step 1: Sample-level processing

```bash
python isoseq_run.py hifi_to_mapped_bam
```

HiFi reads were clustered at the sample level, aligned to the reference genome, and split by chromosome to generate chromosome-specific BAM files for each sample.

### Step 2: Chromosome-level isoform collapsing

```bash
python isoseq_run.py bam_to_isoform
```

For each chromosome, BAM files from all samples were merged and collapsed to generate non-redundant isoform models.

### Step 3: Merging chromosome-level annotations

```bash
python isoseq_run.py merge
```

All chromosome-specific GTF files were merged into a final genome-wide Iso-Seq transcript annotation.
