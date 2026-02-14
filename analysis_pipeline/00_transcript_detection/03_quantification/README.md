# Long-read RNA-seq transcript and gene quantification (read count and TPM)

Long-read RNA-seq quantification was performed using two transcript reference annotations.  
First, quantification was conducted using the long-read–derived reference generated in this study, and these results were used for the primary long-read–based analyses.  
Second, quantification was performed using the GENCODE v47 annotation to enable direct comparison with short-read RNA-seq quantification, ensuring consistency in the reference annotation across sequencing platforms.

---

```bash
cd scripts
```

To improve computational efficiency, transcript quantification was performed in parallel at the sample level, followed by aggregation across samples.

## 1. Sample-level quantification using FLAIR

```bash
python flair_quant_run.py flair_quant
```

## 2. Merging quantification results and calculating gene-level quantification

Transcript-level read counts from all samples were first merged. Transcript TPM values were then calculated according to the standard TPM definition. Gene-level TPM and read counts were obtained by summing the TPM values or read counts of all transcripts belonging to each gene.

```bash
python flair_quant_run.py combine
```

---
