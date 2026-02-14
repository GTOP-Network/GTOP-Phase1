## QTL and GWAS colocalization
This directory contains a series of R scripts designed for integrative fine-mapping of genetic associations, combining results from Genome-Wide Association Studies (GWAS) and QTL. The pipeline supports causal variant prioritization, gene-based annotation, colocalization analysis, and result integration across multiple datasets.

## Script Overview
- a1.finemapping_for_GWAS.R
Performs fine-mapping on GWAS summary statistics to identify likely causal variants within associated loci. Uses Bayesian or likelihood-based methods (e.g., FINEMAP, SuSiE) to estimate posterior probabilities of causality.
- a1.GWAS_loci.R
Identifies and annotates significant GWAS loci based on p-values, LD structure, and genomic features. Outputs locus-level summaries for downstream analysis.
- a2.finemapping_prepare_genes.R
Prepares gene-based annotations and regulatory information (e.g., eQTLs, promoter regions, enhancers) to support gene-centric fine-mapping and interpretation.
- a3.finemapping_for_QTL.R
Conducts fine-mapping specifically for QTL signals (e.g., eQTLs, sQTLs), identifying the most probable causal SNPs influencing gene expression or splicing.
- a4.coloc.R
Performs colocalization analysis between GWAS and QTL signals using tools like coloc or MASHR, testing whether the same causal variant underlies both trait association and gene regulation.
- a5.merge_result.R
Integrates results from all previous steps (GWAS fine-mapping, QTL fine-mapping, colocalization) into a unified output table, enabling comprehensive interpretation of shared and distinct causal mechanisms.

## Workflow Summary
- GWAS Fine-Mapping → a1.finemapping_for_GWAS.R
- GWAS Loci Annotation → a1.GWAS_loci.R
- Gene Annotation Preparation → a2.finemapping_prepare_genes.R
- QTL Fine-Mapping → a3.finemapping_for_QTL.R
- Colocalization Analysis → a4.coloc.R
- Result Integration → a5.merge_result.R
