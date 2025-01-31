# Differential expression analysis from bulk RNA-seq data
This project aims to determine the effect of toxoplasma infection on gene expression in mouse blood and lung. The goal is to produce lists of genes that are differentially expressed between two experimental groups, and identify gene ontology (GO) terms enriched for DE genes.

The dataset includes 3–5 biological replicates from two distinct mouse tissues, encompassing both Toxoplasma-infected and uninfected control samples, derived from the study by Singhania et al., 2019. Libraries were prepared using a strand-specific protocol and sequenced in high-resolution paired-end mode with Illumina HiSeq 4000.

|  | Blood | Lung tissue |
| --- | --- | --- |
| Uninfected mouse | 3 | 3 |
| Infected mouse | 5 | 5 |

# Overview of the workflow
This workflow is based on the guidelines given in the context of the RNA sequencing course of the University of Bern (2024).
## 1 Quality checks
Scripts 01 and 02  
This is the first step of the workflow. Quality results may determine the following steps of the downstream analysis.
## 2 Mapping to the reference genome
Scripts 03, 04, 05, 06, 07, 08 and 12  
Aligning of the reads to the reference genome of *Mus Musculus* and then ensuring for the quality of the mapping.
## 3 Production of a table of counts with the number of reads per gene
Scripts 09 and 11  
Use of the files previously generated to count the number of reads per gene in each sample.
## 4 Differential expression analysis and overrepresentation analysis
Script 10  
After the exploratory data analysis, pairwise contrasts of infected and uninfected within lungs samples were made in order to identify differentially expressed genes. An overrepresentation analysis was also performed to identify more differentially expressed genes than expected.
