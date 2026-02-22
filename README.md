Human COX1 Gene Analysis
This repository contains a bioinformatics workflow for analyzing the Human Cytochrome c Oxidase Subunit I (MT-CO1) mitochondrial gene using R.

Project Overview
The COX1 gene plays a critical role in the mitochondrial electron transport chain.
This project performs computational analysis of the human COX1 sequence to examine its structural and functional characteristics.

Analyses Performed
Sequence loading and preprocessing from FASTA file
Basic sequence statistics (length, overview)
GC content calculation
Nucleotide frequency analysis (A, T, G, C)
Motif identification (Start/Stop codons, CpG sites)
ORF detection and protein translation
Visualization of nucleotide and GC distribution

📊 Results
The nucleotide distribution analysis shows a higher frequency of A, T, and G compared to C in the mitochondrial COX1 sequence.
GC content percentage was calculated to assess sequence stability.
Plots generated:
GC Content Distribution
Nucleotide Frequency Bar Plot

Tools & Libraries Used
R Programming
Biostrings (Bioconductor)
Base R plotting

Files Included
human_cox1_analysis.R
COX1_Human.fasta
Generated PNG plots

README.md
LICENSE

How to Run
Install Biostrings package: Copy code

BiocManager::install("Biostrings")
Set working directory to project folder.
Run human_cox1_analysis.R