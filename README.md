# Hashimoto-Thyroiditis-RNAseq-Pipeline
An in silico RNA-seq pipeline for differential gene expression analysis, functional enrichment analysis and systems biology modeling in Hashimoto's Thyroiditis (HT)

The pipeline transitions from raw sequencing data to systems-level analysis. The initial preprocessing phase was executed as follows:
Quality Control: Evaluated raw reads using `FastQC` and compiled global reports with `MultiQC`.
Adapter Trimming: Performed quality filtering, poly(G) tail and adapter removal using `fastp`.
Quantification: Conducted transcript-level pseudo-alignment and quantification using `Salmon`

The raw FASTQ sequencing data utilized in this project were originally generated and published by Zeng 
et al. (2025). 

If you use or reference the underlying biological data, please cite the original study: Zeng H, Chen Y, Hu T, et al. *Integrative Transcriptome and Machine Learning Analysis Uncovers Critical STAT3/GREM2 Signaling Mechanisms in Dexamethasone Treatment of Hashimoto's Thyroiditis.* ACS Omega. 2025;10(46):55377-55391.
Published 2025 Nov 13.
[10.1021/acsomega.5c04845](https://doi.org/10.1021/acsomega.5c04845)


THIS REPOSITORY IS STILL UNDER DEVELOPMENT
