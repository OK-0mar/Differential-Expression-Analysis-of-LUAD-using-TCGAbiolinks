# Differential-Expression-Analysis-of-LUAD-using-TCGAbiolinks
The primary goal of this project was to replicate the integrative analytical workflow of the TCGAbiolinks R/Bioconductor package on a novel dataset

Methodology & Computational Pipeline:

Using R, I established an end-to-end bioinformatics pipeline encompassing data retrieval, preprocessing, and statistical analysis:

1. Data Retrieval: Queried the Genomic Data Commons (GDC) to download TCGA-LUAD transcriptome profiling data (STAR - raw counts) for both Primary Tumor and Solid Tissue Normal samples.

2. Data Optimization: To ensure efficient pipeline testing and manage memory, I downsized the working dataset to 10 tumor and 10 normal samples before compiling them into a SummarizedExperiment object.

3. Preprocessing: Applied quantile filtering to remove the lowest 25% of expressed genes, reducing statistical noise and improving the power of the subsequent analysis.

4. Differential Expression Analysis (DEA): Utilized the edgeR exactTest pipeline built into TCGAbiolinks to compare the tumor and normal matrices. Significant differentially expressed genes (DEGs) were defined by an FDR threshold of < 0.05 and a Log2 Fold-Change threshold of 1.

5. Enrichment Analysis (EA): Mapped Ensembl IDs to standard HUGO gene symbols to prepare clean gene lists for pathway analysis. I then executed GO (Biological Process, Cellular Component, Molecular Function) and KEGG pathway enrichment analyses on the upregulated and downregulated gene clusters.

Key Visualizations & Outputs:

The pipeline successfully generated robust visual interpretations of the LUAD dataset, including:

1. Volcano Plots: Visualizing the significance state and highlighting the top 20 most significantly altered genes.

2. Enrichment Barplots: Multi-panel PDF charts detailing the top 10 most significantly enriched biological pathways for both the up-regulated and down-regulated gene sets.
