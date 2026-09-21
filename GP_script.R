########################################
# DAV-R course graduation project.
# Title: How to retrieve data from TCGA and perform Differential Expression Analysis
# and Pathway Enrichment Analysis on LUAD data.
# R version: 4.6.1
# Script by: OK
########################################


# 1. INSTALL AND LOAD REQUIRED PACKAGES

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks", force = TRUE)
BiocManager::install("SummarizedExperiment", force = TRUE)
BiocManager::install("edgeR", force = TRUE)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

# 2. QUERY THE GENOMIC DATA COMMONS (GDC)

# Requesting STAR count RNA-seq data for both LUAD Tumors and Normal tissues
query_LUAD <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# Just a little modification to downsize the samples number

query_LUAD$results[[1]] <- rbind(
   head(query_LUAD$results[[1]][query_LUAD$results[[1]]$sample_type == "Primary Tumor", ], 10),
   head(query_LUAD$results[[1]][query_LUAD$results[[1]]$sample_type == "Solid Tissue Normal", ], 10)
 )


# 3. DOWNLOAD AND PREPARE THE DATA

GDCdownload(query_LUAD)
data_LUAD <- GDCprepare(query_LUAD)

# Extract the unstranded raw count matrix from the SummarizedExperiment object
counts_matrix <- assay(data_LUAD, "unstranded")


# 4. PRE-PROCESSING & FILTERING

# Remove the lowest 25% of expressed genes to reduce noise and improve 
# the statistical power of the DEA.
dataFilt <- TCGAanalyze_Filtering(
  tabDF = counts_matrix,
  method = "quantile",
  qnt.cut = 0.25 
)


# 5. DIFFERENTIAL EXPRESSION ANALYSIS (DEA)

# Separate the column names into Normal and Tumor groups
samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = "NT") # Solid Tissue Normal
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = "TP") # Primary Tumor

# Run the DEA using the edgeR exactTest pipeline built into TCGAbiolinks
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNT],
  mat2 = dataFilt[, samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  pipeline = "edgeR",
  method = "exactTest",
  fdr.cut = 0.05,     # False Discovery Rate threshold
  logFC.cut = 1       # Log2 Fold-Change threshold
)


# 6. VISUALIZE RESULTS: VOLCANO PLOT

# Extract the top 20 most significant genes (lowest FDR) to highlight
top_genes <- rownames(dataDEGs[order(dataDEGs$FDR), ])[1:20]


# Generate and save a volcano plot highlighting the top differentially expressed genes
TCGAVisualize_volcano(
  x = dataDEGs$logFC,
  y = dataDEGs$FDR,
  filename = "TCGA_LUAD_Volcano_Plot.pdf",
  x.cut = 1,
  y.cut = 0.05,
  names = rownames(dataDEGs),
  show.names = "highlighted", 
  highlight = top_genes,      
  title = "Differential Expression: TCGA-LUAD (Tumor vs. Normal)",
  legend = "Significance State"
)

# 7. PREPARE GENE LISTS (UP AND DOWN REGULATED)

# Filter DEGs into upregulated and downregulated groups
up_genes <- subset(dataDEGs, logFC > 1 & FDR < 0.05)
down_genes <- subset(dataDEGs, logFC < -1 & FDR < 0.05)

# TCGA raw counts use Ensembl IDs (e.g., ENSG0000...), but pathway enrichment 
# databases require standard Gene Symbols (e.g., TP53, EGFR). We can extract 
# the mapped symbols directly from our original SummarizedExperiment object.

up_symbols <- rowData(data_LUAD)[rownames(up_genes), "gene_name"]
down_symbols <- rowData(data_LUAD)[rownames(down_genes), "gene_name"]

# Remove any NA values and duplicates to ensure clean input
up_symbols <- na.omit(unique(up_symbols))
down_symbols <- na.omit(unique(down_symbols))

# 8. PERFORM ENRICHMENT ANALYSIS

# Run GO (Biological Process, Cellular Component, Molecular Function) 
# and KEGG pathway analysis for Upregulated genes
enrichment_UP <- TCGAanalyze_EAcomplete(
  TFname = "Upregulated_LUAD", 
  RegulonList = up_symbols
)

# Run the same analysis for Downregulated genes
enrichment_DOWN <- TCGAanalyze_EAcomplete(
  TFname = "Downregulated_LUAD", 
  RegulonList = down_symbols
)

# 9. VISUALIZE RESULTS: ENRICHMENT BARPLOTS

# Generate a multi-panel PDF showing the top 10 enriched pathways for UP genes
TCGAvisualize_EAbarplot(
  tf = rownames(enrichment_UP$ResBP),
  GOBPTab = enrichment_UP$ResBP,
  GOCCTab = enrichment_UP$ResCC,
  GOMFTab = enrichment_UP$ResMF,
  PathTab = enrichment_UP$ResPat,
  nRGTab = up_symbols,
  nBar = 10,
  filename = "TCGA_LUAD_Enrichment_UP.pdf"
)

# Generate the same plot for DOWN genes
TCGAvisualize_EAbarplot(
  tf = rownames(enrichment_DOWN$ResBP),
  GOBPTab = enrichment_DOWN$ResBP,
  GOCCTab = enrichment_DOWN$ResCC,
  GOMFTab = enrichment_DOWN$ResMF,
  PathTab = enrichment_DOWN$ResPat,
  nRGTab = down_symbols,
  nBar = 10,
  filename = "TCGA_LUAD_Enrichment_DOWN.pdf"
)