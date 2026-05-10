
# This script must be run at the beginning of every downstream script

# --------------------------------------------------------
# Kindly note that my code is still under development
# --------------------------------------------------------



# 1. Package Installation

# Bioconductor packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "tximport",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "pathview",
  "fgsea",
  "EnhancedVolcano",
  "ComplexHeatmap",
  "WGCNA",
  "STRINGdb",
  "apeglm",
  "rtracklayer"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# CRAN packages installation
cran_pkgs <- c(
  "tidyverse",
  "ggrepel",
  "RColorBrewer",
  "scales",
  "pheatmap",
  "patchwork",
  "ggplotify",
  "circlize",
  "writexl",
  "knitr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg)
  }
}

# Importing all packages
all_pkgs <- c(bioc_pkgs, cran_pkgs)
invisible(lapply(all_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

message("All libraries have been imported successfully")



# 2. Setting input and output using relative paths.

# Note: Relative paths require the .Proj file to be located at the project root 
# Folder structure
#     data/quants/ -> one directory per sample (eg SRR29759248_quant), as it was created by Salmon
#     data/reference/ -> GENCODE GTF annotation file

paths <- list(
  # Input data
  quants = file.path("data", "quants"),
  reference = file.path("data", "reference"),
  
  # Output directories
  results = "results",
  plots = "plots",
  scripts = "scripts"
)

# Output directories are automatically created by running the following
for (d in c(paths$results, paths$plots)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message(paste("Created directory:", d))
  }
}

# The script will warn you if the GTF file is not located in the path 
# that is should be 
paths$gtf <- file.path(paths$reference, "gencode.v49.annotation.gtf.gz")
if (!file.exists(paths$gtf))
  warning("The GTF file is not in the following path: ", paths$gtf)



# 3. Sample Metadata

# The sample size consists of 4 patients (disease, Hashimoto's Thyroiditis) 
# and 5 healthy individuals

sample_ids <- c(
  # Disease — Hashimoto's Thyroiditis
  "SRR29759248", "SRR29759249", "SRR29759250", "SRR29759251",
  # Healthy controls
  "SRR29759252", "SRR29759253", "SRR29759254", "SRR29759262", "SRR29759263"
)

condition <- factor(
  c(rep("disease", 4), rep("healthy", 5)),
  levels = c("healthy", "disease")   # Healthy individuals are set as reference
)

# Creating a metadata table that will be used by the following scripts. 
sample_info <- data.frame(
  sample = sample_ids,
  condition = condition,
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

# Paths to the quant.sf files
quant_files <- file.path(paths$quants,
                         paste0(sample_ids, "_quant"),
                         "quant.sf")
names(quant_files) <- sample_ids

# The script will notify you for the quant.sf files
missing_files <- quant_files[!file.exists(quant_files)]
if (length(missing_files) > 0) {
  warning("The following quant.sf files could not be located:\n",
          paste(missing_files, collapse = "\n"))
} else {
  message("All quant.sf files were located (", length(quant_files), " samples).")
}



# 4. Parameters (thresholds for the next scripts)

cfg <- list(
  # DESeq2
  padj_cutoff = 0.05,
  lfc_cutoff = 1.0
  )



# 5. Color Palette (for the plots, to be decided)
palette <- list(
  condition = c(
    healthy = "#2196F3",
    disease = "#F44336"
  ),
  heatmap = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  volcano_up = "#E74C3C",   # upregulated genes
  volcano_down = "#3498DB",   # downregulated genes
  volcano_ns = "grey70"     # non-significant genes
)



# 6. Transcript to Gene mapping

# Neccessary for tximport 

# Converting Transcript IDs (as created by Salmon) to Gene IDs
# It is only executed once. It is stored for future use (optimizing RAM usage)
tx2gene_path <- file.path(paths$results, "tx2gene.rds")

if (!file.exists(tx2gene_path)) {
  message("Creating tx2gene mapping from GTF")
  
  gtf_gr <- rtracklayer::import(paths$gtf, feature.type = "transcript") # Creating an object (GRanges) representing gene loci
  gtf_df <- as.data.frame(gtf_gr)[, c("transcript_id", "gene_id", "gene_name", "gene_type")] # Converting GRanges to a data frame
  
  # Removing the suffix from Transcript IDs, ensuring that the IDs are
  # written without their version
  # Note: the --gencode flag was used with Salmon. Thus, only base IDs have been kept
  gtf_df$transcript_id_base <- sub("\\..*", "", gtf_df$transcript_id)
  gtf_df$gene_id_base <- sub("\\..*", "", gtf_df$gene_id)
  
  tx2gene <- dplyr::distinct(gtf_df[, c("transcript_id_base", "gene_id_base", "gene_name", "gene_type")])
  colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name", "gene_type") # The first two columns have to be the transcript IDs and gene IDs for the tximport library that will follow
  
  # Saving the object (GRanges) for later use 
  saveRDS(tx2gene, tx2gene_path)
  message("tx2gene has been saved: ", tx2gene_path)
} else {
  tx2gene <- readRDS(tx2gene_path)
  message("tx2gene loaded from cache: ", tx2gene_path)
}



# 7. Session info 
session_info_path <- file.path(paths$results, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_path)
message("Session info has been saved: ", session_info_path)

