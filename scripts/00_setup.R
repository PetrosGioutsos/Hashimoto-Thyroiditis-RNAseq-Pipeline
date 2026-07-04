# runs first — everything else in the pipeline sources this


# 1. Packages

# Bioconductor packages
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

# CRAN packages
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
  "knitr",
  "msigdbr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg)
  }
}

# load everything
all_pkgs <- c(bioc_pkgs, cran_pkgs)
invisible(lapply(all_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

message("All libraries have been imported successfully")



# 2. Paths

# relative paths, so this only works if the .Rproj is at the project root
# folder structure:
#     data/quants/ -> one folder per sample (e.g. SRR29759248_quant), as created by Salmon
#     data/reference/ -> GENCODE GTF file

paths <- list(
  # input data
  quants = file.path("data", "quants"),
  reference = file.path("data", "reference"),
  
  # output
  results = "results",
  plots = "plots",
  scripts = "scripts"
)

# create the output folders if they don't already exist
for (d in c(paths$results, paths$plots)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message(paste("Created directory:", d))
  }
}

# heads-up if the GTF isn't where it should be
paths$gtf <- file.path(paths$reference, "gencode.v49.annotation.gtf.gz")
if (!file.exists(paths$gtf))
  warning("The GTF file is not in the following path: ", paths$gtf)



# 3. Sample metadata

# 4 HT patients + 4 healthy controls (originally 5 healthy samples were
# sequenced, but SRR29759263 turned out to be an outlier in the clustering
# and was dropped, so it's not in the list below)

sample_ids <- c(
  # disease — Hashimoto's Thyroiditis
  "SRR29759248", "SRR29759249", "SRR29759250", "SRR29759251",
  # healthy controls (SRR29759263 excluded — outlier in hierarchical clustering)
  "SRR29759252", "SRR29759253", "SRR29759254", "SRR29759262"
)

condition <- factor(
  c(rep("disease", 4), rep("healthy", 4)),
  levels = c("healthy", "disease")   # healthy = reference level
)

# metadata table used by all the scripts downstream
sample_info <- data.frame(
  sample = sample_ids,
  condition = condition,
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

# paths to the quant.sf files
quant_files <- file.path(paths$quants,
                         paste0(sample_ids, "_quant"),
                         "quant.sf")
names(quant_files) <- sample_ids

# check they're actually all there
missing_files <- quant_files[!file.exists(quant_files)]
if (length(missing_files) > 0) {
  warning("The following quant.sf files could not be located:\n",
          paste(missing_files, collapse = "\n"))
} else {
  message("All quant.sf files were located (", length(quant_files), " samples).")
}



# 4. Parameters (used as thresholds in the other scripts)

cfg <- list(
  # DESeq2
  padj_cutoff = 0.05,
  lfc_cutoff = 1.0,
  
  # ORA / Pathway / GSEA
  pvalue_cutoff  = 0.05,   # p-value
  qvalue_cutoff  = 0.20   # q-value (FDR)
)


# 5. Color palette (for the plots)
palette <- list(
  condition = c(
    healthy = "#2196F3",
    disease = "#F44336"
  ),
  heatmap = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  volcano_up = "#E74C3C",   # up
  volcano_down = "#3498DB",   # down
  volcano_ns = "grey70"     # not significant
)



# 6. Transcript-to-gene mapping

# tximport needs this — maps transcript IDs (from Salmon) to gene IDs.
# Only built once, then cached (saves time and RAM on reruns)
tx2gene_path <- file.path(paths$results, "tx2gene.rds")

if (!file.exists(tx2gene_path)) {
  message("Creating tx2gene mapping from GTF")
  
  gtf_gr <- rtracklayer::import(paths$gtf, feature.type = "transcript") # GRanges object, transcript-level
  gtf_df <- as.data.frame(gtf_gr)[, c("transcript_id", "gene_id", "gene_name", "gene_type")] # easier to work with as a df
  
  # strip the version suffix from the IDs
  # (Salmon was run with --gencode, so this mostly just tidies things up)
  gtf_df$transcript_id_base <- sub("\\..*", "", gtf_df$transcript_id)
  gtf_df$gene_id_base <- sub("\\..*", "", gtf_df$gene_id)
  
  tx2gene <- dplyr::distinct(gtf_df[, c("transcript_id_base", "gene_id_base", "gene_name", "gene_type")])
  colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name", "gene_type") # tximport expects these two columns first
  
  # cache it
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