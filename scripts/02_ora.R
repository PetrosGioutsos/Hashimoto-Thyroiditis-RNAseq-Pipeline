
# --------------------------------------------------------
# Kindly note that my code is still under development
# --------------------------------------------------------

# Input from 00_setup.R
source("scripts/00_setup.R")



# 1. Load DESeq2 results
dds <- readRDS(file.path(paths$results, "dds.rds"))
sig_df <- readRDS(file.path(paths$results, "deseq2_results_sig.rds"))
message("Total DEGs loaded: ", nrow(sig_df))



# 2. Separate protein-coding from non-coding DEGs

# Go databases are limited for covering non-coding genes
# Therefore, ORA is applied exclusively for protein-coding genes.
sig_coding    <- sig_df[!is.na(sig_df$gene_type) & sig_df$gene_type == "protein_coding", ]

# Non-coding DEGs are saved separately for reference and potential future analysis.
sig_noncoding <- sig_df[ is.na(sig_df$gene_type) | sig_df$gene_type != "protein_coding", ]

message("Protein-coding DEGs: ", nrow(sig_coding))
message("Non-coding DEGs (saved separately): ", nrow(sig_noncoding))

saveRDS(sig_noncoding, file.path(paths$results, "deseq2_results_sig_noncoding.rds"))
if (nrow(sig_noncoding) > 0) {
  writexl::write_xlsx(sig_noncoding,
                      path = file.path(paths$results, "deseq2_results_sig_noncoding.xlsx"))
}
message("Non-coding DEGs saved to results/")



# 3. Creating gene lists for ORA

# Three separate analyses to preserve biological directionality:
# all_sig: all protein-coding DEGs
# up: up-regulated in disease
# down: down-regulated in disease

genes_all <- sig_coding$gene_id
genes_up <- sig_coding$gene_id[sig_coding$direction == "up_in_disease"]
genes_down <- sig_coding$gene_id[sig_coding$direction == "down_in_disease"]

message("Genes for ORA — All: ",  length(genes_all))
message("Up: ", length(genes_up))
message("Down: ", length(genes_down))



# 4. Background universe

# The universe must consist only of genes that were actually tested by DESeq2
# (i.e. genes that passed the rowSums >= 10 pre-filtering step), and consisting 
# of protein-coding genes only
coding_gene_ids <- unique(tx2gene$gene_id[tx2gene$gene_type == "protein_coding"])
universe <- intersect(rownames(dds), coding_gene_ids)

message("Universe: ", length(universe), " genes")



# 5. LFC vector for cnetplot coloring

# Named vector of log2FoldChanges for coloring gene nodes
# Named by gene symbol (gene_name), because readable = TRUE in enrichGO
# converts ENSEMBL IDs to gene symbols in the result — the names must match.
# Genes without a symbol fall back to their gene_id to avoid NA names.

lfc_vector <- sig_coding$log2FoldChange
names(lfc_vector) <- ifelse(
  is.na(sig_coding$gene_name) | sig_coding$gene_name == "",
  sig_coding$gene_id,
  sig_coding$gene_name
)



# 6. ORA helper function

# Run enrichGO for a given gene list and GO ontology
# Returns NULL with a warning if the gene list is too small or the call fails.
run_ora <- function(gene_ids, ontology, universe, label) {
  
  if (length(gene_ids) < 5) {
    warning("Skipping ORA [", label, " / ", ontology, "]: ",
            "fewer than 5 genes (n = ", length(gene_ids), ")")
    return(NULL)
  }
  
  result <- tryCatch(clusterProfiler::enrichGO(
    gene = gene_ids,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENSEMBL",   # ENSEMBL IDs as produced by DESeq2/tximport
      ont = ontology,    # "BP", "MF", or "CC"
      pAdjustMethod = "BH",
      pvalueCutoff = cfg$pvalue_cutoff,
      qvalueCutoff = cfg$qvalue_cutoff,
      readable = TRUE    # maps ENSEMBL IDs to gene symbols in the output
    ),
    error = function(e) {
      warning("enrichGO failed [", label, " / ", ontology, "]: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    message("No enriched terms found: [", label, " / ", ontology, "]")
    return(NULL)
  }
  
  # simplify() removes redundant GO terms that share large gene overlaps,
  result <- tryCatch(
    clusterProfiler::simplify(result, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) {
      warning("simplify() failed [", label, " / ", ontology, "]: ", conditionMessage(e),
              " — returning unsimplified result.")
      result
    }
  )
  
  message("Enriched terms [", label, " / ", ontology, "] after simplify(): ",
          nrow(as.data.frame(result)))
  return(result)
}



# 7. Run ORA — all combinations (3 gene lists x 3 ontologies = 9 runs)

ontologies <- c("BP", "MF", "CC")
gene_lists <- list(all_sig = genes_all, up = genes_up, down = genes_down)

# Results stored as a nested list: ora_results[[label]][[ontology]]
ora_results <- list()

for (label in names(gene_lists)) {
  ora_results[[label]] <- list()
  message("ORA: ", label, " (n = ", length(gene_lists[[label]]), ")")
  for (ont in ontologies) {
    ora_results[[label]][[ont]] <- run_ora(
      gene_ids = gene_lists[[label]],
      ontology = ont,
      universe = universe,
      label = label
    )
  }
}

saveRDS(ora_results, file.path(paths$results, "ora_results.rds"))
message("ORA results saved: results/ora_results.rds")



# 8. Export to Excel

# One sheet per non-null result, named as "<gene_list>_<ontology>" (e.g. "up_BP")
excel_sheets <- list()

for (label in names(ora_results)) {
  for (ont in ontologies) {
    res <- ora_results[[label]][[ont]]
    if (!is.null(res)) {
      excel_sheets[[paste0(label, "_", ont)]] <- as.data.frame(res)
    }
  }
}

if (length(excel_sheets) > 0) {
  writexl::write_xlsx(excel_sheets,
                      path = file.path(paths$results, "ora_results.xlsx"))
  message("ORA Excel saved: results/ora_results.xlsx")
} else {message("No enriched terms found across all comparisons. Excel not created.")}



# 9. Plots

save_plot <- function(filename, plot_expr, width = 8, height = 7) {
  filepath <- file.path(paths$plots, filename)
  png(filepath, width = width, height = height, units = "in", res = 300)
  tryCatch(force(plot_expr), finally = dev.off())
  message("Saved: ", filepath)
}

plot_labels <- list(
  all_sig = "All significant DEGs",
  up = "Up-regulated in disease",
  down = "Down-regulated in disease"
)

for (label in names(ora_results)) {
  for (ont in ontologies) {
    res <- ora_results[[label]][[ont]]
    if (is.null(res)) next
    n_terms <- nrow(as.data.frame(res))
    plot_title <- paste0("GO-", ont, " ORA \u2014 ", plot_labels[[label]])
    file_prefix <- paste0("02_ORA_", label, "_", ont)
    
    # --- Dotplot ---
    # Shows enrichment ratio, gene count, and adjusted p-value simultaneously.
    save_plot(
      paste0(file_prefix, "_dotplot.png"),
      print(
        clusterProfiler::dotplot(res,
                                 showCategory = min(20, n_terms),
                                 title = plot_title,
                                 font.size = 10
        ) +
          theme(plot.title = element_text(face = "bold", size = 11))
      ),
      width = 9, height = 8
    )
    
    # --- Barplot ---
    # Built manually with ggplot2 for compatibility across enrichplot versions.
    df_bar <- as.data.frame(res)
    df_bar <- head(df_bar[order(df_bar$p.adjust), ], min(20, n_terms))
    df_bar$Description <- factor(df_bar$Description, levels = rev(df_bar$Description))
    
    save_plot(paste0(file_prefix, "_barplot.png"),
      print(
        ggplot(df_bar, aes(x = Count, y = Description, fill = p.adjust)) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "#E74C3C", high = "#AEC6E8", name = "p.adjust") +
          labs(title = plot_title, x = "Gene Count", y = NULL) +
          theme_bw(base_size = 10) +
          theme(plot.title = element_text(face = "bold", size = 11))
      ),
      width = 9, height = 8
    )
    
    # --- cnetplot ---
    # Shows the gene-to-term network: which genes drive which enriched terms.
    # Limited to top 5 terms
    if (n_terms >= 3 && length(gene_lists[[label]]) >= 5) {
      save_plot(
        paste0(file_prefix, "_cnetplot.png"),
        print(
          clusterProfiler::cnetplot(res,
                                    foldChange = lfc_vector,
                                    showCategory = min(5, n_terms),
                                    node_label = "all"
          ) +
            labs(title = paste0("cnetplot \u2014 ", plot_title)) +
            theme(plot.title = element_text(face = "bold", size = 11))
        ),
        width = 12, height = 10
      )
    }
    
  }
}



# 10. Summary table

summary_rows <- lapply(names(ora_results), function(label) {
  lapply(ontologies, function(ont) {
    res <- ora_results[[label]][[ont]]
    data.frame(
      gene_list = label,
      ontology = ont,
      n_input = length(gene_lists[[label]]),
      n_terms = if (is.null(res)) 0L else nrow(as.data.frame(res)),
      stringsAsFactors = FALSE
    )
  })
})

summary_df <- do.call(rbind, do.call(c, summary_rows))
rownames(summary_df) <- NULL

saveRDS(summary_df, file.path(paths$results, "ora_summary.rds"))
writexl::write_xlsx(summary_df, path = file.path(paths$results, "ora_summary.xlsx"))

print(summary_df)
message("ORA analysis complete.")
