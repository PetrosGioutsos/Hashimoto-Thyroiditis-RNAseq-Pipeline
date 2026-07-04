# pulls in everything set up in 00_setup.R
source("scripts/00_setup.R")



# 1. Load DESeq2 results

# shrunk results table, all tested genes (not just the significant ones) —
# GSEA needs the full ranked list, no cutoff applied here.
# Subramanian et al. (2005, PNAS) are explicit about this: GSEA
# "considers all of the genes in an experiment, not only those above an
# arbitrary cutoff in terms of fold-change or significance."
# LFC and lfcSE both come from lfcShrink(dds, type = "apeglm"), consistent
# with the Wald test formulation in Love et al. (2014, Genome Biology).
res_df <- readRDS(file.path(paths$results, "deseq2_results_full.rds"))
message("Genes loaded for ranking: ", nrow(res_df))



# 2. Build the ranked list

# genes are ranked by the Wald statistic: shrunk_stat = LFC_shrunk / SE_shrunk.
# both come from lfcShrink(dds, type = "apeglm"), following Love et al.
# (2014, Genome Biology) — this factors in estimation uncertainty while
# still being grounded in a peer-reviewed shrinkage framework.
# genes with no valid shrunk_stat or gene symbol get dropped.
# duplicate gene symbols are resolved by keeping whichever entry has the
# larger absolute shrunk_stat (standard practice)

res_df$shrunk_stat <- res_df$log2FoldChange / res_df$lfcSE

rank_df <- res_df[!is.na(res_df$shrunk_stat) & !is.na(res_df$gene_name), ]

# resolve duplicate gene symbols
rank_df <- rank_df[order(-abs(rank_df$shrunk_stat)), ]
rank_df <- rank_df[!duplicated(rank_df$gene_name), ]

ranked_list <- rank_df$shrunk_stat
names(ranked_list) <- rank_df$gene_name
ranked_list <- sort(ranked_list, decreasing = TRUE)

message("Genes in ranked list: ", length(ranked_list))



# 3. Gene set collections via msigdbr

# using two collections:
#   H  — Hallmark: 50 well-defined, low-noise biological states
#   C2 — KEGG: lets us compare directly against the ORA-based KEGG in 03_pathway.R

get_gene_sets <- function(collection, subcollection = NULL) {
  df <- msigdbr(species = "Homo sapiens",
                collection = collection,
                subcollection = subcollection)
  split(df$gene_symbol, df$gs_name)
}

gene_sets <- list(
  Hallmark = get_gene_sets("H"),
  KEGG = get_gene_sets("C2", "CP:KEGG_LEGACY")
)

message("Gene set collection sizes — Hallmark: ", length(gene_sets$Hallmark),
        " | KEGG: ", length(gene_sets$KEGG))



# 4. Save plot helper (same as the earlier scripts)

save_plot <- function(filename, plot_expr, width = 10, height = 8) {
  filepath <- file.path(paths$plots, filename)
  png(filepath, width = width, height = height, units = "in", res = 300)
  tryCatch(force(plot_expr), finally = dev.off())
  message("  Saved: ", filepath)
}



# 5. GSEA helper

# runs fgsea::fgseaMultilevel on one gene set collection.
# minSize = 15 and maxSize = 500 follow Subramanian et al. (2005), who
# recommend focusing on robust signals rather than tiny/huge gene sets.
# returns NULL with a warning if the call fails

run_gsea <- function(ranked_list, pathways, label) {
  
  result <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = ranked_list,
      minSize = 15,
      maxSize = 500,
      eps = 0          # lets very low p-values be estimated properly
    ),
    error = function(e) {
      warning("fgseaMultilevel failed [", label, "]: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(result) || nrow(result) == 0) {
    message("No results returned by fgsea [", label, "]")
    return(NULL)
  }
  
  # sort by NES, add direction label
  result <- result[order(result$padj, -abs(result$NES)), ]
  result$direction <- ifelse(result$NES > 0, "up_in_disease", "down_in_disease")
  
  n_sig <- sum(result$padj < cfg$pvalue_cutoff, na.rm = TRUE)
  message("Significant gene sets [", label, "] (padj < ", cfg$pvalue_cutoff, "): ", n_sig)
  result
}



# 6. Run GSEA for both collections

set.seed(42)

gsea_results <- list(
  Hallmark = run_gsea(ranked_list, gene_sets$Hallmark, "Hallmark"),
  KEGG = run_gsea(ranked_list, gene_sets$KEGG, "KEGG")
)

saveRDS(gsea_results, file.path(paths$results, "gsea_results.rds"))
message("GSEA results saved: results/gsea_results.rds")



# 7. Plots

collection_labels <- list(
  Hallmark = "Hallmark",
  KEGG = "KEGG"
)

for (collection in names(gsea_results)) {
  
  res <- gsea_results[[collection]]
  if (is.null(res)) next
  
  gs <- gene_sets[[collection]]
  label <- collection_labels[[collection]]
  prefix <- paste0("04_GSEA_", collection)
  
  # significant gene sets only, for most of the plots below
  res_sig <- res[!is.na(res$padj) & res$padj < cfg$pvalue_cutoff, ]
  n_sig <- nrow(res_sig)
  
  if (n_sig == 0) {
    message("No significant gene sets to plot [", collection, "]")
    next
  }
  
  
  # Dotplot (NES vs pathway, coloured by padj)
  # direction (NES), significance (padj) and gene set size all in one view.
  # top 20 significant gene sets
  top_sets <- head(res_sig[order(res_sig$padj), ], 20)
  
  save_plot(
    paste0(prefix, "_dotplot.png"),
    {
      df_dot <- as.data.frame(top_sets)
      df_dot$pathway <- factor(df_dot$pathway, levels = rev(df_dot$pathway))
      print(
        ggplot(df_dot, aes(x = NES, y = pathway, colour = padj, size = size)) +
          geom_point() +
          scale_colour_gradient(low = "#E74C3C", high = "#AEC6E8", name = "padj") +
          scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
          geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
          labs(
            title = paste0("GSEA \u2014 ", label),
            subtitle = paste0("Top ", nrow(df_dot), " significant gene sets (padj < ",
                              cfg$pvalue_cutoff, ")"),
            x = "Normalized Enrichment Score (NES)",
            y = NULL
          ) +
          theme_bw(base_size = 10) +
          theme(plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 9, colour = "grey40"))
      )
    },
    width = 11, height = max(5, 0.35 * nrow(top_sets) + 3)
  )
  
  
  # Barplot (NES, coloured by direction)
  # signed bars: up-regulated sets to the right, down-regulated to the left
  save_plot(
    paste0(prefix, "_barplot.png"),
    {
      df_bar <- as.data.frame(top_sets)
      df_bar$pathway <- factor(df_bar$pathway, levels = df_bar$pathway[order(df_bar$NES)])
      print(
        ggplot(df_bar, aes(x = NES, y = pathway,
                           fill = ifelse(NES > 0, "up_in_disease", "down_in_disease"))) +
          geom_bar(stat = "identity") +
          scale_fill_manual(
            values = c(up_in_disease = palette$volcano_up,
                       down_in_disease = palette$volcano_down),
            name = "Direction"
          ) +
          geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
          labs(
            title = paste0("GSEA \u2014 ", label),
            subtitle = paste0("Top ", nrow(df_bar), " significant gene sets"),
            x = "Normalized Enrichment Score (NES)",
            y = NULL
          ) +
          theme_bw(base_size = 10) +
          theme(plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 9, colour = "grey40"),
                legend.position = "bottom")
      )
    },
    width = 11, height = max(5, 0.35 * nrow(top_sets) + 3)
  )
  
  
  # Enrichment plots (top 3 up, top 3 down)
  # the running-sum enrichment plot is the classic GSEA visual
  # (Subramanian et al. 2005, Fig. 1B) — shows where the gene set members
  # fall across the ranked list
  top_up <- head(res_sig$pathway[res_sig$direction == "up_in_disease"],   3)
  top_down <- head(res_sig$pathway[res_sig$direction == "down_in_disease"], 3)
  top_for_nes <- c(top_up, top_down)
  
  for (pw in top_for_nes) {
    direction_tag <- ifelse(pw %in% top_up, "up", "down")
    safe_name <- gsub("[^A-Za-z0-9_]", "_", pw)
    fname <- paste0(prefix, "_enrichmentplot_", direction_tag, "_", safe_name, ".png")
    
    save_plot(
      fname,
      {
        p <- fgsea::plotEnrichment(gs[[pw]], ranked_list) +
          labs(
            title = pw,
            subtitle = paste0(label, " \u2014 ",
                              ifelse(direction_tag == "up",
                                     "Up-regulated in disease",
                                     "Down-regulated in disease")),
            x = "Gene rank",
            y = "Enrichment score"
          ) +
          theme_bw(base_size = 11) +
          theme(plot.title = element_text(face = "bold", size = 10),
                plot.subtitle = element_text(size = 9, colour = "grey40"))
        print(p)
      },
      width = 9, height = 5
    )
  }
}


# 8. Export to Excel (one sheet per collection)

excel_sheets <- list()

for (collection in names(gsea_results)) {
  res <- gsea_results[[collection]]
  if (is.null(res)) next
  
  # leadingEdge is a list column — flatten it to a semicolon-separated
  # string so Excel can actually display it
  df_out <- as.data.frame(res)
  df_out$leadingEdge <- vapply(
    df_out$leadingEdge,
    function(x) paste(x, collapse = "; "),
    character(1)
  )
  excel_sheets[[collection]] <- df_out
}

if (length(excel_sheets) > 0) {
  writexl::write_xlsx(excel_sheets,
                      path = file.path(paths$results, "gsea_results.xlsx"))
  message("GSEA Excel saved: results/gsea_results.xlsx")
} else {
  message("No GSEA results to export.")
}



# 9. Summary table

summary_df <- do.call(rbind, lapply(names(gsea_results), function(nm) {
  res <- gsea_results[[nm]]
  if (is.null(res)) {
    data.frame(collection = nm, n_tested = 0L,
               n_sig = 0L, n_up = 0L, n_down = 0L,
               stringsAsFactors = FALSE)
  } else {
    sig <- res[!is.na(res$padj) & res$padj < cfg$pvalue_cutoff, ]
    data.frame(
      collection = nm,
      n_tested = nrow(res),
      n_sig = nrow(sig),
      n_up = sum(sig$direction == "up_in_disease"),
      n_down = sum(sig$direction == "down_in_disease"),
      stringsAsFactors = FALSE
    )
  }
}))
rownames(summary_df) <- NULL

saveRDS(summary_df, file.path(paths$results, "gsea_summary.rds"))
writexl::write_xlsx(summary_df,
                    path = file.path(paths$results, "gsea_summary.xlsx"))

print(summary_df)
message("GSEA analysis complete.")