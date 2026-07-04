# pulls in everything set up in 00_setup.R
source("scripts/00_setup.R")



# 1. Load DESeq2 results

dds <- readRDS(file.path(paths$results, "dds.rds"))
sig_df <- readRDS(file.path(paths$results, "deseq2_results_sig.rds"))
message("Total DEGs loaded: ", nrow(sig_df))



# 2. Protein-coding DEGs — up and down lists

# KEGG coverage for non-coding genes is limited/inconsistent, so this only
# runs on protein-coding DEGs
sig_coding <- sig_df[!is.na(sig_df$gene_type) & sig_df$gene_type == "protein_coding", ]
message("Protein-coding DEGs: ", nrow(sig_coding))

genes_up <- sig_coding$gene_id[sig_coding$direction == "up_in_disease"]
genes_down <- sig_coding$gene_id[sig_coding$direction == "down_in_disease"]

message("Up-regulated: ", length(genes_up))
message("Down-regulated: ", length(genes_down))



# 3. ENSEMBL to Entrez ID conversion

# enrichKEGG() only takes Entrez Gene IDs, so bitr() does the mapping.
# anything without a valid Entrez entry just gets dropped

convert_to_entrez <- function(ensembl_ids, label) {
  mapped <- clusterProfiler::bitr(
    ensembl_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  n_lost <- length(ensembl_ids) - nrow(mapped)
  if (n_lost > 0)
    message("  [", label, "] ", n_lost,
            " IDs could not be mapped to Entrez and were dropped.")
  message("  [", label, "] ", nrow(mapped), " IDs successfully mapped.")
  mapped
}

map_up <- convert_to_entrez(genes_up,   "up")
map_down <- convert_to_entrez(genes_down, "down")

entrez_up <- map_up$ENTREZID
entrez_down <- map_down$ENTREZID



# 4. Background universe

# same idea as in the ORA script: protein-coding genes that passed DESeq2's
# rowSums >= 10 filter, converted to Entrez IDs. keeping the universe to
# genes that were actually tested avoids inflating the enrichment stats
# with genes that were never measured in this experiment

coding_gene_ids <- unique(tx2gene$gene_id[tx2gene$gene_type == "protein_coding"])
universe_ensembl <- intersect(rownames(dds), coding_gene_ids)

universe_map <- clusterProfiler::bitr(
  universe_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
universe_entrez <- universe_map$ENTREZID
message("Universe (Entrez): ", length(universe_entrez), " genes")



# 5. LFC vectors (Entrez-named) — for pathview and cnetplot

# both pathview and cnetplot want a named numeric vector:
#   names  = Entrez Gene IDs
#   values = log2FoldChange

make_lfc_entrez <- function(mapping_df, sig_coding_df) {
  merged <- merge(
    mapping_df,
    sig_coding_df[, c("gene_id", "log2FoldChange")],
    by.x = "ENSEMBL", by.y = "gene_id"
  )
  lfc <- merged$log2FoldChange
  names(lfc) <- merged$ENTREZID
  lfc
}

lfc_up_entrez <- make_lfc_entrez(map_up, sig_coding)
lfc_down_entrez <- make_lfc_entrez(map_down, sig_coding)

# combined vector (all sig DEGs) for pathview coloring — this way the whole
# expression context shows up on the pathway diagram, not just the subset
# of genes that drove the enrichment
lfc_all_entrez <- c(lfc_up_entrez, lfc_down_entrez)



# 6. Gene-symbol-named LFC vectors — for cnetplot

# after setReadable(), the enrichment tables show gene symbols instead of
# Entrez IDs, so cnetplot needs symbol-named LFCs to match

make_lfc_symbol <- function(mapping_df, sig_coding_df) {
  merged <- merge(
    mapping_df,
    sig_coding_df[, c("gene_id", "log2FoldChange", "gene_name")],
    by.x = "ENSEMBL", by.y = "gene_id"
  )
  lfc <- merged$log2FoldChange
  # fall back to Entrez ID if there's no gene symbol
  names(lfc) <- ifelse(
    is.na(merged$gene_name) | merged$gene_name == "",
    merged$ENTREZID,
    merged$gene_name
  )
  lfc
}

lfc_up_symbol <- make_lfc_symbol(map_up, sig_coding)
lfc_down_symbol <- make_lfc_symbol(map_down, sig_coding)



# 7. Save plot helper (same as in 01 and 02)

save_plot <- function(filename, plot_expr, width = 10, height = 8) {
  filepath <- file.path(paths$plots, filename)
  png(filepath, width = width, height = height, units = "in", res = 300)
  tryCatch(force(plot_expr), finally = dev.off())
  message("  Saved: ", filepath)
}



# 8. KEGG enrichment helper

run_kegg <- function(entrez_ids, universe_entrez, label) {
  
  if (length(entrez_ids) < 5) {
    warning("Skipping KEGG [", label, "]: fewer than 5 genes (n = ",
            length(entrez_ids), ")")
    return(NULL)
  }
  
  result <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = entrez_ids,
      universe = universe_entrez,
      organism = "hsa",           # human
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = cfg$pvalue_cutoff,
      qvalueCutoff = cfg$qvalue_cutoff
    ),
    error = function(e) {
      warning("enrichKEGG failed [", label, "]: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    message("No enriched KEGG pathways found: [", label, "]")
    return(NULL)
  }
  
  # setReadable() swaps in gene symbols for the Entrez IDs in the result table
  result <- tryCatch(
    clusterProfiler::setReadable(result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) {
      warning("setReadable() failed [", label, "]: ", conditionMessage(e),
              " — keeping Entrez IDs.")
      result
    }
  )
  
  message("Enriched KEGG pathways [", label, "]: ", nrow(as.data.frame(result)))
  result
}



# 9. Run KEGG enrichment

kegg_up <- run_kegg(entrez_up,   universe_entrez, "up")
kegg_down <- run_kegg(entrez_down, universe_entrez, "down")

# named list, so the loops below can plot/export both at once
pathway_results <- list(
  up_KEGG = kegg_up,
  down_KEGG = kegg_down
)

saveRDS(pathway_results, file.path(paths$results, "pathway_results.rds"))
message("Pathway results saved: results/pathway_results.rds")



# 10. Plots: dotplot, barplot, cnetplot

plot_labels <- list(
  up_KEGG = "Up-regulated in disease \u2014 KEGG",
  down_KEGG = "Down-regulated in disease \u2014 KEGG"
)

lfc_symbol_list <- list(
  up_KEGG = lfc_up_symbol,
  down_KEGG = lfc_down_symbol
)

for (result_name in names(pathway_results)) {
  res <- pathway_results[[result_name]]
  if (is.null(res)) next
  
  n_terms <- nrow(as.data.frame(res))
  plot_title <- plot_labels[[result_name]]
  prefix <- paste0("03_Pathway_", result_name)
  lfc_vec <- lfc_symbol_list[[result_name]]
  
  # Dotplot
  # enrichment ratio, gene count and adjusted p-value in one view
  save_plot(
    paste0(prefix, "_dotplot.png"),
    print(clusterProfiler::dotplot(res,
                                   showCategory = min(20, n_terms),
                                   title = plot_title,
                                   font.size = 10
    ) +
      theme(plot.title = element_text(face = "bold", size = 11))
    )
  )
  
  # Barplot
  df_bar <- as.data.frame(res)
  df_bar <- head(df_bar[order(df_bar$p.adjust), ], min(20, n_terms))
  df_bar$Description <- factor(df_bar$Description, levels = rev(df_bar$Description))
  
  save_plot(
    paste0(prefix, "_barplot.png"),
    print(ggplot(df_bar, aes(x = Count, y = Description, fill = p.adjust)) +
            geom_bar(stat = "identity") +
            scale_fill_gradient(low = "#E74C3C",
                                high = "#AEC6E8",
                                name = "p.adjust") +
            labs(title = plot_title, x = "Gene Count", y = NULL) +
            theme_bw(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 11))
    )
  )
  
  # cnetplot
  # which genes drive which pathways (top 5). node color = log2FC,
  # node size = gene count
  if (n_terms >= 3 && length(lfc_vec) >= 5) {
    save_plot(
      paste0(prefix, "_cnetplot.png"),
      print(clusterProfiler::cnetplot(res,
                                      foldChange = lfc_vec,
                                      showCategory = n_terms,
                                      node_label = "all"
      ) +
        labs(title = paste0("cnetplot \u2014 ", plot_title)) +
        theme(plot.title = element_text(face = "bold", size = 11))
      ),
      width = 13, height = 11
    )
  }
}



# 11. Pathview diagrams (top 5 KEGG pathways each, up and down)

# pathview() pulls the KEGG pathway image straight from the internet,
# overlays the LFC values as node colors, and writes out a PNG.
# NOTE: pathview always writes into the current working directory,
# not paths$plots — this function just moves the files afterwards

run_pathview <- function(kegg_result, lfc_entrez_vec, label, n_top = 5) {
  
  if (is.null(kegg_result)) {
    message("Skipping pathview [", label, "]: no KEGG results available.")
    return(invisible(NULL))
  }
  
  df_kegg <- as.data.frame(kegg_result)
  top_pathways <- head(df_kegg[order(df_kegg$p.adjust), "ID"], n_top)
  
  message("Running pathview [", label, "] for ", length(top_pathways), " pathways...")
  
  for (pid in top_pathways) {
    tryCatch({
      
      pathview::pathview(
        gene.data = lfc_entrez_vec,  # named vector: Entrez ID -> LFC
        pathway.id = pid,
        species = "hsa",
        out.suffix = label,           # gets appended to the output filename
        kegg.native = TRUE,            # keeps the original KEGG pathway layout
        
        # blue (down) -> white (no change) -> red (up)
        low = list(gene = palette$volcano_down),
        mid = list(gene = "white"),
        high = list(gene = palette$volcano_up),
        key.pos = "topright"
      )
      
      # pathview names its output "<pid>.<label>.png" (sometimes a PDF too)
      # and drops it in the working directory — move it into plots/ instead
      generated_files <- list.files(
        pattern = paste0("^", pid, "\\.", label),
        full.names = TRUE
      )
      
      if (length(generated_files) == 0) {
        warning("No output files found for pathway: ", pid)
      }
      
      for (f in generated_files) {
        dest <- file.path(paths$plots, basename(f))
        file.rename(f, dest)
        message("  Moved: ", basename(f), " \u2192 plots/")
      }
      
    }, error = function(e) {
      warning("pathview failed [", pid, " / ", label, "]: ", conditionMessage(e))
    })
  }
}

run_pathview(kegg_up, lfc_all_entrez, label = "up", n_top = 5)
run_pathview(kegg_down, lfc_all_entrez, label = "down", n_top = 5)



# 12. Export to Excel (one sheet per non-null result)

excel_sheets <- list()

for (result_name in names(pathway_results)) {
  res <- pathway_results[[result_name]]
  if (!is.null(res))
    excel_sheets[[result_name]] <- as.data.frame(res)
}

if (length(excel_sheets) > 0) {
  writexl::write_xlsx(
    excel_sheets,
    path = file.path(paths$results, "pathway_results.xlsx")
  )
  message("Pathway Excel saved: results/pathway_results.xlsx")
} else {
  message("No enriched pathways found across all comparisons. Excel not created.")
}



# 13. Summary table
summary_df <- do.call(rbind, lapply(names(pathway_results), function(nm) {
  res <- pathway_results[[nm]]
  data.frame(
    analysis = nm,
    n_terms = if (is.null(res)) 0L else nrow(as.data.frame(res)),
    stringsAsFactors = FALSE
  )
}))
rownames(summary_df) <- NULL

saveRDS(summary_df, file.path(paths$results, "pathway_summary.rds"))
writexl::write_xlsx(summary_df,
                    path = file.path(paths$results, "pathway_summary.xlsx"))

print(summary_df)
message("Pathway analysis complete")