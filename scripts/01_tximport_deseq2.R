
# --------------------------------------------------------
# Kindly note that my code is still under development
# --------------------------------------------------------

# Input from 00_setup.R
source("scripts/00_setup.R")



# 1. tximport 

# tximport summarizes transcript-level TPMs and estimated counts to gene level
# quant.sf formats were created by Salmon. 
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene[, c("tx_id", "gene_id")],
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"
)
message("Genes quantified: ", nrow(txi$counts))
message("Samples: ", ncol(txi$counts))



# 2. DESeqDataSet construction

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_info,
  design = ~ condition   # healthy = reference
)
message("Genes before filtering: ", nrow(dds))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

message("Genes after filtering: ", nrow(dds), "  (Removed: ", sum(!keep), ")")



# 3. DESeq2 analysis

set.seed(42)
dds <- DESeq(dds)



# 4. VST transformation for visualization (PCA, heatmap)

vst <- vst(dds, blind = FALSE)  # uses dispersion estimates from the fitted model



# 5. LFC shrinkage 

# apeglm shrinkage reduces noise on low-count genes

res_names <- resultsNames(dds)
message("Available coefficients: ", paste(res_names, collapse = ", "))

coef_name <- "condition_disease_vs_healthy"

if (!coef_name %in% res_names) {
  stop("Coefficient '", coef_name, "' not found in resultsNames(dds). ",
       "Check that the reference level in sample_info is set to 'healthy'.")
}

res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm"
)

# un-shrunk results 
res_raw <- results(dds,
  contrast = c("condition", "disease", "healthy"),
  alpha = cfg$padj_cutoff
)

message("Genes with padj < ", cfg$padj_cutoff, " : ",
        sum(res_shrunk$padj < cfg$padj_cutoff, na.rm = TRUE))
message("... whith an |LFC| >= ", cfg$lfc_cutoff, " : ",
        sum(res_shrunk$padj < cfg$padj_cutoff &
              abs(res_shrunk$log2FoldChange) >= cfg$lfc_cutoff, na.rm = TRUE))



# 6. Results annotation 

gene_lookup <- dplyr::distinct(tx2gene[, c("gene_id", "gene_name", "gene_type")])

# converting to data frame and annotation
res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
res_df <- dplyr::left_join(res_df, gene_lookup, by = "gene_id")

# selecting only columns that exist
desired_order <- c("gene_id", "gene_name", "gene_type", "baseMean", "log2FoldChange", 
                   "lfcSE", "svalue", "pvalue", "padj")
res_df <- res_df[, desired_order[desired_order %in% colnames(res_df)]]

# sort by padj, then by |LFC|
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]
rownames(res_df) <- NULL

# Significant DEGs
sig_df <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < cfg$padj_cutoff &
    abs(res_df$log2FoldChange) >= cfg$lfc_cutoff, ]

sig_df$direction <- ifelse(sig_df$log2FoldChange > 0, "up_in_disease", "down_in_disease")

message("Total significant DEGs: ", nrow(sig_df))
message("Up-regulated in disease: ", sum(sig_df$direction == "up_in_disease"))
message("Down-regulated in disease: ", sum(sig_df$direction == "down_in_disease"))

# save RDS
saveRDS(dds, file.path(paths$results, "dds.rds"))
saveRDS(vst, file.path(paths$results, "vst.rds"))
saveRDS(res_df, file.path(paths$results, "deseq2_results_full.rds"))
saveRDS(sig_df, file.path(paths$results, "deseq2_results_sig.rds"))

message("RDS files saved to results/")

# save excel
writexl::write_xlsx(
  list(All_genes = res_df, Significant_DEGs = sig_df),
  path = file.path(paths$results, "deseq2_results.xlsx")
)

# 7. Plots

# save plot as PDF
save_plot <- function(filename, plot_expr, width = 7, height = 6) {   
  filepath <- file.path(paths$plots, filename)
  pdf(filepath, width = width, height = height)
  tryCatch(force(plot_expr), finally = dev.off())
  message("  Saved: ", filepath)
}

# PCA plot
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition, label  = name)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.2, show.legend = FALSE, max.overlaps = 20) +
  scale_colour_manual(values = palette$condition) +
  labs(
    title = "PCA — VST-normalised counts",
    subtitle = "Top 500 most variable genes",
    x = paste0("PC1  (", pct_var[1], "% variance)"),
    y = paste0("PC2  (", pct_var[2], "% variance)"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

save_plot("01_PCA_plot.pdf", print(pca_plot), width = 7, height = 6)

# MA plot on shrunk LFCs
save_plot("01_MA_plot.pdf",
  {
    DESeq2::plotMA(
      res_shrunk,
      alpha = cfg$padj_cutoff,
      main = paste0("MA plot  (apeglm shrinkage)\n",
                      "padj < ", cfg$padj_cutoff, "  |LFC| \u2265 ", cfg$lfc_cutoff),
      colSig = "#E74C3C", colNS = "grey60", ylim = c(-6, 6)
    )
    abline(h = c(-cfg$lfc_cutoff, cfg$lfc_cutoff), lty = 2, col = "steelblue", lwd = 1.2)
  },
  width = 7, height = 5
)

# Volcano plot
ev_df <- as.data.frame(res_shrunk)
ev_df$gene_name <- gene_lookup$gene_name[match(rownames(ev_df), gene_lookup$gene_id)]
ev_df$gene_id <- rownames(ev_df)

# Use gene_name as label; if duplicated, append gene_id to disambiguate
ev_df$label <- ifelse(is.na(ev_df$gene_name), ev_df$gene_id, ev_df$gene_name)
dup_labels <- ev_df$label[duplicated(ev_df$label)]
ev_df$label <- ifelse(
  ev_df$label %in% dup_labels,
  paste0(ev_df$label, "_", ev_df$gene_id),
  ev_df$label
)
rownames(ev_df) <- ev_df$label

top_labels <- head(ev_df$label[order(ev_df$padj)], 20)

save_plot(
  "01_Volcano_plot.pdf",
  {
    print(
      EnhancedVolcano::EnhancedVolcano(
        ev_df,
        lab = rownames(ev_df),
        x = "log2FoldChange",
        y = "padj",
        selectLab = top_labels,
        xlab = bquote(Log[2]~fold~change~"(disease / healthy)"),
        ylab = bquote(-Log[10]~italic(p)[adj]),
        pCutoff = cfg$padj_cutoff,
        FCcutoff = cfg$lfc_cutoff,
        col = c(
          palette$volcano_ns,
          palette$volcano_ns,
          palette$volcano_down,
          palette$volcano_up),
        colAlpha = 0.7,
        pointSize = 2.0,
        labSize = 3.2,
        drawConnectors = TRUE,
        widthConnectors = 0.4,
        title = "Volcano plot — DESeq2 (apeglm shrinkage)",
        subtitle = paste0("Hashimoto's Thyroiditis vs Healthy Controls\n",
                                 "padj < ", cfg$padj_cutoff,
                                 "  |LFC| \u2265 ", cfg$lfc_cutoff),
        caption = paste0("Total DEGs: ", nrow(sig_df),
                                 "  (up: ", sum(sig_df$direction == "up_in_disease"),
                                 "  |  down: ", sum(sig_df$direction == "down_in_disease"), ")")
      )
    )
  },
  width = 9, height = 8
)

# Heatmap (for top 50 significant DEGs)
n_top_heatmap <- min(50, nrow(sig_df))
if (n_top_heatmap > 0) {
  top_gene_ids <- sig_df$gene_id[seq_len(n_top_heatmap)]
  vst_mat <- assay(vst)
  common_ids <- intersect(top_gene_ids, rownames(vst_mat))
  heat_mat <- vst_mat[common_ids, , drop = FALSE]
  
  # Replace gene_id rownames with gene_name where available
  gene_name_map <- setNames(sig_df$gene_name, sig_df$gene_id)
  rownames(heat_mat) <- ifelse(
    is.na(gene_name_map[rownames(heat_mat)]) | gene_name_map[rownames(heat_mat)] == "",
    rownames(heat_mat),
    gene_name_map[rownames(heat_mat)]
  )
  # Z-score scaling
  heat_mat_z <- t(scale(t(heat_mat)))
  
  annot_col_heat <- data.frame(Condition = vst$condition, row.names = colnames(vst))
  annot_colours <- list(Condition = palette$condition)
  
  save_plot(
    "01_Heatmap_50.pdf",
    {
      pheatmap::pheatmap(
        heat_mat_z,
        annotation_col = annot_col_heat,
        annotation_colors = annot_colours,
        color = palette$heatmap,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 7,
        fontsize_col = 9,
        border_color = NA,
        main = paste0("Top ", n_top_heatmap, " significant DEGs\n",
                                   "(z-score of VST counts)"),
        breaks = seq(-2, 2, length.out = 101)
      )
    },
    width = 9, height = 10
  )
} else 
  {message("No significant DEGs to plot in heatmap.")}

# plotCounts (top 6)
n_top_counts <- min(6, nrow(sig_df))
if (n_top_counts > 0) {
  top6_ids <- sig_df$gene_id[seq_len(n_top_counts)]
  top6_names <- sig_df$gene_name[seq_len(n_top_counts)]
  top6_names <- ifelse(is.na(top6_names) | top6_names == "", top6_ids, top6_names)
  count_plots <- lapply(seq_len(n_top_counts), function(i) {
    cnt_df <- DESeq2::plotCounts(dds, gene = top6_ids[i],
      intgroup = "condition",
      returnData = TRUE,
      normalized = TRUE)
    cnt_df$condition <- factor(cnt_df$condition, levels = c("healthy", "disease"))
    
    ggplot(cnt_df, aes(x = condition, y = count, colour = condition)) +
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.4, linewidth = 0.6, colour = "black") +
      scale_colour_manual(values = palette$condition) +
      scale_y_log10(labels = scales::comma) +
      labs(
        title = top6_names[i],
        subtitle = paste0("padj = ", formatC(sig_df$padj[i], digits = 2, format = "e"),
                          "   LFC = ", round(sig_df$log2FoldChange[i], 2)),
        x = NULL,
        y = "Normalised counts (log10)"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey40")
      )
  })
  
  save_plot("01_plotcounts_top6.pdf",
    print(patchwork::wrap_plots(count_plots, ncol = 3)),
    width = 12, height = 8
  )
  
} else {message("No significant DEGs to plot in plotCounts.")}

message("Genes tested: ", nrow(res_shrunk))
message("Significant DEGs: ", nrow(sig_df),
        "(padj < ", cfg$padj_cutoff,
        "|LFC| \u2265 ", cfg$lfc_cutoff, ")")
message("Up in disease: ", sum(sig_df$direction == "up_in_disease"))
message("Down in disease: ", sum(sig_df$direction == "down_in_disease"))

