#!/usr/bin/env Rscript
# Colocalization plotting with PIP on y-axis and larger fonts

library(data.table)
library(tidyverse)
library(patchwork)

# Colors
gwas_color <- "#A24443"  # Red for GWAS
eqtl_base_color <- "#CCCCCC"  # Grey for background variants

# Color palette for credible sets (colorblind-friendly)
cs_colors <- c(
  "CS_0.95" = "#E69F00",  # Orange - high confidence CS
  "CS_0.70" = "#56B4E9",  # Sky blue - medium confidence CS
  "CS_0.50" = "#009E73"   # Green - lower confidence CS
)

plot_colocalization_pip <- function(gene_id,
                                     eqtl_bed_file,
                                     gwas_file,
                                     gene_name = NULL,
                                     zoom_window = NULL,
                                     gwas_threshold = 5e-8,
                                     pip_threshold = 0.01,
                                     output_file = NULL,
                                     plot_width = NULL,
                                     plot_height = 12) {

  # Read eQTL data for the gene
  cat(sprintf("Reading eQTL data for gene: %s\n", gene_id))

  # Read the full bed file and filter by gene
  eqtl_data <- fread(cmd = sprintf("zcat %s", eqtl_bed_file), header = TRUE)
  setnames(eqtl_data, "#chr", "chr")

  # Filter for this gene
  gene_data <- eqtl_data[gene_ID == gene_id]

  if (nrow(gene_data) == 0) {
    stop(sprintf("No data found for gene %s", gene_id))
  }

  cat(sprintf("Found %d variants for gene %s\n", nrow(gene_data), gene_id))

  # Get chromosome (should be consistent)
  region_chr <- unique(gene_data$chr)
  if (length(region_chr) > 1) {
    warning("Multiple chromosomes found, using first one")
    region_chr <- region_chr[1]
  }
  region_chr <- as.character(region_chr)

  # Determine plotting region
  if (!is.null(zoom_window)) {
    plot_start <- as.integer(zoom_window[1])
    plot_end <- as.integer(zoom_window[2])
    cat(sprintf("Using zoom window: chr%s:%d-%d\n", region_chr, plot_start, plot_end))
  } else {
    # Use range of variants with PIP above threshold
    sig_variants <- gene_data[PIP > pip_threshold]

    if (nrow(sig_variants) == 0) {
      # If no variants above threshold, use all
      sig_variants <- gene_data
    }

    plot_start <- as.integer(min(sig_variants$start))
    plot_end <- as.integer(max(sig_variants$end))

    # Add 10% padding
    padding <- as.integer((plot_end - plot_start) * 0.1)
    plot_start <- max(1, plot_start - padding)
    plot_end <- plot_end + padding

    cat(sprintf("Using variant range: chr%s:%d-%d\n", region_chr, plot_start, plot_end))
  }

  # Filter to plotting region
  eqtl_plot_data <- gene_data[start >= plot_start & end <= plot_end]

  # Assign credible set status (priority: 0.95 > 0.70 > 0.50)
  eqtl_plot_data[, cs := NA_character_]
  eqtl_plot_data[cs_coverage_0.5 == 1, cs := "CS_0.50"]
  eqtl_plot_data[cs_coverage_0.7 == 1, cs := "CS_0.70"]
  eqtl_plot_data[cs_coverage_0.95 == 1, cs := "CS_0.95"]

  # Use midpoint for plotting
  eqtl_plot_data[, pos := (start + end) / 2]

  cat(sprintf("\nPlotting %d total variants\n", nrow(eqtl_plot_data)))

  # Count CS variants
  cs_counts <- eqtl_plot_data[!is.na(cs), .N, by = cs]
  if (nrow(cs_counts) > 0) {
    cat("Credible set variants:\n")
    print(cs_counts)
  } else {
    cat("No credible set variants found\n")
  }

  # Separate background and CS variants
  bg_variants <- eqtl_plot_data[is.na(cs)]
  cs_variants <- eqtl_plot_data[!is.na(cs)]

  cat(sprintf("  Background: %d variants\n", nrow(bg_variants)))
  cat(sprintf("  In credible sets: %d variants\n", nrow(cs_variants)))

  # Read GWAS data
  cat("\nReading GWAS data...\n")
  gwas_query <- sprintf("tabix %s %s:%d-%d", gwas_file, region_chr, plot_start, plot_end)
  cat("Query:", gwas_query, "\n")

  gwas_data <- fread(cmd = gwas_query, header = FALSE,
                     col.names = c("chr", "pos", "A1", "A2", "SNP", "effect_allele_frequency",
                                   "INFO", "beta", "se", "p", "n_sample"))

  if (nrow(gwas_data) == 0) {
    stop(sprintf("No GWAS data found in region chr%s:%d-%d", region_chr, plot_start, plot_end))
  }

  cat(sprintf("GWAS: %d variants in region\n", nrow(gwas_data)))

  # Calculate -log10(p) for GWAS
  gwas_data[, neg_log_p := -log10(pmax(p, 1e-300))]
  gwas_threshold_y <- -log10(gwas_threshold)

  # LARGE FONT SIZES (increased from original)
  base_font_size <- 18  # Was 12
  axis_text_size <- 16  # Was 10
  axis_title_size <- 20  # Was 11
  plot_title_size <- 22  # Was 14
  legend_title_size <- 18  # Was 10
  legend_text_size <- 16  # Was 9
  threshold_text_size <- 5  # Was 3

  # GWAS plot (top panel)
  p_gwas <- ggplot(gwas_data, aes(x = pos / 1e6, y = neg_log_p)) +
    geom_point(alpha = 0.6, size = 2.5, color = gwas_color) +  # Larger points
    geom_hline(yintercept = gwas_threshold_y, linetype = "dashed",
               color = "#8B0000", linewidth = 1.2) +  # Thicker line
    annotate("text", x = plot_start / 1e6, y = gwas_threshold_y + 0.5,
             label = sprintf("p = %.0e", gwas_threshold),
             hjust = 0, vjust = 0, size = threshold_text_size,
             color = "#8B0000", fontface = "bold") +
    labs(title = sprintf("GWAS Signal: chr%s:%.2f-%.2f Mb",
                         region_chr, plot_start/1e6, plot_end/1e6),
         x = NULL,
         y = expression(bold(-log[10] * " " * italic(p)))) +
    theme_bw(base_size = base_font_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = plot_title_size),
      axis.title.y = element_text(face = "bold", size = axis_title_size),
      axis.text = element_text(size = axis_text_size, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 15, r = 15, b = 5, l = 15)
    )

  # eQTL plot (bottom panel) - PIP on y-axis
  p_eqtl <- ggplot() +
    # Background variants (grey)
    geom_segment(data = bg_variants,
                 aes(x = pos / 1e6, xend = pos / 1e6, y = 0, yend = PIP),
                 alpha = 0.3, linewidth = 0.8, color = eqtl_base_color) +
    geom_point(data = bg_variants,
               aes(x = pos / 1e6, y = PIP),
               alpha = 0.3, size = 1.5, color = eqtl_base_color)

  # Add credible set variants if present
  if (nrow(cs_variants) > 0) {
    # Order CS by confidence for plotting
    cs_variants[, cs := factor(cs, levels = c("CS_0.95", "CS_0.70", "CS_0.50"))]

    p_eqtl <- p_eqtl +
      geom_segment(data = cs_variants,
                   aes(x = pos / 1e6, xend = pos / 1e6, y = 0, yend = PIP, color = cs),
                   alpha = 0.9, linewidth = 1.5) +  # Thicker lines
      geom_point(data = cs_variants,
                 aes(x = pos / 1e6, y = PIP, color = cs),
                 alpha = 0.9, size = 3.5) +  # Larger points
      scale_color_manual(values = cs_colors,
                         name = "Credible Set",
                         labels = c("CS 95%", "CS 70%", "CS 50%"),
                         na.translate = FALSE)
  }

  # Finish eQTL plot
  # Create title with gene name if provided
  eqtl_title <- if (!is.null(gene_name)) {
    sprintf("%s (%s) - DLPFC eQTL", gene_name, gene_id)
  } else {
    sprintf("%s (DLPFC eQTL)", gene_id)
  }

  p_eqtl <- p_eqtl +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(title = eqtl_title,
         x = sprintf("Position on Chromosome %s (Mb)", region_chr),
         y = "PIP (Posterior Inclusion Probability)") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    theme_bw(base_size = base_font_size) +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = plot_title_size),
      axis.title = element_text(face = "bold", size = axis_title_size),
      axis.text = element_text(size = axis_text_size, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = legend_title_size),
      legend.text = element_text(size = legend_text_size, face = "bold"),
      legend.key.size = unit(1.5, "cm"),  # Larger legend keys
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 15, b = 15, l = 15)
    )

  # Combine plots
  p_combined <- p_gwas / p_eqtl + plot_layout(heights = c(1, 1))

  # Calculate plot width
  if (is.null(plot_width)) {
    region_mb <- (plot_end - plot_start) / 1e6
    plot_width <- min(18, max(10, region_mb * 2))  # Slightly wider
  }

  # Save plot
  if (!is.null(output_file)) {
    cat(sprintf("\nSaving plot to: %s\n", output_file))
    ggsave(output_file, plot = p_combined,
           width = plot_width, height = plot_height, dpi = 300, bg = "white")
  }

  cat("\nPlot complete!\n")
  return(p_combined)
}


# Test if run as script
if (!interactive()) {
  cat("==========================================\n")
  cat("Testing PIP-based colocalization plot\n")
  cat("==========================================\n\n")

  p1 <- plot_colocalization_pip(
    gene_id = "ENSG00000188338",
    eqtl_bed_file = "DLPFC_DeJager_eQTL.exported.toploci.bed.gz",
    gwas_file = "GWAS/unadjusted_wp_ukb.sumstats_hg38_sorted.tsv.gz",
    output_file = "colocalization_ENSG00000188338_PIP.png"
  )

  cat("\n==========================================\n")
  cat("Test completed successfully!\n")
  cat("==========================================\n")
}
