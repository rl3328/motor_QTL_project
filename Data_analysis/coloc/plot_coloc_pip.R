#!/usr/bin/env Rscript
# Colocalization data extraction and plotting with PIP on y-axis
# Separated into two functions: extract_coloc_data() and plot_coloc_from_data()

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

################################################################################
# FUNCTION 1: Extract and save colocalization data
################################################################################

extract_coloc_data <- function(gene_id,
                                eqtl_data,
                                gwas_file_path,
                                gene_name = NULL,
                                zoom_window = NULL,
                                gwas_threshold = 5e-8,
                                pip_threshold = 0.01,
                                data_dir = "./data") {

  cat(sprintf("Extracting data for gene: %s\n", gene_id))

  # Rename #chr to chr if needed
  if ("#chr" %in% names(eqtl_data)) {
    setnames(eqtl_data, "#chr", "chr")
  }

  # Filter for this gene
  gene_data <- eqtl_data[gene_ID == gene_id]
  cat(sprintf("Found %d variants\n", nrow(gene_data)))

  # Get chromosome
  region_chr <- as.character(unique(gene_data$chr)[1])

  # Determine plotting region
  if (!is.null(zoom_window)) {
    plot_start <- as.integer(zoom_window[1])
    plot_end <- as.integer(zoom_window[2])
  } else {
    sig_variants <- gene_data[PIP > pip_threshold]
    if (nrow(sig_variants) == 0) sig_variants <- gene_data

    plot_start <- as.integer(min(sig_variants$start))
    plot_end <- as.integer(max(sig_variants$end))

    # Add 10% padding
    padding <- as.integer((plot_end - plot_start) * 0.1)
    plot_start <- max(1, plot_start - padding)
    plot_end <- plot_end + padding
  }

  cat(sprintf("Region: chr%s:%d-%d\n", region_chr, plot_start, plot_end))

  # Filter to plotting region
  eqtl_plot_data <- gene_data[start >= plot_start & end <= plot_end]

  # Assign credible set status (priority: 0.95 > 0.70 > 0.50)
  eqtl_plot_data[, cs := NA_character_]
  eqtl_plot_data[cs_coverage_0.5 == 1, cs := "CS_0.50"]
  eqtl_plot_data[cs_coverage_0.7 == 1, cs := "CS_0.70"]
  eqtl_plot_data[cs_coverage_0.95 == 1, cs := "CS_0.95"]

  # Use midpoint for plotting
  eqtl_plot_data[, pos := (start + end) / 2]

  # Count CS variants
  cs_counts <- eqtl_plot_data[!is.na(cs), .N, by = cs]
  bg_variants <- eqtl_plot_data[is.na(cs)]
  cs_variants <- eqtl_plot_data[!is.na(cs)]

  cat(sprintf("eQTL: %d variants (%d in CS)\n", nrow(eqtl_plot_data), nrow(cs_variants)))

  # Query GWAS data - try both chr formats
  gwas_chr_formats <- if (grepl("^chr", region_chr)) {
    c(region_chr, gsub("^chr", "", region_chr))
  } else {
    c(region_chr, paste0("chr", region_chr))
  }

  gwas_data <- NULL
  for (gwas_chr in gwas_chr_formats) {
    gwas_query <- sprintf("tabix %s %s:%d-%d", gwas_file_path, gwas_chr, plot_start, plot_end)
    gwas_data <- tryCatch({
      fread(cmd = gwas_query, header = FALSE,
            col.names = c("chr", "pos", "A1", "A2", "SNP", "effect_allele_frequency",
                          "INFO", "beta", "se", "p", "n_sample"))
    }, error = function(e) NULL)

    if (!is.null(gwas_data) && nrow(gwas_data) > 0) break
  }

  cat(sprintf("GWAS: %d variants\n", nrow(gwas_data)))

  # Calculate -log10(p) for GWAS
  gwas_data[, neg_log_p := -log10(pmax(p, 1e-300))]

  # Save data
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  finemap_data <- list(
    gene_id = gene_id,
    gene_name = gene_name,
    region = list(chr = region_chr, start = plot_start, end = plot_end),
    eqtl_data = eqtl_plot_data,
    gwas_data = gwas_data,
    thresholds = list(gwas_threshold = gwas_threshold, pip_threshold = pip_threshold),
    metadata = list(
      n_variants_total = nrow(eqtl_plot_data),
      n_variants_bg = nrow(bg_variants),
      n_variants_cs = nrow(cs_variants),
      cs_counts = if(nrow(cs_counts) > 0) as.data.frame(cs_counts) else NULL
    )
  )

  data_file <- file.path(data_dir, sprintf("%s_finemap.rds", gene_id))
  saveRDS(finemap_data, data_file)
  cat(sprintf("Saved: %s\n", data_file))

  invisible(finemap_data)
}

################################################################################
# FUNCTION 2: Load saved data and create plots
################################################################################

plot_coloc_from_data <- function(gene_id,
                                  data_dir = "./data",
                                  output_file = NULL,
                                  plot_width = NULL,
                                  plot_height = 12) {

  # Load saved data
  data_file <- file.path(data_dir, sprintf("%s_finemap.rds", gene_id))
  finemap_data <- readRDS(data_file)

  # Extract components
  gene_name <- finemap_data$gene_name
  region_chr <- finemap_data$region$chr
  plot_start <- finemap_data$region$start
  plot_end <- finemap_data$region$end
  eqtl_plot_data <- finemap_data$eqtl_data
  gwas_data <- finemap_data$gwas_data
  gwas_threshold <- finemap_data$thresholds$gwas_threshold
  bg_variants <- eqtl_plot_data[is.na(cs)]
  cs_variants <- eqtl_plot_data[!is.na(cs)]

  cat(sprintf("Plotting %s (%d eQTL, %d GWAS variants)\n",
              gene_id, nrow(eqtl_plot_data), nrow(gwas_data)))


  # Calculate threshold
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
    plot_width <- min(18, max(10, region_mb * 2))
  }

  # Save plot
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p_combined,
           width = plot_width, height = plot_height, dpi = 300, bg = "white")
    cat(sprintf("Saved: %s\n", output_file))
  }

  p_combined
}
