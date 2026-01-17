################################################################################
# eQTL-GWAS Comparison Functions
#
# This file contains functions for extracting and visualizing eQTL and GWAS data.
# Source this file in your analysis scripts or notebooks.
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

################################################################################
# Configuration
################################################################################

# Data paths
# EQTL_DIR <- "/mnt/lustre/home/rl3328/rl3328/motor_qtl/eQTL"
# GWAS_FILE <- "/mnt/lustre/home/rl3328/rl3328/motor_qtl/GWAS/unadjusted_wp_ukb.sumstats_hg38_sorted.tsv.gz"

# Tissue names
TISSUES <- c("Quad", "VM", "SMA")

# Color palette for tissues and GWAS
TISSUE_COLORS <- c(
  "Quad" = "#1D96C3",   # Blue
  "SMA" = "#A24443",    # Red
  "VM" = "#73A58B",     # Green
  "GWAS" = "#35464A"    # Dark gray
)

# Significance thresholds
GWAS_THRESHOLD <- 5e-8
EQTL_THRESHOLD <- 1e-5

################################################################################
# Data Loading Functions
################################################################################

#' Load eQTL data for a specific tissue
#'
#' @param tissue Tissue name (Quad, SMA, or VM)
#' @param chromosome Chromosome number (e.g., 1, 2, ..., 22, 'X', 'Y')
#' @param gene_id Gene ID (e.g., 'ENSG00000101557')
#' @param start_pos Start position for region filter
#' @param end_pos End position for region filter
#' @return data.table with eQTL data
load_eqtl_data <- function(tissue,
                           chromosome = NULL,
                           gene_id = NULL,
                           start_pos = NULL,
                           end_pos = NULL) {

  # Convert to proper types
  if (!is.null(start_pos)) start_pos <- as.integer(start_pos)
  if (!is.null(end_pos)) end_pos <- as.integer(end_pos)

  tissue_dir <- file.path(EQTL_DIR, tissue)

  # Find relevant files
  if (!is.null(chromosome)) {
    chr_str <- paste0("chr", chromosome)
    pattern <- paste0("*", chr_str, "_", chr_str, ".cis_qtl.pairs.tsv.gz")
    files <- list.files(tissue_dir, pattern = pattern, full.names = TRUE)
  } else {
    files <- list.files(tissue_dir, pattern = "*.cis_qtl.pairs.tsv.gz", full.names = TRUE)
  }

  if (length(files) == 0) {
    message(sprintf("No eQTL files found for %s, chromosome %s", tissue, chromosome))
    return(data.table())
  }

  # Load and concatenate data
  dfs <- list()
  for (file in files) {
    tryCatch({
      df <- fread(file, sep = "\t", showProgress = FALSE)

      # Filter by gene if specified
      if (!is.null(gene_id)) {
        df <- df[molecular_trait_id == gene_id]
      }

      # Filter by position if specified
      if (!is.null(start_pos) && !is.null(end_pos)) {
        df <- df[pos >= start_pos & pos <= end_pos]
      }

      if (nrow(df) > 0) {
        df[, tissue := tissue]
        dfs[[length(dfs) + 1]] <- df
      }
    }, error = function(e) {
      message(sprintf("Error loading %s: %s", file, e$message))
    })
  }

  if (length(dfs) == 0) {
    return(data.table())
  }

  result <- rbindlist(dfs, fill = TRUE)

  # Clean chromosome format
  result[, chrom := gsub("chr", "", chrom)]

  return(result)
}


#' Load GWAS data
#'
#' @param chromosome Chromosome number
#' @param start_pos Start position for region filter
#' @param end_pos End position for region filter
#' @return data.table with GWAS data
load_gwas_data <- function(chromosome = NULL,
                          start_pos = NULL,
                          end_pos = NULL) {

  # Convert to proper types
  if (!is.null(start_pos)) start_pos <- as.integer(start_pos)
  if (!is.null(end_pos)) end_pos <- as.integer(end_pos)

  # Check if tabix is available
  use_tabix <- FALSE
  if (!is.null(chromosome) && system("which tabix", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
    use_tabix <- TRUE
  }

  if (use_tabix && !is.null(chromosome)) {
    chr_str <- as.character(chromosome)

    # Use tabix for indexed access
    if (!is.null(start_pos) && !is.null(end_pos)) {
      region <- sprintf("%s:%d-%d", chr_str, start_pos, end_pos)
    } else {
      region <- chr_str
    }

    tryCatch({
      cmd <- sprintf("tabix %s %s", GWAS_FILE, region)
      df <- fread(cmd = cmd, sep = "\t", showProgress = FALSE)

      # Get header
      header_cmd <- sprintf("zcat %s | head -1", GWAS_FILE)
      header <- system(header_cmd, intern = TRUE)
      header <- strsplit(header, "\t")[[1]]

      setnames(df, header)

      # Fix column names if needed
      if ("#chrom" %in% names(df)) {
        setnames(df, "#chrom", "chrom")
      }

    }, error = function(e) {
      message(sprintf("Tabix error: %s", e$message))
      message("Falling back to full file reading...")
      use_tabix <<- FALSE
    })
  }

  if (!use_tabix) {
    # Fallback: read full file with filtering
    if (is.null(chromosome)) {
      message("Loading entire GWAS file... this may take a while")
    }
    df <- fread(GWAS_FILE, sep = "\t", showProgress = FALSE)

    # Fix column names if needed
    if ("#chrom" %in% names(df)) {
      setnames(df, "#chrom", "chrom")
    }

    if (!is.null(chromosome)) {
      chr_str <- as.character(chromosome)
      df <- df[chrom == chr_str]

      if (!is.null(start_pos) && !is.null(end_pos)) {
        df <- df[pos >= start_pos & pos <= end_pos]
      }
    }
  }

  # Clean chromosome format
  df[, chrom := gsub("chr", "", chrom)]

  return(df)
}


################################################################################
# Data Extraction Function
################################################################################

#' Extract GWAS and eQTL data for a genomic region
#'
#' @param chromosome Chromosome number
#' @param center_pos Center position for region
#' @param window_size Size of window in bp (default: 500kb)
#' @param gene_id Gene ID to filter (optional)
#' @param tissues List of tissues to include
#' @return List containing GWAS data, eQTL data by tissue, and metadata
extract_region_data <- function(chromosome,
                                center_pos,
                                window_size = 500000,
                                gene_id = NULL,
                                tissues = c("Quad", "VM", "SMA")) {

  # Ensure proper types
  chromosome <- as.character(chromosome)
  center_pos <- as.integer(center_pos)
  window_size <- as.integer(window_size)

  # Calculate region boundaries (convert to integers)
  start_pos <- as.integer(max(0, center_pos - window_size / 2))
  end_pos <- as.integer(center_pos + window_size / 2)

  message(sprintf("Extracting data for Chr%s: %s - %s",
                 chromosome,
                 format(start_pos, big.mark = ","),
                 format(end_pos, big.mark = ",")))

  # Load GWAS data
  message("  Loading GWAS data...")
  gwas_df <- load_gwas_data(chromosome, start_pos, end_pos)
  message(sprintf("    Found %s variants", format(nrow(gwas_df), big.mark = ",")))

  # Load eQTL data for each tissue
  eqtl_list <- list()
  for (tissue in tissues) {
    message(sprintf("  Loading %s eQTL data...", tissue))
    df <- load_eqtl_data(tissue, chromosome, gene_id, start_pos, end_pos)
    if (nrow(df) > 0) {
      eqtl_list[[tissue]] <- df
      message(sprintf("    Found %s variant-gene pairs", format(nrow(df), big.mark = ",")))
    } else {
      message(sprintf("    No data found"))
    }
  }

  # Create metadata
  metadata <- list(
    chromosome = chromosome,
    start_pos = start_pos,
    end_pos = end_pos,
    center_pos = center_pos,
    window_size = window_size,
    gene_id = gene_id,
    tissues = tissues,
    extraction_time = Sys.time()
  )

  # Return combined data
  result <- list(
    gwas = gwas_df,
    eqtl = eqtl_list,
    metadata = metadata
  )

  message("\nData extraction complete!")

  return(result)
}


#' Save extracted data to RDS file
#'
#' @param data_list Data list from extract_region_data()
#' @param file_path Path to save RDS file
save_extracted_data <- function(data_list, file_path) {
  # Create directory if it doesn't exist
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)

  # Save to RDS
  saveRDS(data_list, file_path)
  message(sprintf("Data saved to: %s", file_path))

  # Print summary
  message(sprintf("  GWAS variants: %s", format(nrow(data_list$gwas), big.mark = ",")))
  for (tissue in names(data_list$eqtl)) {
    message(sprintf("  %s eQTL: %s pairs", tissue, format(nrow(data_list$eqtl[[tissue]]), big.mark = ",")))
  }
}


################################################################################
# Visualization Function
################################################################################

#' Plot multi-panel figure from extracted data
#'
#' @param data_list Data list from extract_region_data() or loaded from RDS
#' @param panel_heights Relative heights for panels (GWAS, then tissues)
#' @param total_width Total figure width in inches
#' @param total_height Total figure height in inches
#' @param show_thresholds Whether to show significance thresholds
#' @param plot_title Custom plot title
#' @param save_path Path to save figure
#' @param gene_name_map Named vector mapping ENSG IDs to gene names (e.g., c("ENSG00000141458" = "NPC1"))
#' @return Combined plot object (patchwork)
plot_extracted_data <- function(data_list,
                                panel_heights = c(1, 0.8, 0.8, 0.8),
                                total_width = 12,
                                total_height = 10,
                                show_thresholds = TRUE,
                                plot_title = NULL,
                                save_path = NULL,
                                gene_name_map = NULL) {

  # Load patchwork for combining plots
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }
  library(patchwork)

  # Extract metadata
  meta <- data_list$metadata
  gwas_df <- data_list$gwas
  eqtl_list <- data_list$eqtl

  # Create title with gene name
  if (is.null(plot_title)) {
    if (!is.null(meta$gene_id)) {
      # Use gene name if available in mapping, otherwise use ENSG ID
      if (!is.null(gene_name_map) && meta$gene_id %in% names(gene_name_map)) {
        gene_display <- gene_name_map[meta$gene_id]
      } else {
        gene_display <- meta$gene_id
      }
      title_parts <- paste("Colocalization of GWAS and eQTL for ", gene_display)
    }
    plot_title <- title_parts
  }

  # Common theme for all panels
  common_theme <- theme_bw(base_size = 11) +
    theme(
      axis.text = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.margin = margin(5, 10, 5, 10)
    )

  plots <- list()

  # --- GWAS Panel ---
  if (nrow(gwas_df) > 0) {
    gwas_df[, neg_log_p := -log10(pmax(p, 1e-300))]

    p_gwas <- ggplot(gwas_df, aes(x = pos, y = neg_log_p)) +
      geom_point(
        color = TISSUE_COLORS["GWAS"],
        size = 1.2,
        alpha = 0.6
      ) +
      labs(
        y = expression(bold(-log[10](italic(p)))),
        title = "GWAS (Gait Speed)"
      ) +
      scale_x_continuous(labels = comma) +
      common_theme +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0, color = TISSUE_COLORS["GWAS"])
      )

    if (show_thresholds) {
      p_gwas <- p_gwas +
        geom_hline(
          yintercept = -log10(GWAS_THRESHOLD),
          color = "#8B0000",
          linetype = "dashed",
          linewidth = 0.6
        )
    }

    plots[["GWAS"]] <- p_gwas
  }

  # --- eQTL Panels (one per tissue) ---
  for (tissue in meta$tissues) {
    if (tissue %in% names(eqtl_list)) {
      tissue_data <- eqtl_list[[tissue]]
      tissue_data[, neg_log_p := -log10(pmax(pvalue, 1e-300))]

      p_eqtl <- ggplot(tissue_data, aes(x = pos, y = neg_log_p)) +
        geom_point(
          color = TISSUE_COLORS[tissue],
          size = 2,
          alpha = 0.7
        ) +
        labs(
          y = expression(bold(-log[10](italic(p)))),
          title = paste0(tissue, " eQTL")
        ) +
        scale_x_continuous(labels = comma) +
        common_theme +
        theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0, color = TISSUE_COLORS[tissue])
        )

      if (show_thresholds) {
        p_eqtl <- p_eqtl +
          geom_hline(
            yintercept = -log10(EQTL_THRESHOLD),
            color = "#CD5C5C",
            linetype = "dotted",
            linewidth = 0.6
          )
      }

      plots[[tissue]] <- p_eqtl
    }
  }

  # Combine plots with patchwork
  # Add x-axis label only to bottom plot
  n_plots <- length(plots)
  plot_names <- names(plots)

  # Update bottom plot to have x-axis label
  plots[[n_plots]] <- plots[[n_plots]] +
    labs(x = sprintf("Chromosome %s Position (bp)", meta$chromosome)) +
    theme(axis.title.x = element_text(size = 10, face = "bold"))

  # Combine with specified heights
  combined_plot <- wrap_plots(plots, ncol = 1, heights = panel_heights[1:n_plots]) +
    plot_annotation(
      title = plot_title,
      theme = theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    )

  # Save if requested
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = combined_plot,
      width = total_width,
      height = total_height,
      dpi = 300,
      bg = "white"
    )
    message(sprintf("\nFigure saved to: %s", save_path))
    message(sprintf("Dimensions: %.1f x %.1f inches", total_width, total_height))
  }

  return(combined_plot)
}


################################################################################
# Helper Functions
################################################################################

#' Find all genes with eQTL data on a specific chromosome
#'
#' @param chromosome Chromosome number
#' @param tissues Tissues to search (default: just Quad for speed)
#' @return data.frame with gene IDs and their genomic positions
find_genes_on_chromosome <- function(chromosome, tissues = c("Quad")) {

  dfs <- list()
  for (tissue in tissues) {
    df <- load_eqtl_data(tissue, chromosome)
    if (nrow(df) > 0) {
      dfs[[length(dfs) + 1]] <- df
    }
  }

  if (length(dfs) == 0) {
    return(data.frame())
  }

  combined <- rbindlist(dfs)

  # Get unique genes with position info
  gene_info <- combined[, .(
    min_pos = min(pos),
    max_pos = max(pos),
    n_variants = .N
  ), by = molecular_trait_id]

  gene_info[, center_pos := as.integer((min_pos + max_pos) / 2)]
  setnames(gene_info, "molecular_trait_id", "gene_id")

  gene_info <- gene_info[order(min_pos)]

  return(as.data.frame(gene_info))
}


#' Check which tissues have data for a specific gene
#'
#' @param gene_id Gene ID
#' @param chromosome Chromosome number
check_gene_tissues <- function(gene_id, chromosome) {
  cat(sprintf("Checking which tissues have data for %s on chr%s...\n\n", gene_id, chromosome))

  for (tissue in c("Quad", "VM", "SMA")) {
    df <- load_eqtl_data(tissue, chromosome, gene_id = gene_id)
    if (nrow(df) > 0) {
      cat(sprintf("✓ %s: %s variants (p-value range: %.2e - %.2e)\n",
                 tissue,
                 format(nrow(df), big.mark = ","),
                 min(df$pvalue),
                 max(df$pvalue)))
    } else {
      cat(sprintf("✗ %s: NO DATA\n", tissue))
    }
  }
}


message("eQTL-GWAS comparison functions loaded successfully!")
