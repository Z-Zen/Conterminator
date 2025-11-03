#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_gc_content.R <input_gc_coverage.tsv> <output_prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# Read data with all columns as character first
cat("Reading data from", input_file, "\n")
raw_data <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE, 
                        colClasses = "character")

cat("Raw columns:", paste(colnames(raw_data), collapse = ", "), "\n")
cat("Number of raw columns:", ncol(raw_data), "\n")

# The paste command creates columns:
# 1:chrom, 2:start, 3:end, 4:pct_at, 5:pct_gc, 6:chrom2, 7:start2, 8:end2, 9:mean_coverage_value, 10:header_text
# We need columns 4, 5, and 9

if (ncol(raw_data) >= 9) {
  data <- raw_data %>%
    select(
      chrom = 1,
      start = 2,
      end = 3,
      pct_at = 4,
      pct_gc = 5,
      mean_coverage = 9  # The actual coverage value from bedtools
    )
} else {
  stop("Unexpected number of columns in input file. Expected at least 9 columns, got ", ncol(raw_data))
}

# Convert to numeric and clean
data <- data %>%
  mutate(
    chrom = as.character(chrom),
    start = as.numeric(start),
    end = as.numeric(end),
    pct_at = as.numeric(pct_at),
    pct_gc = as.numeric(pct_gc),
    mean_coverage = as.numeric(mean_coverage)
  )

# Check for data validity before filtering
cat("Total windows:", nrow(data), "\n")
cat("Windows with non-zero GC:", sum(data$pct_gc > 0, na.rm = TRUE), "\n")
cat("Windows with valid coverage:", sum(!is.na(data$mean_coverage)), "\n")

# Now filter
data <- data %>%
  filter(
    !is.na(start),
    !is.na(end),
    !is.na(mean_coverage),
    !is.na(pct_gc)
  )

# Check if we have any non-zero GC content
if (all(data$pct_gc == 0, na.rm = TRUE)) {
  cat("\n")
  cat("ERROR: All GC content values are 0!\n")
  cat("This means your reference genome has no actual sequence data in these regions.\n")
  cat("Possible causes:\n")
  cat("  1. The reference genome is masked (all N's)\n")
  cat("  2. The windows are in telomeric/centromeric regions with no sequence\n")
  cat("  3. The reference file is corrupted or incorrect\n")
  cat("\nPlease:\n")
  cat("  1. Use an unmasked reference genome (not *.masked.fa)\n")
  cat("  2. Check: samtools faidx <ref.fa> chr1:1000000-1000500\n")
  cat("  3. Verify you see ACGT bases, not all N's\n")
  quit(status = 1)
}

# Filter out zero GC regions for plotting
data <- data %>%
  filter(pct_gc > 0)

# Check if we have enough data
if (nrow(data) == 0) {
  cat("\n")
  cat("ERROR: No windows with GC content > 0 found!\n")
  quit(status = 1)
}

cat("Processed", nrow(data), "valid windows\n")
cat("GC range:", round(min(data$pct_gc, na.rm = TRUE), 4), "-", 
    round(max(data$pct_gc, na.rm = TRUE), 4), "\n")
cat("Coverage range:", round(min(data$mean_coverage, na.rm = TRUE), 2), "-", 
    round(max(data$mean_coverage, na.rm = TRUE), 2), "\n")

# Subsample for faster plotting if dataset is large
max_points_for_plotting <- 100000
if (nrow(data) > max_points_for_plotting) {
  cat("Subsampling to", max_points_for_plotting, "windows for faster plotting...\n")
  set.seed(42)
  data_plot <- data %>% sample_n(max_points_for_plotting)
  # Keep full data for summary stats
  data_full <- data
} else {
  data_plot <- data
  data_full <- data
}

# Set theme
theme_set(theme_bw(base_size = 12))

# Plot 1: GC content distribution with density line
p1 <- ggplot(data_plot, aes(x = pct_gc * 100)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", 
                 color = "black", alpha = 0.7) +
  geom_density(color = "red", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = median(data_plot$pct_gc) * 100, 
             linetype = "dashed", color = "darkred", linewidth = 0.5, alpha = 0.7) +
  annotate("text", x = median(data_plot$pct_gc) * 100, 
           y = max(density(data_plot$pct_gc * 100)$y) * 0.9,
           label = paste0("Median: ", round(median(data_plot$pct_gc) * 100, 1), "%"),
           hjust = -0.2, color = "darkred", size = 3) +
  labs(
    title = "GC Content Distribution",
    x = "GC Content (%)",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 2: Coverage vs GC content (use hexbin for speed)
# Filter out zero coverage for better visualization
data_plot_nonzero <- data_plot %>% filter(mean_coverage > 0)

if (nrow(data_plot_nonzero) > 100) {
  # Subsample even more for the smooth curve to avoid memory issues
  set.seed(123)
  data_for_smooth <- data_plot_nonzero %>% sample_n(min(5000, nrow(data_plot_nonzero)))
  
  p2 <- ggplot(data_plot_nonzero, aes(x = pct_gc * 100, y = mean_coverage)) +
    geom_bin2d(bins = 50) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") +
    scale_y_log10(limits = c(max(0.01, min(data_plot_nonzero$mean_coverage)), 
                             max(data_plot_nonzero$mean_coverage))) +
    geom_smooth(data = data_for_smooth, method = "loess", color = "red", 
                se = FALSE, linewidth = 1, formula = y ~ x) +
    labs(
      title = "Coverage vs GC Content",
      x = "GC Content (%)",
      y = "Mean Coverage (log scale)",
      fill = "Count (log10)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
} else {
  # Fallback if not enough non-zero data
  p2 <- ggplot(data_plot, aes(x = pct_gc * 100, y = mean_coverage)) +
    geom_point(alpha = 0.3) +
    labs(
      title = "Coverage vs GC Content",
      x = "GC Content (%)",
      y = "Mean Coverage",
      fill = "Count"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Plot 3: Coverage distribution
p3 <- ggplot(data_plot, aes(x = mean_coverage + 0.1)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black", alpha = 0.7) +
  scale_x_log10() +
  labs(
    title = "Coverage Distribution",
    x = "Mean Coverage (log10)",
    y = "Number of Windows"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 4: GC content density by coverage quantiles
# Handle various edge cases with coverage distributions
create_plot4 <- function(data_to_plot) {
  tryCatch({
    # Filter to windows with coverage > 0 for meaningful quartiles
    data_with_coverage <- data_to_plot %>% filter(mean_coverage > 0)
    
    if (nrow(data_with_coverage) > 100) {
      # Calculate quantile breaks
      quant_breaks <- unique(quantile(data_with_coverage$mean_coverage, 
                                       probs = seq(0, 1, 0.25)))
      
      # Check if we have enough unique breaks for quartiles
      if (length(quant_breaks) >= 5) {
        # Add small jitter to ensure unique breaks if needed
        if (length(unique(quant_breaks)) < length(quant_breaks)) {
          quant_breaks <- quant_breaks + seq_along(quant_breaks) * .Machine$double.eps
        }
        
        data_quantiles <- data_with_coverage %>%
          mutate(coverage_group = cut(mean_coverage, 
                                       breaks = quant_breaks,
                                       labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                                       include.lowest = TRUE))
        
        # Calculate medians for each group
        group_medians <- data_quantiles %>%
          group_by(coverage_group) %>%
          summarise(median_gc = median(pct_gc * 100), .groups = "drop")
        
        ggplot(data_quantiles, aes(x = pct_gc * 100, fill = coverage_group)) +
          geom_density(alpha = 0.5) +
          geom_vline(data = group_medians, aes(xintercept = median_gc, color = coverage_group),
                     linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
          labs(
            title = "GC Content by Coverage Quartiles",
            subtitle = "Dashed lines show median GC% for each quartile",
            x = "GC Content (%)",
            y = "Density",
            fill = "Coverage"
          ) +
          scale_fill_brewer(palette = "Set1") +
          scale_color_brewer(palette = "Set1") +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 9))
      } else {
        # Fall back to median split
        median_cov <- median(data_with_coverage$mean_coverage)
        data_split <- data_with_coverage %>%
          mutate(coverage_group = ifelse(mean_coverage <= median_cov, 
                                          "Below Median", "Above Median"))
        
        # Calculate medians for each group
        group_medians <- data_split %>%
          group_by(coverage_group) %>%
          summarise(median_gc = median(pct_gc * 100), .groups = "drop")
        
        ggplot(data_split, aes(x = pct_gc * 100, fill = coverage_group)) +
          geom_density(alpha = 0.5) +
          geom_vline(data = group_medians, aes(xintercept = median_gc, color = coverage_group),
                     linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
          labs(
            title = "GC Content by Coverage Groups",
            subtitle = "Dashed lines show median GC% for each group",
            x = "GC Content (%)",
            y = "Density",
            fill = "Coverage"
          ) +
          scale_fill_manual(values = c("Below Median" = "#E41A1C", 
                                        "Above Median" = "#377EB8")) +
          scale_color_manual(values = c("Below Median" = "#E41A1C", 
                                         "Above Median" = "#377EB8")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 9))
      }
    } else {
      # Not enough data, use simple bin2d
      ggplot(data_to_plot, aes(x = pct_gc * 100, y = mean_coverage + 0.1)) +
        geom_bin2d(bins = 30) +
        scale_fill_viridis_c(trans = "log10") +
        scale_y_log10() +
        labs(
          title = "GC Content vs Coverage",
          x = "GC Content (%)",
          y = "Mean Coverage (log10)",
          fill = "Count"
        ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }
  }, error = function(e) {
    # If all else fails, make a simple bin2d plot
    cat("Note: Using bin2d plot due to coverage distribution\n")
    ggplot(data_to_plot, aes(x = pct_gc * 100, y = mean_coverage + 0.1)) +
      geom_bin2d(bins = 30) +
      scale_fill_viridis_c(trans = "log10") +
      scale_y_log10() +
      labs(
        title = "GC Content vs Coverage",
        x = "GC Content (%)",
        y = "Mean Coverage (log10)",
        fill = "Count"
      ) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })
}

p4 <- create_plot4(data_plot)

# Combine plots
cat("Generating combined plot...\n")
combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# Save combined plot
output_file <- paste0(output_prefix, "_gc_plots.png")
ggsave(output_file, combined, width = 12, height = 10, dpi = 300)
cat("Saved combined plot to", output_file, "\n")

# Save individual plots
ggsave(paste0(output_prefix, "_gc_distribution.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(paste0(output_prefix, "_coverage_vs_gc.png"), p2, width = 7, height = 5, dpi = 300)
ggsave(paste0(output_prefix, "_coverage_distribution.png"), p3, width = 6, height = 5, dpi = 300)
ggsave(paste0(output_prefix, "_gc_by_coverage_quartiles.png"), p4, width = 7, height = 5, dpi = 300)

# Generate summary statistics (use full dataset)
summary_stats <- data_full %>%
  summarise(
    n_windows = n(),
    mean_gc = mean(pct_gc, na.rm = TRUE),
    median_gc = median(pct_gc, na.rm = TRUE),
    sd_gc = sd(pct_gc, na.rm = TRUE),
    mean_coverage = mean(mean_coverage, na.rm = TRUE),
    median_coverage = median(mean_coverage, na.rm = TRUE)
  )

# Save summary statistics
summary_file <- paste0(output_prefix, "_gc_summary.txt")
write.table(summary_stats, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved summary statistics to", summary_file, "\n")

cat("\nDone!\n")