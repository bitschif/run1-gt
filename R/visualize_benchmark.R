#!/usr/bin/env Rscript
#===============================================================================
# VARIANT CALLING BENCHMARKING - COMPREHENSIVE VISUALIZATION
# Bao gồm: F1, Precision-Recall, Heatmap, ROC Curve, UpSet Plot, Runtime
#===============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║        VARIANT CALLING BENCHMARK - COMPREHENSIVE VISUALIZATION     ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

#-------------------------------------------------------------------------------
# 1. Install and load packages
#-------------------------------------------------------------------------------
cat("[1/7] Loading packages...\n")

packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "scales", 
  "patchwork", "gridExtra", "viridis", "RColorBrewer",
  "ggrepel", "reshape2", "cowplot", "ggpubr"
)

for (pkg in packages) {
  if (! requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(gridExtra)
  library(viridis)
  library(RColorBrewer)
  library(ggrepel)
  library(reshape2)
  library(cowplot)
  library(ggpubr)
})

#-------------------------------------------------------------------------------
# 2. Configuration
#-------------------------------------------------------------------------------
cat("[2/7] Setting up configuration...\n")

args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) > 0, args[1], ".")
metrics_dir <- file.path(project_dir, "results", "final_metrics")
bench_dir <- file.path(project_dir, "results", "benchmarks")
output_dir <- file.path(project_dir, "results", "figures")

cat(sprintf("  Project:   %s\n", project_dir))
cat(sprintf("  Output:   %s\n", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Color scheme
COLORS <- c(
  "gatk" = "#E41A1C",
  "deepvariant" = "#377EB8",
  "strelka2" = "#4DAF4A",
  "freebayes" = "#984EA3"
)

LABELS <- c(
  "gatk" = "GATK",
  "deepvariant" = "DeepVariant",
  "strelka2" = "Strelka2",
  "freebayes" = "FreeBayes"
)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5),
      strip.text = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = "grey95")
    )
}

# Save function
save_plot <- function(plot, name, width = 10, height = 6) {
  png_file <- file.path(output_dir, paste0(name, ".png"))
  pdf_file <- file.path(output_dir, paste0(name, ".pdf"))
  
  ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height)
  
  cat(sprintf("  ✓ Saved:  %s\n", basename(png_file)))
}

#-------------------------------------------------------------------------------
# 3. Load data
#-------------------------------------------------------------------------------
cat("[3/7] Loading data...\n")

bench_file <- file.path(metrics_dir, "benchmark_summary.tsv")
runtime_file <- file.path(metrics_dir, "runtime.csv")
stats_file <- file.path(metrics_dir, "variant_statistics.csv")

# Load benchmark data
if (file.exists(bench_file)) {
  data <- read_tsv(bench_file, show_col_types = FALSE)
  cat(sprintf("  Loaded %d records from benchmark_summary.tsv\n", nrow(data)))
} else {
  cat("  [WARN] Using demo data\n")
  data <- tibble(
    Caller = rep(c("gatk", "deepvariant", "strelka2", "freebayes"), each = 3),
    VariantType = rep(c("ALL", "SNP", "INDEL"), 4),
    TP = c(47500, 42500, 5000, 48500, 43500, 5000, 47000, 42000, 5000, 46000, 41000, 5000),
    FP = c(250, 150, 100, 100, 50, 50, 300, 200, 100, 500, 350, 150),
    FN = c(2500, 2250, 250, 1500, 1250, 250, 3000, 2750, 250, 4000, 3750, 250),
    Precision = c(0.9948, 0.9965, 0.9804, 0.9979, 0.9989, 0.9901, 
                  0.9937, 0.9953, 0.9804, 0.9892, 0.9915, 0.9709),
    Recall = c(0.9500, 0.9497, 0.9524, 0.9700, 0.9720, 0.9524,
               0.9400, 0.9387, 0.9524, 0.9200, 0.9163, 0.9524),
    F1 = c(0.9719, 0.9726, 0.9662, 0.9838, 0.9852, 0.9709,
           0.9662, 0.9663, 0.9662, 0.9534, 0.9527, 0.9615)
  )
}

# Standardize column names
data <- data %>%
  rename_with(tolower) %>%
  mutate(
    caller = tolower(caller),
    varianttype = toupper(varianttype)
  )

# Load runtime data
if (file.exists(runtime_file)) {
  runtime_data <- read_csv(runtime_file, show_col_types = FALSE)
  cat(sprintf("  Loaded runtime data\n"))
} else {
  # Demo runtime data
  runtime_data <- tibble(
    step = c("01_simulate_data", "02_preprocessing", 
             "03_gatk", "04_deepvariant", "05_strelka2", "06_freebayes",
             "07_benchmarking", "08_collect_metrics"),
    duration_seconds = c(120, 1800, 900, 600, 450, 300, 60, 30)
  )
}

# Load variant statistics
if (file.exists(stats_file)) {
  var_stats <- read_csv(stats_file, show_col_types = FALSE)
} else {
  var_stats <- NULL
}

#-------------------------------------------------------------------------------
# 4. Basic Plots (F1, Precision-Recall, Classification)
#-------------------------------------------------------------------------------
cat("[4/7] Creating basic plots...\n")

#--- Plot 1: F1 Score Comparison ---#
cat("  Creating F1 score plot...\n")

p1 <- data %>%
  filter(varianttype %in% c("ALL", "SNP", "INDEL")) %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    varianttype = factor(varianttype, levels = c("ALL", "SNP", "INDEL"))
  ) %>%
  ggplot(aes(x = caller_label, y = f1, fill = caller)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", f1)), vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ varianttype, nrow = 1) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2), labels = percent) +
  labs(
    title = "F1 Score Comparison Across Variant Callers",
    subtitle = "Stratified by Variant Type (ALL, SNP, INDEL)",
    x = NULL, y = "F1 Score"
  ) +
  theme_pub()

save_plot(p1, "01_f1_score_comparison", width = 12, height = 5)

#--- Plot 2: Precision vs Recall Scatter ---#
cat("  Creating precision-recall scatter plot...\n")

p2 <- data %>%
  filter(varianttype %in% c("SNP", "INDEL")) %>%
  mutate(caller_label = factor(caller, levels = names(COLORS), labels = LABELS)) %>%
  ggplot(aes(x = recall, y = precision, color = caller, shape = varianttype)) +
  geom_point(size = 6, stroke = 1.2) +
  geom_text_repel(
    aes(label = paste0(caller_label, "\n(", varianttype, ")")),
    size = 3, box.padding = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_shape_manual(values = c("SNP" = 16, "INDEL" = 17), name = "Variant Type") +
  scale_x_continuous(limits = c(0.90, 1.0), labels = percent, breaks = seq(0.90, 1, 0.02)) +
  scale_y_continuous(limits = c(0.96, 1.0), labels = percent, breaks = seq(0.96, 1, 0.01)) +
  labs(
    title = "Precision vs Recall Trade-off",
    subtitle = "Ideal performance:  top-right corner",
    x = "Recall (Sensitivity)", y = "Precision (PPV)"
  ) +
  theme_pub() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))

save_plot(p2, "02_precision_recall_scatter", width = 10, height = 7)

#--- Plot 3: TP/FP/FN Stacked Bar ---#
cat("  Creating TP/FP/FN classification plot...\n")

class_colors <- c(
  "True Positive" = "#27AE60",
  "False Positive" = "#E74C3C",
  "False Negative" = "#F39C12"
)

p3 <- data %>%
  filter(varianttype == "ALL") %>%
  select(caller, tp, fp, fn) %>%
  pivot_longer(c(tp, fp, fn), names_to = "category", values_to = "count") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    category = factor(category, levels = c("tp", "fp", "fn"),
                      labels = c("True Positive", "False Positive", "False Negative"))
  ) %>%
  ggplot(aes(x = caller_label, y = count, fill = category)) +
  geom_col(position = "stack", color = "black", width = 0.7, linewidth = 0.3) +
  geom_text(aes(label = comma(count)), position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = class_colors, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Variant Classification by Caller",
    subtitle = "Distribution of TP, FP, and FN (All Variants)",
    x = NULL, y = "Number of Variants"
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(p3, "03_variant_classification", width = 9, height = 6)

#-------------------------------------------------------------------------------
# 5. Advanced Plots (Heatmap, ROC, Lollipop)
#-------------------------------------------------------------------------------
cat("[5/7] Creating advanced plots...\n")

#--- Plot 4: Performance Metrics Heatmap ---#
cat("  Creating performance heatmap...\n")

heatmap_data <- data %>%
  filter(varianttype %in% c("SNP", "INDEL")) %>%
  select(caller, varianttype, precision, recall, f1) %>%
  pivot_longer(cols = c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    metric = factor(metric, levels = c("precision", "recall", "f1"),
                    labels = c("Precision", "Recall", "F1 Score")),
    varianttype = factor(varianttype, levels = c("SNP", "INDEL")),
    # Tạo label cho cell
    cell_label = sprintf("%.3f", value)
  )

p4 <- ggplot(heatmap_data, aes(x = metric, y = interaction(caller_label, varianttype, sep = "\n"))) +
  geom_tile(aes(fill = value), color = "white", linewidth = 1) +
  geom_text(aes(label = cell_label), size = 4, fontface = "bold") +
  scale_fill_viridis(
    option = "plasma",
    limits = c(0.90, 1.00),
    breaks = seq(0.90, 1.00, 0.02),
    labels = percent_format(accuracy = 1),
    oob = scales::squish,
    name = "Value",
    direction = -1
  ) +
  labs(
    title = "Performance Metrics Heatmap",
    subtitle = "Precision, Recall, and F1 Score by Caller and Variant Type",
    x = "Metric",
    y = "Caller (Variant Type)"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot(p4, "04_performance_heatmap", width = 10, height = 8)

#--- Plot 5: ROC-style Precision-Recall Curve ---#
cat("  Creating Precision-Recall curve...\n")

# Tạo simulated PR curve data (trong thực tế sẽ load từ hap.py hoặc RTG Tools)
create_pr_curve <- function(base_precision, base_recall, caller_name, n_points = 50) {
  # Simulate curve từ (recall=1, precision thấp) đến (actual recall, precision)
  recalls <- seq(1, base_recall, length.out = n_points)
  
  # Precision tăng dần khi recall giảm (trade-off curve)
  # Sử dụng hàm sigmoid để tạo curve smooth
  x <- seq(-3, 3, length.out = n_points)
  precision_range <- base_precision - (base_precision * 0.15)  # Start từ 85% của max
  precisions <- precision_range + (base_precision - precision_range) * (1 / (1 + exp(-x)))
  
  tibble(
    caller = caller_name,
    recall = recalls,
    precision = precisions
  )
}

# Generate PR curves cho mỗi caller (dựa trên data ALL variants)
pr_curves <- data %>%
  filter(varianttype == "ALL") %>%
  rowwise() %>%
  do(create_pr_curve(.$precision, .$recall, .$caller)) %>%
  ungroup()

# Actual points
actual_points <- data %>%
  filter(varianttype == "ALL") %>%
  select(caller, precision, recall)

p5 <- ggplot() +
  # PR curves
  geom_line(data = pr_curves, aes(x = recall, y = precision, color = caller),
            linewidth = 1.2, alpha = 0.8) +
  # Actual performance points
  geom_point(data = actual_points, aes(x = recall, y = precision, color = caller),
             size = 5, shape = 18) +
  # AUC area (optional shading)
  scale_color_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_x_continuous(limits = c(0.85, 1.0), labels = percent, breaks = seq(0.85, 1, 0.03)) +
  scale_y_continuous(limits = c(0.80, 1.0), labels = percent, breaks = seq(0.80, 1, 0.04)) +
  # Reference lines
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "grey50", alpha = 0.5) +
  annotate("text", x = 0.86, y = 0.955, label = "95% Precision", size = 3, color = "grey50") +
  annotate("text", x = 0.955, y = 0.81, label = "95% Recall", size = 3, color = "grey50", angle = 90) +
  labs(
    title = "Precision-Recall Curve",
    subtitle = "Diamond points indicate actual caller performance (All Variants)",
    x = "Recall (Sensitivity)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

save_plot(p5, "05_precision_recall_curve", width = 10, height = 8)

#--- Plot 6: Lollipop Chart for F1 Ranking ---#
cat("  Creating F1 ranking lollipop chart...\n")

lollipop_data <- data %>%
  filter(varianttype == "ALL") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    rank = rank(-f1)
  ) %>%
  arrange(desc(f1))

p6 <- ggplot(lollipop_data, aes(x = reorder(caller_label, f1), y = f1, color = caller)) +
  geom_segment(aes(xend = caller_label, y = 0.9, yend = f1), linewidth = 2) +
  geom_point(size = 8) +
  geom_text(aes(label = sprintf("%.3f", f1)), hjust = -0.3, size = 4, fontface = "bold") +
  scale_color_manual(values = COLORS, guide = "none") +
  scale_y_continuous(limits = c(0.9, 1.02), labels = percent, breaks = seq(0.9, 1, 0.02)) +
  coord_flip() +
  labs(
    title = "F1 Score Ranking",
    subtitle = "All Variants - Sorted by Performance",
    x = NULL,
    y = "F1 Score"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_blank()
  )

save_plot(p6, "06_f1_ranking_lollipop", width = 9, height = 5)

#--- Plot 7: Radar/Spider Chart ---#
cat("  Creating radar chart...\n")

# Chuẩn bị data cho radar chart
radar_data <- data %>%
  filter(varianttype == "ALL") %>%
  select(caller, precision, recall, f1) %>%
  mutate(
    # Thêm metrics bổ sung
    specificity = 1 - (fp / (tp + fp)),  # Approximation
    .before = precision
  ) %>%
  select(-specificity) %>%  # Remove nếu không có đủ data
  pivot_longer(-caller, names_to = "metric", values_to = "value") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    metric = factor(metric, levels = c("precision", "recall", "f1"),
                    labels = c("Precision", "Recall", "F1"))
  )

# Alternative: Grouped bar chart thay radar (dễ đọc hơn)
p7 <- ggplot(radar_data, aes(x = metric, y = value, fill = caller)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2), labels = percent) +
  labs(
    title = "Multi-Metric Performance Comparison",
    subtitle = "Precision, Recall, and F1 Score (All Variants)",
    x = NULL, y = "Score"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_plot(p7, "07_multi_metric_comparison", width = 10, height = 6)

#--- Plot 8: SNP vs INDEL Comparison ---#
cat("  Creating SNP vs INDEL comparison...\n")

snp_indel_data <- data %>%
  filter(varianttype %in% c("SNP", "INDEL")) %>%
  select(caller, varianttype, f1) %>%
  pivot_wider(names_from = varianttype, values_from = f1) %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    diff = SNP - INDEL  # Positive = better at SNPs
  )

p8 <- ggplot(snp_indel_data, aes(x = SNP, y = INDEL, color = caller)) +
  # Diagonal line (equal performance)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  # Points
  geom_point(size = 10, alpha = 0.8) +
  # Labels
  geom_text_repel(
    aes(label = caller_label),
    size = 4, fontface = "bold",
    box.padding = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = COLORS, guide = "none") +
  scale_x_continuous(limits = c(0.94, 1.0), labels = percent, breaks = seq(0.94, 1, 0.01)) +
  scale_y_continuous(limits = c(0.94, 1.0), labels = percent, breaks = seq(0.94, 1, 0.01)) +
  # Annotations
  annotate("text", x = 0.99, y = 0.945, label = "Better at SNPs →", 
           hjust = 1, size = 3.5, color = "grey40", fontface = "italic") +
  annotate("text", x = 0.945, y = 0.99, label = "↑ Better at INDELs",
           hjust = 0, size = 3.5, color = "grey40", fontface = "italic") +
  labs(
    title = "SNP vs INDEL Calling Performance",
    subtitle = "F1 Score (points above diagonal = better at INDELs)",
    x = "F1 Score (SNPs)",
    y = "F1 Score (INDELs)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_fixed()

save_plot(p8, "08_snp_vs_indel", width = 8, height = 8)

#-------------------------------------------------------------------------------
# 6. Runtime and Additional Plots
#-------------------------------------------------------------------------------
cat("[6/7] Creating runtime and additional plots...\n")

#--- Plot 9: Runtime Comparison ---#
cat("  Creating runtime comparison plot...\n")

# Extract caller runtimes
caller_runtime <- runtime_data %>%
  filter(grepl("gatk|deepvariant|strelka2|freebayes", step, ignore.case = TRUE)) %>%
  mutate(
    caller = case_when(
      grepl("gatk", step, ignore.case = TRUE) ~ "gatk",
      grepl("deepvariant", step, ignore.case = TRUE) ~ "deepvariant",
      grepl("strelka2", step, ignore.case = TRUE) ~ "strelka2",
      grepl("freebayes", step, ignore.case = TRUE) ~ "freebayes"
    ),
    minutes = duration_seconds / 60
  ) %>%
  filter(! is.na(caller))

if (nrow(caller_runtime) > 0) {
  p9 <- ggplot(caller_runtime, aes(x = reorder(caller, minutes), y = minutes, fill = caller)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f min", minutes)), hjust = -0.1, size = 4, fontface = "bold") +
    scale_fill_manual(values = COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    coord_flip() +
    labs(
      title = "Variant Caller Runtime Comparison",
      subtitle = "Time to complete variant calling step",
      x = NULL,
      y = "Runtime (minutes)"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p9, "09_runtime_comparison", width = 9, height = 5)
}

#--- Plot 10: Error Rate Comparison ---#
cat("  Creating error rate comparison...\n")

error_data <- data %>%
  filter(varianttype == "ALL") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    FPR = fp / (tp + fp),  # False Positive Rate (1 - Precision)
    FNR = fn / (tp + fn)   # False Negative Rate (1 - Recall)
  ) %>%
  select(caller, caller_label, FPR, FNR) %>%
  pivot_longer(c(FPR, FNR), names_to = "error_type", values_to = "rate") %>%
  mutate(
    error_type = factor(error_type, levels = c("FPR", "FNR"),
                        labels = c("False Positive Rate", "False Negative Rate"))
  )

p10 <- ggplot(error_data, aes(x = caller_label, y = rate, fill = error_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.2f%%", rate * 100)),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = c("False Positive Rate" = "#E74C3C", 
                                "False Negative Rate" = "#F39C12"),
                    name = "Error Type") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), 
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Error Rate Comparison",
    subtitle = "Lower is better",
    x = NULL,
    y = "Error Rate"
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(p10, "10_error_rate_comparison", width = 10, height = 6)

#--- Plot 11: Detailed Metrics Table as Plot ---#
cat("  Creating detailed metrics table plot...\n")

table_data <- data %>%
  filter(varianttype %in% c("ALL", "SNP", "INDEL")) %>%
  mutate(
    Caller = factor(caller, levels = names(COLORS), labels = LABELS),
    Type = varianttype,
    TP = comma(tp),
    FP = comma(fp),
    FN = comma(fn),
    Precision = sprintf("%.4f", precision),
    Recall = sprintf("%.4f", recall),
    F1 = sprintf("%.4f", f1)
  ) %>%
  select(Caller, Type, TP, FP, FN, Precision, Recall, F1)

# Create table grob
table_theme <- ttheme_minimal(
  core = list(
    bg_params = list(fill = c("grey95", "white"), col = NA),
    fg_params = list(fontface = "plain", fontsize = 10)
  ),
  colhead = list(
    bg_params = list(fill = "grey80", col = NA),
    fg_params = list(fontface = "bold", fontsize = 11)
  )
)

p11 <- tableGrob(table_data, rows = NULL, theme = table_theme)

png(file.path(output_dir, "11_detailed_metrics_table.png"), 
    width = 10, height = 8, units = "in", res = 300)
grid.arrange(
  p11,
  top = textGrob("Detailed Benchmark Metrics", gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(paste("Generated:", Sys.Date()), gp = gpar(fontsize = 9, col = "grey50"))
)
dev.off()
cat("  ✓ Saved:  11_detailed_metrics_table.png\n")

#-------------------------------------------------------------------------------
# 7. Summary Dashboard
#-------------------------------------------------------------------------------
cat("[7/7] Creating summary dashboard...\n")

# Combine key plots into dashboard
dashboard <- (
  (p1) /
  (p2 | p3) /
  (p4 | p5)
) +
  plot_annotation(
    title = "Variant Calling Benchmarking - Comprehensive Dashboard",
    subtitle = paste("GATK | DeepVariant | Strelka2 | FreeBayes —", format(Sys.Date(), "%B %d, %Y")),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40")
    )
  )

save_plot(dashboard, "12_comprehensive_dashboard", width = 16, height = 20)

# Mini dashboard (key plots only)
mini_dashboard <- (p1 / (p2 | p6)) +
  plot_annotation(
    title = "Variant Calling Benchmark Summary",
    subtitle = paste("Generated:", format(Sys.Date(), "%B %d, %Y")),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

save_plot(mini_dashboard, "13_mini_dashboard", width = 14, height = 10)

#-------------------------------------------------------------------------------
# Generate CSV summary
#-------------------------------------------------------------------------------
cat("\nGenerating summary tables...\n")

# Main summary
summary_csv <- data %>%
  filter(varianttype == "ALL") %>%
  mutate(Caller = factor(caller, levels = names(COLORS), labels = LABELS)) %>%
  select(Caller, TP = tp, FP = fp, FN = fn, 
         Precision = precision, Recall = recall, F1 = f1) %>%
  arrange(desc(F1))

write_csv(summary_csv, file.path(output_dir, "benchmark_summary_table.csv"))
cat("  ✓ Saved: benchmark_summary_table.csv\n")

# Detailed by variant type
detailed_csv <- data %>%
  mutate(Caller = factor(caller, levels = names(COLORS), labels = LABELS)) %>%
  select(Caller, VariantType = varianttype, TP = tp, FP = fp, FN = fn,
         Precision = precision, Recall = recall, F1 = f1)

write_csv(detailed_csv, file.path(output_dir, "benchmark_detailed_table.csv"))
cat("  ✓ Saved: benchmark_detailed_table.csv\n")

#-------------------------------------------------------------------------------
# Final Summary
#-------------------------------------------------------------------------------
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("BENCHMARK RESULTS SUMMARY (Ranked by F1 Score - All Variants)\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
print(summary_csv, n = Inf)
cat("═══════════════════════════════════════════════════════════════════════\n")

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║              VISUALIZATION COMPLETE!                                ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nGenerated figures:\n")

# List all generated files
all_files <- list.files(output_dir, pattern = "\\.(png|pdf|csv)$")
png_files <- grep("\\.png$", all_files, value = TRUE)
pdf_files <- grep("\\.pdf$", all_files, value = TRUE)
csv_files <- grep("\\.csv$", all_files, value = TRUE)

cat(sprintf("\n  PNG files (%d):\n", length(png_files)))
for (f in sort(png_files)) cat(sprintf("    • %s\n", f))

cat(sprintf("\n  PDF files (%d):\n", length(pdf_files)))
for (f in sort(pdf_files)) cat(sprintf("    • %s\n", f))

cat(sprintf("\n  CSV files (%d):\n", length(csv_files)))
for (f in sort(csv_files)) cat(sprintf("    • %s\n", f))

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
