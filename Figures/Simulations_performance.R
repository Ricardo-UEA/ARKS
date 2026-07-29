#!/usr/bin/env Rscript
# ==============================================================================
# make_fig3.R - ARKS pre- vs post-depletion CAMISIM benchmark, Figure 3
# ==============================================================================
# Rebuilds Fig. 3a-d from the corrected evaluate_arks_depletion.py outputs
# (arks_depletion_eval_*.csv), replacing the earlier version of this figure
# that used numbers from the buggy threshold-mismatched pipeline.
#
#   Panel a: precision / recall / F1 (+-1 SD) across filter configurations,
#            pre-depletion, ordered by F1.
#   Panel b: precision-recall space for the same configurations, with iso-F1
#            contours; the selected cov >= 0.005 threshold highlighted in red.
#   Panel c: per-sample precision / recall / F1, pre-depletion vs post-ARKS,
#            at the cov >= 0.005 threshold.
#   Panel d: per-sample false-positive read fraction, pre-depletion vs
#            post-ARKS (log scale; dashed line at the 1% reference used in
#            the text).
#
# NOTE ON PANEL A/B: the original figure had 8 filter configurations,
# including "reads >= 100, cov >= 0.001". The current threshold-sweep CSVs
# only have 7 - evaluate_arks_depletion.py's sweep varies min_kmers/min_cov
# but never tried a min_reads-based config. Panels a/b below show the 7 that
# exist. Say the word if you want that 8th config added to the Python script
# and rerun - it's a one-line addition (apply_filter() already supports
# min_reads).
#
# Requires: ggplot2, dplyr, tidyr, readr, patchwork, scales, ggrepel
#   install.packages(c("ggplot2","dplyr","tidyr","readr","patchwork","scales","ggrepel"))
#
# Usage:
#   Rscript make_fig3.R
# (edit base_dir below if your files live somewhere else)
# ==============================================================================

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales", "ggrepel")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with: install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

base_dir <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Simultations/arks_depletion_eval"

message(
  "Note: panel a/b show 7 of the original 8 filter configs - ",
  "'reads >= 100, cov >= 0.001' isn't in the current threshold-sweep data. ",
  "See the header comment for how to add it back."
)

# ------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------
sweep_pre <- read_csv(
  file.path(base_dir, "arks_depletion_eval_threshold_comparison_pre.csv"),
  show_col_types = FALSE
)
pre_metrics <- read_csv(
  file.path(base_dir, "arks_depletion_eval_pre_per_sample_metrics.csv"),
  show_col_types = FALSE
)
post_metrics <- read_csv(
  file.path(base_dir, "arks_depletion_eval_post_per_sample_metrics.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------------------
# Shared config labels / categories for panels a & b
# ------------------------------------------------------------------------
config_labels <- c(
  "No filter"                          = "No filter",
  "cov >= 0.001 only"                  = "cov ≥ 0.001 only",
  "kmers >= 100 only"                  = "k-mers ≥ 100 only",
  "kmers >= 100, cov >= 0.001"          = "k-mers ≥ 100, cov ≥ 0.001",
  "kmers >= 500, cov >= 0.001"          = "k-mers ≥ 500, cov ≥ 0.001",
  "kmers >= 100, cov >= 0.002"          = "k-mers ≥ 100, cov ≥ 0.002",
  "cov >= 0.005 (manuscript threshold)" = "cov ≥ 0.005 (optimal F1)"
)

sweep_pre <- sweep_pre %>%
  mutate(
    Label = unname(config_labels[Config]),
    FilterType = case_when(
      Config == "No filter" ~ "No filter",
      min_kmers > 0 & min_cov > 0 ~ "Combined",
      TRUE ~ "Single criterion"
    ),
    IsOptimal = Config == "cov >= 0.005 (manuscript threshold)"
  )

label_order <- sweep_pre$Label[order(sweep_pre$F1)]

# ------------------------------------------------------------------------
# Panel a: precision / recall / F1 (+-1 SD) across filter configs
# ------------------------------------------------------------------------
panel_a_data <- sweep_pre %>%
  select(Label, Precision, Precision_std, Recall, Recall_std, F1, F1_std) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1),
    names_to = "Metric", values_to = "Score"
  ) %>%
  mutate(
    SD = case_when(
      Metric == "Precision" ~ Precision_std,
      Metric == "Recall"    ~ Recall_std,
      Metric == "F1"        ~ F1_std
    ),
    Metric = factor(Metric, levels = c("Precision", "Recall", "F1")),
    Label = factor(Label, levels = label_order)
  )

metric_colors <- c(Precision = "#E41A1C", Recall = "#377EB8", F1 = "#4DAF4A")
metric_shapes <- c(Precision = 16, Recall = 17, F1 = 15)

panel_a <- ggplot(panel_a_data, aes(x = Score, y = Label, color = Metric, shape = Metric)) +
  geom_segment(
    aes(x = Score - SD, xend = Score + SD, y = Label, yend = Label),
    position = position_dodge(width = 0.6)
  ) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = metric_shapes) +
  labs(x = "Score (± 1 SD)", y = NULL, color = "Metric", shape = "Metric") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# ------------------------------------------------------------------------
# Panel b: precision-recall space with iso-F1 contours
# ------------------------------------------------------------------------
iso_f1 <- lapply(c(0.1, 0.2, 0.4, 0.6, 0.8), function(f1) {
  r_start <- f1 / (2 - f1)  # recall where precision = 1
  r <- seq(r_start, 1, length.out = 200)
  p <- f1 * r / (2 * r - f1)
  data.frame(Recall = r, Precision = p, F1level = factor(f1))
}) %>% bind_rows()

point_colors_b <- c(
  "No filter" = "grey50", "Single criterion" = "#377EB8", "Combined" = "#4DAF4A"
)

panel_b <- ggplot() +
  geom_line(
    data = iso_f1, aes(x = Recall, y = Precision, group = F1level),
    linetype = "dashed", color = "grey75"
  ) +
  geom_point(
    data = sweep_pre %>% filter(!IsOptimal),
    aes(x = Recall, y = Precision, color = FilterType), size = 3
  ) +
  geom_point(
    data = sweep_pre %>% filter(IsOptimal),
    aes(x = Recall, y = Precision), color = "#E41A1C", size = 3.6
  ) +
  geom_text_repel(
    data = sweep_pre, aes(x = Recall, y = Precision, label = Label),
    size = 2.6, max.overlaps = 20, seed = 1
  ) +
  scale_color_manual(values = point_colors_b, name = "Filter type") +
  coord_cartesian(xlim = c(0.8, 1), ylim = c(0, 1)) +
  labs(x = "Recall (sensitivity)", y = "Precision (PPV)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# ------------------------------------------------------------------------
# Panels c & d: per-sample pre- vs post-ARKS comparison
# ------------------------------------------------------------------------
combined_metrics <- bind_rows(
  pre_metrics %>% mutate(Condition = "Pre-depletion"),
  post_metrics %>% mutate(Condition = "Post-ARKS")
) %>%
  mutate(Condition = factor(Condition, levels = c("Pre-depletion", "Post-ARKS")))

condition_colors <- c("Pre-depletion" = "#7FA8C9", "Post-ARKS" = "#2E8B57")

# Panel c: precision / recall / F1
panel_c_data <- combined_metrics %>%
  select(Sample, Condition, Precision_PPV, Recall_Sensitivity, F1) %>%
  pivot_longer(
    cols = c(Precision_PPV, Recall_Sensitivity, F1),
    names_to = "Metric", values_to = "Score"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "Precision_PPV" ~ "Precision",
      Metric == "Recall_Sensitivity" ~ "Recall",
      TRUE ~ "F1"
    ),
    Metric = factor(Metric, levels = c("Precision", "Recall", "F1"))
  )

panel_c <- ggplot(panel_c_data, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(alpha = 0.6, trim = TRUE, color = NA) +
  geom_jitter(width = 0.08, size = 0.6, alpha = 0.4, color = "grey30") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 1.5, color = "black") +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# Panel d: false-positive read fraction
panel_d <- ggplot(combined_metrics, aes(x = Condition, y = FP_Read_Fraction_Pct, fill = Condition)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_violin(alpha = 0.6, trim = TRUE, color = NA) +
  geom_jitter(width = 0.08, size = 0.6, alpha = 0.4, color = "grey30") +
  scale_fill_manual(values = condition_colors) +
  scale_y_log10(labels = label_number(suffix = "%")) +
  labs(
    x = NULL, y = "False-positive read fraction (%, log scale)",
    caption = "Reads assigned to genera absent from CAMISIM ground truth (n = 100 simulations)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# ------------------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------------------
fig3 <- (panel_a + panel_b) / (panel_c + panel_d) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(base_dir, "Fig3_updated.pdf"), fig3, width = 13, height = 10, dpi = 400)
ggsave(file.path(base_dir, "Fig3_updated.png"), fig3, width = 13, height = 10, dpi = 400)

cat("Saved Fig3_updated.pdf and Fig3_updated.png to", base_dir, "\n")