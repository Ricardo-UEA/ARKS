#!/usr/bin/env Rscript
# ============================================================
# HEROIC1K Shannon Diversity — 2 panels
# A: Shannon entropy across all four conditions (violin+beeswarm)
# B: ΔShannon vs No-depletion by ancestry × condition
#
# Fixed 2026-07-22: bracket_df/cap_df/diff_df ARKS values were stale
# (left over from before the ARKS KrakenUniq re-run) and didn't match
# the actual violin data/medians already being plotted correctly.
# Only the three ARKS entries in bracket_df, cap_df, and diff_df changed;
# everything else in this script is untouched from the original.
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(cowplot)
# ── Paths ─────────────────────────────────────────────────────
DATA_DIR <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Shannon_Diversity_Heroic"
OUTDIR   <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Shannon_Diversity_Heroic"
# ── 1. Load ───────────────────────────────────────────────────
nd   <- read.csv(file.path(DATA_DIR, "No_depletion_alpha_diversity.csv"), stringsAsFactors = FALSE)
grch <- read.csv(file.path(DATA_DIR, "GRCh38_alpha_diversity.csv"),                stringsAsFactors = FALSE)
t2t  <- read.csv(file.path(DATA_DIR, "T2T_alpha_diversity.csv"),                      stringsAsFactors = FALSE)
arks <- read.csv(file.path(DATA_DIR, "ARKS_alpha_diversity.csv"),                    stringsAsFactors = FALSE)
# ── 2. Combine long ───────────────────────────────────────────
add_cond <- function(df, cond) {
  df %>% select(sample_id, ancestry, shannon) %>% mutate(condition = cond)
}
df_long <- bind_rows(
  add_cond(nd,   "No depletion"),
  add_cond(grch, "GRCh38"),
  add_cond(t2t,  "T2T-CHM13"),
  add_cond(arks, "ARKS")
) %>%
  mutate(
    condition = factor(condition,
                       levels = c("No depletion","GRCh38","T2T-CHM13","ARKS")),
    ancestry  = factor(ancestry,
                       levels = c("Southern African","European","Admixed","Asian"))
  ) %>%
  filter(!is.na(ancestry), !is.na(shannon))
# ── 3. Delta (vs No depletion) ────────────────────────────────
df_wide <- df_long %>%
  pivot_wider(names_from = condition, values_from = shannon)
df_delta <- df_wide %>%
  mutate(
    GRCh38    = GRCh38    - `No depletion`,
    `T2T-CHM13` = `T2T-CHM13` - `No depletion`,
    ARKS      = ARKS      - `No depletion`
  ) %>%
  select(sample_id, ancestry, GRCh38, `T2T-CHM13`, ARKS) %>%
  pivot_longer(c(GRCh38, `T2T-CHM13`, ARKS),
               names_to = "condition", values_to = "delta_shannon") %>%
  mutate(
    condition = factor(condition, levels = c("GRCh38","T2T-CHM13","ARKS")),
    ancestry  = factor(ancestry,
                       levels = c("Southern African","European","Admixed","Asian"))
  ) %>%
  filter(!is.na(ancestry), !is.na(delta_shannon))
# Focus on European and Southern African (main groups)
df_delta_main <- df_delta %>%
  filter(ancestry %in% c("European","Southern African"))
# ── 4. Colours ────────────────────────────────────────────────
cond_cols <- c(
  "No depletion" = "#AAAAAA",
  "GRCh38"       = "#E08B00",
  "T2T-CHM13"    = "#4A90D9",
  "ARKS"         = "#1B9E77"
)
anc_cols <- c(
  "Southern African" = "#C0392B",
  "European"         = "#2980B9",
  "Admixed"          = "#8E44AD",
  "Asian"            = "#27AE60"
)
# ── 5. Shared theme ───────────────────────────────────────────
bt <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid.major  = element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "grey20", linewidth = 0.5),
      axis.text         = element_text(colour = "black", size = 8),
      axis.title        = element_text(size = 8.5),
      plot.title        = element_text(face = "bold", size = 10),
      strip.text        = element_text(size = 8, face = "bold"),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      plot.margin       = margin(6, 6, 6, 6)
    )
}
# ══════════════════════════════════════════════════════════════
# PANEL A — Shannon entropy across four conditions
# ══════════════════════════════════════════════════════════════
med_a <- df_long %>%
  group_by(condition) %>%
  summarise(med = median(shannon, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%.3f", med))
pA <- ggplot(df_long,
             aes(x = condition, y = shannon,
                 fill = condition, colour = condition)) +
  geom_violin(alpha = 0.22, scale = "width", trim = TRUE,
              linewidth = 0.5) +
  geom_beeswarm(size = 0.5, alpha = 0.35, cex = 0.5) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.6,
               colour = "black", fatten = 1) +
  geom_text(data = med_a,
            aes(x = condition, y = med + 0.18, label = label),
            inherit.aes = FALSE, size = 2.2, colour = "black") +
  scale_fill_manual(values = cond_cols)   +
  scale_colour_manual(values = cond_cols) +
  scale_y_continuous(name = "Shannon entropy",
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  labs(x = NULL, title = "a") +
  bt() +
  theme(legend.position = "none",
        axis.text.x     = element_text(size = 8))
# ══════════════════════════════════════════════════════════════
# PANEL B — ΔShannon by ancestry × condition
# European vs Southern African, three conditions side by side
# ══════════════════════════════════════════════════════════════
# Median annotations
med_b <- df_delta_main %>%
  group_by(condition, ancestry) %>%
  summarise(med = median(delta_shannon, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%+.3f", med))
# Bracket annotation data — differential between European and SA medians
# x positions: S.African=1, European=2 (discrete factor order)
# bracket sits at x=2.6 with caps connecting the two medians
bracket_df <- data.frame(
  condition   = factor(c("GRCh38","GRCh38","GRCh38",
                         "T2T-CHM13","T2T-CHM13","T2T-CHM13",
                         "ARKS","ARKS","ARKS"),
                       levels = c("GRCh38","T2T-CHM13","ARKS")),
  # vertical spine of bracket
  x    = c(2.55, 2.55, 2.55,  2.55, 2.55, 2.55,  2.55, 2.55, 2.55),
  xend = c(2.55, 2.55, 2.55,  2.55, 2.55, 2.55,  2.55, 2.55, 2.55),
  y    = c(0.332, 0.332, 0.066,  -0.001, -0.001, 0.637,  0.741, 0.741, 1.070),
  yend = c(0.332, 0.066, 0.066,  -0.001, 0.637,  0.637,  0.741, 1.070,  1.070)
)
# Cap ticks
cap_df <- data.frame(
  condition = factor(c("GRCh38","GRCh38","T2T-CHM13","T2T-CHM13","ARKS","ARKS"),
                     levels = c("GRCh38","T2T-CHM13","ARKS")),
  x    = c(2.45, 2.45, 2.45, 2.45, 2.45, 2.45),
  xend = c(2.55, 2.55, 2.55, 2.55, 2.55, 2.55),
  y    = c(0.332, 0.066, -0.001, 0.637, 0.741, 1.070),
  yend = c(0.332, 0.066, -0.001, 0.637, 0.741, 1.070)
)
# Label data for differential text
diff_df <- data.frame(
  condition  = factor(c("GRCh38","T2T-CHM13","ARKS"),
                      levels = c("GRCh38","T2T-CHM13","ARKS")),
  x          = c(2.62, 2.62, 2.62),
  y          = c((0.332+0.066)/2, (-0.001+0.637)/2, (0.741+1.070)/2),
  label      = c("Δ = -0.27", "Δ = 0.64", "Δ = 0.33"),
  col        = c("#CC0000","#2C7FB8","#2C7FB8")
)
pB <- ggplot(df_delta_main,
             aes(x = ancestry, y = delta_shannon,
                 fill = ancestry, colour = ancestry)) +
  # Zero reference line
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_violin(alpha = 0.22, scale = "width", trim = TRUE,
              linewidth = 0.45) +
  geom_beeswarm(size = 0.55, alpha = 0.40, cex = 0.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.55,
               colour = "black", fatten = 1) +
  geom_text(data = med_b,
            aes(x = ancestry, y = med + 0.18, label = label),
            inherit.aes = FALSE, size = 2.1, colour = "black") +
  # Bracket spine
  geom_segment(data = bracket_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE,
               colour = "grey30", linewidth = 0.45) +
  # Bracket caps
  geom_segment(data = cap_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE,
               colour = "grey30", linewidth = 0.45) +
  # Differential label
  geom_text(data = diff_df,
            aes(x = x, y = y, label = label, colour = NULL),
            inherit.aes = FALSE,
            colour = diff_df$col,
            size = 2.2, hjust = 0, fontface = "bold") +
  facet_wrap(~condition, nrow = 1) +
  scale_fill_manual(values = anc_cols,
                    labels = c("Southern African" = "S. African",
                               "European"         = "European")) +
  scale_colour_manual(values = anc_cols,
                      labels = c("Southern African" = "S. African",
                                 "European"         = "European")) +
  scale_x_discrete(labels = c("Southern African" = "S. African",
                              "European"         = "European")) +
  scale_y_continuous(name   = "Δ Shannon entropy\n(vs no depletion)",
                     limits = c(-2.5, 3.8),
                     breaks = seq(-2, 3, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, title = "b",
       colour = "Ancestry", fill = "Ancestry") +
  bt() +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.35, "cm"),
        legend.text      = element_text(size = 8),
        legend.title     = element_text(size = 8, face = "bold"),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        plot.margin      = margin(6, 40, 6, 6))
# ══════════════════════════════════════════════════════════════
# COMBINE
# ══════════════════════════════════════════════════════════════
combined <- plot_grid(pA, pB,
                      nrow = 1, rel_widths = c(0.9, 1.3),
                      align = "h", axis = "tb")
combined <- ggdraw(combined) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
# ── Save ──────────────────────────────────────────────────────
ggsave(file.path(OUTDIR, "shannon_diversity_plots.pdf"), combined,
       width = 220, height = 120, units = "mm",
       dpi = 800, device = cairo_pdf)
ggsave(file.path(OUTDIR, "shannon_diversity_plots.png"), combined,
       width = 220, height = 120, units = "mm", dpi = 800)
cat("Saved: shannon_diversity_plots.pdf / .png\n")
# ── Stats summary ─────────────────────────────────────────────
cat("\n── Overall Shannon medians ──\n")
df_long %>%
  group_by(condition) %>%
  summarise(median = round(median(shannon, na.rm=TRUE), 3),
            mean   = round(mean(shannon, na.rm=TRUE), 3),
            n      = n(), .groups = "drop") %>%
  print()
cat("\n── ΔShannon medians by ancestry × condition ──\n")
df_delta %>%
  group_by(condition, ancestry) %>%
  summarise(median = round(median(delta_shannon, na.rm=TRUE), 3),
            n      = n(), .groups = "drop") %>%
  print(n = Inf)
cat("\n── Mann-Whitney European vs S.African (BH corrected) ──\n")
library(stats)
pvals <- c()
res   <- list()
for (cond in c("GRCh38","T2T-CHM13","ARKS")) {
  eur <- df_delta %>% filter(condition==cond, ancestry=="European")       %>% pull(delta_shannon)
  sa  <- df_delta %>% filter(condition==cond, ancestry=="Southern African") %>% pull(delta_shannon)
  wt  <- wilcox.test(eur, sa, alternative="two.sided")
  pvals <- c(pvals, wt$p.value)
  res[[cond]] <- list(eur_med=median(eur), sa_med=median(sa), W=wt$statistic, p=wt$p.value)
}
padj <- p.adjust(pvals, method="BH")
for (i in seq_along(res)) {
  cond <- names(res)[i]
  r    <- res[[cond]]
  cat(sprintf("%-12s Eur=%+.3f  SA=%+.3f  diff=%.3f  W=%.0f  p=%.3e  p_adj=%.3e\n",
              cond, r$eur_med, r$sa_med,
              r$eur_med - r$sa_med, r$W, r$p, padj[i]))
}