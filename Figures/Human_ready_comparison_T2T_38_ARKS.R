# ── HEROIC1K Human Read Depletion Plots — v3 ──────────────────────────────────
# Changes from v2:
#   - GRCh38 now orange (#E08B00)
#   - Tick labels staggered to prevent overlap on panel a x-axis
#   - Median labels confirmed (not means)
#   - Median label positions adjusted to avoid crossbar overlap
#   - Panel b median labels right-justified and truncated to fit

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(cowplot)

df <- read.csv(
  "/Users/hce24xau/Desktop/Genomics_Benchmarking/human_reads_summary/human_reads_summary.csv",
  stringsAsFactors = FALSE
)

df$condition <- factor(df$condition,
                       levels = c("No_depletion", "GRCh38", "T2T", "ARKS"),
                       labels = c("Pre-depletion", "GRCh38", "T2T-CHM13", "ARKS"))

# ── Colours — GRCh38 now orange ───────────────────────────────────────────────
cols <- c(
  "Pre-depletion" = "#AAAAAA",
  "GRCh38"        = "#E08B00",
  "T2T-CHM13"     = "#4A90D9",
  "ARKS"          = "#1B9E77"
)

# ── Medians ───────────────────────────────────────────────────────────────────
med_df <- df %>%
  group_by(condition) %>%
  summarise(med = median(human_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%.3f%%", med))

# ── PANEL A ───────────────────────────────────────────────────────────────────
p1 <- ggplot(df, aes(x = condition, y = human_pct,
                     fill = condition, colour = condition)) +
  geom_violin(alpha = 0.25, scale = "width", trim = TRUE, linewidth = 0.5) +
  geom_beeswarm(size = 0.65, alpha = 0.5, cex = 0.75) +
  stat_summary(fun = median, geom = "crossbar", width = 0.38,
               linewidth = 0.6, colour = "black", fatten = 1.2) +
  # Median labels above crossbar — nudge up enough to clear it
  geom_text(data = med_df,
            aes(x = condition, y = med * 1.8, label = label),
            inherit.aes = FALSE, size = 2.3, colour = "black") +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_y_log10(
    limits = c(0.0004, 25),
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = c(0.02, 0)
  ) +
  annotation_logticks(sides = "l", size = 0.25, colour = "grey60") +
  labs(x = NULL, y = "Human read fraction (%, log scale)", title = "a") +
  theme_cowplot(font_size = 9) +
  background_grid(major = "y", minor = "none",
                  colour.major = "grey93", size.major = 0.3) +
  theme(
    legend.position  = "none",
    # Stagger x-axis labels to prevent overlap
    axis.text.x      = element_text(size = 7.5,
                                    margin = margin(t = 2)),
    axis.title.y     = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  # Stagger by offsetting alternate labels using scale override
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# ── PANEL B ───────────────────────────────────────────────────────────────────
df_wide <- df %>%
  select(sample_id, condition, human_pct) %>%
  pivot_wider(names_from = condition, values_from = human_pct) %>%
  arrange(desc(`Pre-depletion`)) %>%
  mutate(rank = row_number())

df_long <- df_wide %>%
  pivot_longer(cols = c("Pre-depletion","GRCh38","T2T-CHM13","ARKS"),
               names_to = "condition", values_to = "human_pct") %>%
  mutate(condition = factor(condition, levels = names(cols)))

med_lines <- df %>%
  group_by(condition) %>%
  summarise(med = median(human_pct, na.rm = TRUE), .groups = "drop")

# Stagger median label y-positions slightly to avoid overlap (T2T and ARKS close)
med_lines <- med_lines %>%
  mutate(label_y = case_when(
    condition == "T2T-CHM13"    ~ med * 1.25,
    condition == "ARKS"         ~ med * 0.80,
    condition == "GRCh38"       ~ med * 1.08,
    condition == "Pre-depletion"~ med * 1.08,
    TRUE ~ med
  ))

p2 <- ggplot(df_long, aes(x = rank, y = human_pct,
                          colour = condition, shape = condition)) +
  geom_point(size = 0.65, alpha = 0.55) +
  geom_hline(data = med_lines,
             aes(yintercept = med, colour = condition),
             linetype = "dashed", linewidth = 0.45, alpha = 0.85) +
  # Right-margin labels at staggered positions
  geom_text(data = med_lines,
            aes(x = 183, y = label_y,
                label = sprintf("%.3f%%", med),
                colour = condition),
            inherit.aes = FALSE, hjust = 0, size = 2.0) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(
    values = c("Pre-depletion"=16,"GRCh38"=17,"T2T-CHM13"=15,"ARKS"=18),
    name = NULL) +
  scale_y_log10(
    limits = c(0.0004, 25),
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = c(0.02, 0)
  ) +
  scale_x_continuous(limits = c(1, 200), expand = c(0.01, 0)) +
  annotation_logticks(sides = "l", size = 0.25, colour = "grey60") +
  labs(x = "Samples (ranked by pre-depletion human read fraction)",
       y = "Human read fraction (%, log scale)", title = "b") +
  theme_cowplot(font_size = 9) +
  background_grid(major = "y", minor = "none",
                  colour.major = "grey93", size.major = 0.3) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.35, "cm"),
    legend.text      = element_text(size = 7.5),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.x     = element_text(size = 8),
    axis.title.y     = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(alpha=1, size=2,
                                                   linewidth=1.2)),
         shape = guide_legend(nrow = 1))

# ── PANEL C ───────────────────────────────────────────────────────────────────
df_delta <- df_wide %>%
  mutate(
    GRCh38      = `Pre-depletion` - GRCh38,
    `T2T-CHM13` = `Pre-depletion` - `T2T-CHM13`,
    ARKS        = `Pre-depletion` - ARKS
  ) %>%
  pivot_longer(cols = c("GRCh38","T2T-CHM13","ARKS"),
               names_to = "tool", values_to = "removed_pct") %>%
  mutate(tool = factor(tool, levels = c("GRCh38","T2T-CHM13","ARKS")))

delta_cols <- c("GRCh38"="#E08B00","T2T-CHM13"="#4A90D9","ARKS"="#1B9E77")

delta_meds <- df_delta %>%
  group_by(tool) %>%
  summarise(med = median(removed_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%.3f%%", med))

p3 <- ggplot(df_delta, aes(x = tool, y = removed_pct,
                           fill = tool, colour = tool)) +
  geom_violin(alpha = 0.25, scale = "width", linewidth = 0.5) +
  geom_beeswarm(size = 0.65, alpha = 0.5, cex = 0.75) +
  stat_summary(fun = median, geom = "crossbar", width = 0.38,
               linewidth = 0.6, colour = "black", fatten = 1.2) +
  geom_text(data = delta_meds,
            aes(x = tool, y = med * 1.8, label = label),
            inherit.aes = FALSE, size = 2.3, colour = "black") +
  scale_fill_manual(values = delta_cols) +
  scale_colour_manual(values = delta_cols) +
  scale_y_log10(
    limits  = c(0.0002, 5),
    breaks  = c(0.001, 0.01, 0.1, 1),
    labels  = c("0.001%", "0.01%", "0.1%", "1%"),
    expand  = c(0.02, 0)
  ) +
  annotation_logticks(sides = "l", size = 0.25, colour = "grey60") +
  labs(x = NULL, y = "Human reads removed (%, log scale)", title = "c") +
  theme_cowplot(font_size = 9) +
  background_grid(major = "y", minor = "none",
                  colour.major = "grey93", size.major = 0.3) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 7.5),
    axis.title.y     = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# ── Combine ────────────────────────────────────────────────────────────────────
combined <- plot_grid(p1, p2, p3,
                      nrow = 1,
                      rel_widths = c(1.0, 1.6, 1.0),
                      align = "h")

combined <- ggdraw(combined) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# ── Save ───────────────────────────────────────────────────────────────────────
outdir <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/human_reads_summary"
dir.create(outdir, showWarnings = FALSE)

ggsave(file.path(outdir, "heroic1k_fig3_v3.pdf"),
       combined, width = 190, height = 90, units = "mm", dpi = 800)
ggsave(file.path(outdir, "heroic1k_fig3_v3.png"),
       combined, width = 190, height = 90, units = "mm", dpi = 800)
cat("Saved: heroic1k_fig3_v3.pdf / .png\n")

# ── Stats ──────────────────────────────────────────────────────────────────────
cat("\n=== SUMMARY STATISTICS (MEDIANS) ===\n")
df %>%
  group_by(condition) %>%
  summarise(n = n(),
            median = round(median(human_pct), 4),
            mean   = round(mean(human_pct), 4),
            IQR_lo = round(quantile(human_pct, 0.25), 4),
            IQR_hi = round(quantile(human_pct, 0.75), 4),
            .groups = "drop") %>% print()

arks <- df %>% filter(condition == "ARKS")          %>% arrange(sample_id) %>% pull(human_pct)
t2t  <- df %>% filter(condition == "T2T-CHM13")     %>% arrange(sample_id) %>% pull(human_pct)
grch <- df %>% filter(condition == "GRCh38")        %>% arrange(sample_id) %>% pull(human_pct)
pre  <- df %>% filter(condition == "Pre-depletion") %>% arrange(sample_id) %>% pull(human_pct)

cat("\n=== PAIRED WILCOXON TESTS ===\n")
for (pair in list(list(arks, t2t, "ARKS vs T2T-CHM13"),
                  list(arks, grch, "ARKS vs GRCh38"),
                  list(t2t,  grch, "T2T-CHM13 vs GRCh38"),
                  list(arks, pre,  "ARKS vs Pre-depletion"))) {
  wt <- wilcox.test(pair[[1]], pair[[2]], paired = TRUE)
  cat(sprintf("%-30s W = %6.0f,  p = %.2e\n", pair[[3]], wt$statistic, wt$p.value))
}