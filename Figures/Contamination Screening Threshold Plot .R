#!/usr/bin/env Rscript
# ============================================================
# Contamination Screening Threshold Plot — Grid Layout
# One panel per pan-genome reference, plus two summary panels
# for total contamination k-mers removed (Unique + Assembled stages)
# Data from: Unique_kmers_classification.csv
# ============================================================
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────
DATA_DIR    <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Unique_kmers_classifications/Unique Kraken csv"
OUTPUT_PATH <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Unique_kmers_classifications/Contamination_screen_grid.pdf"

KMER_THRESH    <- 100
COV_THRESH     <- 0.01   # applied to all references except 2nd Screening
COV_THRESH_2ND <- 0.005  # refined threshold for 2nd Screening Post Assembly

EXCLUDE_RANKS <- c(
  "no rank", "superkingdom", "kingdom", "subkingdom",
  "clade", "phylum", "subphylum", "class", "subclass",
  "order", "suborder", "infraorder", "family", "superfamily"
)
EXCLUDE_NAMES <- c(
  "unclassified", "root", "cellular organisms",
  "Bacteria", "Eukaryota", "Archaea", "Viruses",
  "unclassified viruses", "unclassified DNA viruses",
  "unclassified dsDNA viruses", "other entries",
  "other sequences", "artificial sequences",
  # Host (Homo sapiens) lineage nodes at ranks NOT already covered by
  # EXCLUDE_RANKS (superclass/superorder/parvorder/subfamily/genus/species).
  # In khoe_san_pan_genome.csv these six rows all carry the full host
  # k-mer/coverage total (kmers = 2,401,285,368; cov = 0.9262) — i.e. they
  # are the assembly's own sequence, not contamination — and without this
  # exclusion they swamp the plot as the top "excluded" hits by several
  # orders of magnitude. Ancestor ranks already excluded via EXCLUDE_RANKS:
  # Chordata (phylum), Mammalia (class), Primates (order), Hominidae (family).
  "Sarcopterygii", "Euarchontoglires", "Catarrhini", "Homininae",
  "Homo", "Homo sapiens",
  # Database/library artefacts, not biological contaminants.
  "synthetic construct"
)

# ── Colour palette — one per reference ───────────────────────
ref_colours <- c(
  "Arab Pan Genome"                        = "#E69F00",
  "African Pan Genome"                     = "#1A9850",
  "Khoe-San Pan Genome"                    = "#D55E00",
  "Chinese Pan Genome"                     = "#F46D43",
  "Human Pan Genome Consortium"            = "#7B2D8B",
  "2nd Screening Post Assembly of K-mers"  = "#0072B2",
  "1st Screening Pre-Assembly of K-mers"   = "#56B4E9"
)

# Display labels
ref_labels <- c(
  "Arab Pan Genome"                        = "Arab Pan Genome",
  "African Pan Genome"                     = "African Pan Genome",
  "Khoe-San Pan Genome"                    = "Khoe-San Pan Genome",
  "Chinese Pan Genome"                     = "Chinese Pan Genome",
  "Human Pan Genome Consortium"            = "Human Pan Genome Consortium",
  "2nd Screening Post Assembly of K-mers"  = "2nd Screening\n(Post-Assembly)"
)

# ── Helper: abbreviate species names ─────────────────────────
short_tax_label <- function(x) {
  sapply(x, function(name) {
    name  <- trimws(name)
    words <- unlist(strsplit(name, "\\s+"))
    if (length(words) >= 2) {
      paste0(substr(words[1], 1, 1), ". ", paste(words[-1], collapse = " "))
    } else {
      name
    }
  })
}

# ── File map: reference label → CSV path ─────────────────────
FILES <- list(
  "Arab Pan Genome"                        = file.path(DATA_DIR, "Arab_pan_genome.csv"),
  "African Pan Genome"                     = file.path(DATA_DIR, "african_pan_genome.csv"),
  "Khoe-San Pan Genome"                    = file.path(DATA_DIR, "khoe_san_pan_genome.csv"),
  "Chinese Pan Genome"                     = file.path(DATA_DIR, "chinese_pan_geome.csv"),
  "Human Pan Genome Consortium"            = file.path(DATA_DIR, "hprc2.csv"),
  "2nd Screening Post Assembly of K-mers"  = file.path(DATA_DIR, "2nd_screen_contamination.csv")
)

read_kraken <- function(path, label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- trimws(colnames(df))
  names(df)[names(df) == "%"] <- "pct"
  df$taxName   <- trimws(df$taxName)
  df$rank      <- trimws(df$rank)
  df$Reference <- label
  df
}

# ── Read & clean data ─────────────────────────────────────────
raw <- bind_rows(mapply(read_kraken, FILES, names(FILES), SIMPLIFY = FALSE))

df_all <- raw %>%
  filter(
    !rank    %in% EXCLUDE_RANKS,
    !taxName %in% EXCLUDE_NAMES,
    kmers > 0,
    cov   > 0
  ) %>%
  distinct()                              # remove duplicate rows

# ── NEW: totals for the two extra summary panels ──────────────
# "Unique" stage — contamination k-mers removed from each individual
# pan-genome reference, prior to combining them into the assembled database.
unique_removed_df <- data.frame(
  Reference = c(
    "African Pan Genome",
    "Chinese Pan Genome",
    "Khoe-San Pan Genome",
    "Human Pan Genome Consortium",
    "Arab Pan Genome"
  ),
  Count = c(43323, 4157, 518, 3890, 579),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Count)) %>%
  mutate(Reference = factor(Reference, levels = Reference))

# Two-stage screening total, split by WHEN each pass happened relative to
# assembling the individual pan-genomes into the combined database:
#   1st Screen = PRE-assembly. This is simply the cumulative total of the
#                five "Unique" per-reference counts above (43,323 + 4,157 +
#                518 + 3,890 + 579 = 52,467) — i.e. contamination caught while
#                each reference was still screened on its own.
#   2nd Screen = POST-assembly. A further, independent pass run on the
#                combined/assembled k-mer database, catching contamination
#                the per-reference screens missed (408,365 additional k-mers).
# 52,467 + 408,365 = 460,832 total contamination k-mers removed.
assembled_removed_df <- data.frame(
  Stage = c("1st Screen\n(Pre-Assembly)", "2nd Screen\n(Post-Assembly)"),
  Count = c(52467, 408365),
  stringsAsFactors = FALSE
) %>%
  mutate(Stage = factor(Stage, levels = Stage))

# ── Helper: build a summary bar-chart panel ───────────────────
summary_bar_panel <- function(df, x_col, fill_values, title, subtitle) {
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = Count, fill = .data[[x_col]])) +
    geom_col(width = 0.65, colour = "grey20", linewidth = 0.25) +
    geom_text(
      aes(label = comma(Count)),
      vjust = -0.5, size = 2.4, colour = "black", fontface = "bold"
    ) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_y_log10(
      name   = "K-mers removed (log10)",
      labels = label_comma(),
      expand = expansion(mult = c(0, 0.28))
    ) +
    labs(title = title, subtitle = subtitle, x = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "grey25", linewidth = 0.45),
      axis.text.x        = element_text(colour = "black", size = 6.3, angle = 20, hjust = 1),
      axis.text.y        = element_text(colour = "black", size = 7),
      axis.title         = element_text(colour = "black", size = 7.5),
      plot.title         = element_text(
        size = 8.5, face = "bold", colour = "grey20", margin = margin(b = 1)
      ),
      plot.subtitle      = element_text(
        size = 6.8, face = "italic", colour = "grey40", margin = margin(b = 3)
      ),
      plot.margin        = margin(5, 7, 5, 5)
    )
  
  p
}

p_unique_summary <- summary_bar_panel(
  df          = unique_removed_df,
  x_col       = "Reference",
  fill_values = ref_colours,
  title       = "Unique Screening (Pre-Assembly)",
  subtitle    = "Contamination k-mers removed per reference, before assembly"
)

p_assembled_summary <- summary_bar_panel(
  df          = assembled_removed_df,
  x_col       = "Stage",
  fill_values = c(
    "1st Screen\n(Pre-Assembly)"  = ref_colours[["1st Screening Pre-Assembly of K-mers"]],
    "2nd Screen\n(Post-Assembly)" = ref_colours[["2nd Screening Post Assembly of K-mers"]]
  ),
  title       = "Screening Summary",
  subtitle    = "1st = cumulative total of the five pre-assembly reference\nscreens (left); 2nd = additional pass post-assembly"
)

# ── Build one plot per reference ──────────────────────────────
ref_order <- c(
  "Arab Pan Genome",
  "African Pan Genome",
  "Khoe-San Pan Genome",
  "Chinese Pan Genome",
  "Human Pan Genome Consortium",
  "2nd Screening Post Assembly of K-mers"
)

plot_list <- lapply(ref_order, function(ref) {
  
  df <- df_all %>% filter(Reference == ref)
  
  if (nrow(df) == 0) return(NULL)
  
  # Per-panel coverage threshold — force into local scalar so ggplot
  # captures the value now, not lazily at render time
  cov_thresh <- local({
    if (ref == "2nd Screening Post Assembly of K-mers") COV_THRESH_2ND else COV_THRESH
  })
  
  # Compute excluded flag using the correct threshold for this panel
  df <- df %>%
    mutate(excluded = kmers > KMER_THRESH | cov > cov_thresh)
  
  # Top-N labels: up to 10 excluded taxa ranked by k-mer count
  n_label <- min(10, sum(df$excluded))
  
  label_df <- df %>%
    filter(excluded == TRUE) %>%
    arrange(desc(kmers)) %>%
    slice_head(n = n_label) %>%
    mutate(short_label = short_tax_label(taxName))
  
  df <- df %>%
    mutate(
      to_label = paste(taxName, kmers, cov) %in%
        paste(label_df$taxName, label_df$kmers, label_df$cov),
      label = if_else(to_label, short_tax_label(taxName), NA_character_)
    )
  
  # Axis limits — generous padding on log scale
  x_min <- pmin(min(df$kmers) * 0.3, KMER_THRESH * 0.15)
  x_max <- max(df$kmers) * 12
  
  y_min <- pmin(min(df$cov) * 0.3, cov_thresh * 0.15)
  y_max <- max(df$cov) * 10
  
  pt_colour <- ref_colours[ref]
  
  df_base    <- df[!df$to_label, ]
  df_flagged <- df[df$to_label,  ]
  
  # Zone rect data — use explicit axis limits rather than Inf/-Inf
  # so rects render correctly on log scales
  zones <- data.frame(
    xmin  = c(KMER_THRESH, x_min,       x_min),
    xmax  = c(x_max,       KMER_THRESH, KMER_THRESH),
    ymin  = c(y_min,       cov_thresh,  y_min),
    ymax  = c(y_max,       y_max,       cov_thresh),
    fill  = c("#CC0000",   "#CC0000",   "#1A9850"),
    alpha = c(0.07,        0.07,        0.06),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot(df, aes(x = kmers, y = cov)) +
    
    # ── Background zones ───────────────────────────────────
    geom_rect(
      data = zones,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = fill, alpha = alpha),
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_alpha_identity() +
    
    # ── Threshold lines ────────────────────────────────────
    geom_vline(xintercept = KMER_THRESH, linetype = "dashed",
               colour = "#CC0000", linewidth = 0.5) +
    geom_hline(yintercept = cov_thresh,  linetype = "dashed",
               colour = "#CC0000", linewidth = 0.5) +
    
    # ── Non-labelled points ────────────────────────────────
    geom_point(
      data  = df_base,
      aes(size = log10(kmers + 1)),
      colour = pt_colour,
      alpha  = 0.45,
      shape  = 16
    ) +
    
    # ── Labelled (excluded) points ─────────────────────────
    geom_point(
      data   = df_flagged,
      aes(size = log10(kmers + 1)),
      shape  = 21, fill = NA,
      colour = "black", stroke = 0.9, alpha = 1
    ) +
    geom_point(
      data   = df_flagged,
      aes(size = log10(kmers + 1)),
      colour = pt_colour,
      shape  = 16, alpha = 0.95
    ) +
    
    # ── Labels ─────────────────────────────────────────────
    geom_label_repel(
      data               = df_flagged,
      aes(label          = label),
      colour             = "black",
      fill               = alpha("white", 0.88),
      size               = 1.9,
      fontface           = "italic",
      label.size         = 0.12,
      label.padding      = unit(0.08, "lines"),
      segment.colour     = "grey45",
      segment.size       = 0.22,
      segment.linetype   = "solid",
      min.segment.length = 0,
      box.padding        = 0.75,
      point.padding      = 0.30,
      force              = 0.15,
      force_pull         = 0.12,
      max.overlaps       = Inf,
      direction          = "both",
      seed               = 42
    ) +
    
    # ── Zone annotations ───────────────────────────────────
    annotate("text",
             x = x_min * 1.8, y = y_min * 1.8,
             label = "Retained", hjust = 0, vjust = 0,
             size = 2.2, colour = "grey55"
    ) +
    annotate("text",
             x = KMER_THRESH * 1.2, y = y_min * 1.8,
             label = "Excluded", hjust = 0, vjust = 0,
             size = 2.2, colour = "#CC0000", fontface = "bold"
    ) +
    
    # ── Threshold annotations ──────────────────────────────
    annotate("text",
             x = KMER_THRESH * 1.08, y = y_max * 0.40,
             label = "k-mers > 100", hjust = 0, vjust = 1,
             size = 1.8, colour = "#CC0000"
    ) +
    annotate("text",
             x = x_min * 1.4, y = cov_thresh * 1.35,
             label = paste0("cov > ", cov_thresh), hjust = 0, vjust = 0,
             size = 1.8, colour = "#CC0000"
    ) +
    
    # ── Scales ─────────────────────────────────────────────
    scale_x_log10(
      name   = "Classified k-mers detected",
      labels = label_comma(),
      limits = c(x_min, x_max),
      oob    = scales::squish
    ) +
    scale_y_log10(
      name   = "Genome coverage",
      labels = label_scientific(digits = 1),
      limits = c(y_min, y_max),
      oob    = scales::squish
    ) +
    scale_size_continuous(range = c(0.5, 4.0), guide = "none") +
    
    # ── Panel title (+ subtitle for 2nd Screening) ────────
    labs(
      title    = ref_labels[ref],
      subtitle = if (ref == "2nd Screening Post Assembly of K-mers")
        "Coverage threshold refined to 0.005"
      else
        NULL
    ) +
    
    # ── Theme ──────────────────────────────────────────────
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major  = element_line(colour = "grey93", linewidth = 0.25),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "grey25", linewidth = 0.45),
      axis.text         = element_text(colour = "black", size = 7),
      axis.title        = element_text(colour = "black", size = 7.5),
      plot.title        = element_text(
        size = 8.5, face = "bold",
        colour = pt_colour,
        margin = margin(b = 1)
      ),
      plot.subtitle     = element_text(
        size = 6.8, face = "italic",
        colour = "grey40",
        margin = margin(b = 3)
      ),
      plot.margin       = margin(5, 7, 5, 5)
    )
  
  p
})
names(plot_list) <- ref_order

# ── Assemble grid: 2 rows of taxon-level panels + 1 row of summaries ──
final_plot <-
  (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]] | plot_list[[6]]) /
  (p_unique_summary | p_assembled_summary) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    caption = paste0(
      "Rows 1–2: labelled points (○) show the top 10 excluded taxa per reference, ranked by classified k-mer count. ",
      "Thresholds: k-mers > ", KMER_THRESH, " or genome coverage > ", COV_THRESH,
      " (2nd Screening: coverage > 0.005). ",
      "Row 3 (left): total contamination k-mers removed per reference during unique, pre-assembly screening. ",
      "Row 3 (right): the same pre-assembly total re-shown as \"1st Screen\", alongside \"2nd Screen\" — ",
      "a further pass run after the references were combined into the assembled k-mer database."
    ),
    theme = theme(
      plot.caption = element_text(size = 6.5, colour = "grey40", hjust = 0)
    )
  )

# ── Save ──────────────────────────────────────────────────────
ggsave(
  OUTPUT_PATH,
  final_plot,
  width  = 270,
  height = 300,   # increased from 200mm to fit the extra summary row
  units  = "mm",
  dpi    = 800,
  device = cairo_pdf
)
ggsave(
  sub("\\.pdf$", ".png", OUTPUT_PATH),
  final_plot,
  width  = 270,
  height = 300,
  units  = "mm",
  dpi    = 800
)

cat("Saved:", OUTPUT_PATH, "\n")
cat("Saved:", sub("\\.pdf$", ".png", OUTPUT_PATH), "\n")

# ── Summary of labelled taxa per panel ────────────────────────
cat("\n── Labelled taxa per reference ──\n")
df_all %>%
  mutate(
    cov_thresh_local = if_else(
      Reference == "2nd Screening Post Assembly of K-mers",
      COV_THRESH_2ND, COV_THRESH
    ),
    excluded = kmers > KMER_THRESH | cov > cov_thresh_local
  ) %>%
  filter(excluded == TRUE) %>%
  group_by(Reference) %>%
  arrange(desc(kmers), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  select(Reference, taxName, rank, kmers, cov) %>%
  as.data.frame() %>%
  print()

# ── Summary of the two new totals panels ──────────────────────
cat("\n── Unique Screening totals per reference ──\n")
print(unique_removed_df)
cat("\n── Assembled Screening totals per pass ──\n")
print(assembled_removed_df)

print(final_plot)