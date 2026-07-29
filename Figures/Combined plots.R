# ============================================================
#  ARKS – Combined Nature Figure
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggrepel)
  library(dplyr)
  library(scales)
  library(cowplot)
})

out_path <- "/Users/hce24xau/Desktop/Genomics_Benchmarking/Map, Unique, Distinct/Unique_distinct_Map.pdf"

pal <- c(
  "GRCh37"              = "#E41A1C",
  "GRCh38"              = "#FF7F00",
  "IPD-IMGT/HLA"        = "#A65628",
  "T2T"                 = "#F781BF",
  "African Pan"         = "#009E73",
  "Arab Pan Genome"     = "#377EB8",
  "Khoe-San Pan Genome" = "#999999",
  "Chinese Pan Genome"  = "#C45000",
  "HPRC Pan Genome"     = "#7B2D8B"
)

# Ordered by ascending distinct k-mers
dataset_order <- c(
  "IPD-IMGT/HLA",
  "African Pan",
  "GRCh37",
  "GRCh38",
  "T2T",
  "Khoe-San Pan Genome",
  "Chinese Pan Genome",
  "Arab Pan Genome",
  "HPRC Pan Genome"
)

unique_tbl <- data.frame(
  Dataset = factor(dataset_order, levels = dataset_order),
  Unique  = c(
    488497,
    3755046,
    291533,
    84476,
    403,
    28674511,
    178268913,
    111221570,
    681691121
  )
)

distinct_tbl <- data.frame(
  Dataset  = factor(dataset_order, levels = dataset_order),
  Distinct = c(
    543558,
    17048468,
    2490864048,
    2505641687,
    2512377669,
    2682436250,
    3003063917,
    3056369047,
    3709552761
  )
)

fold_tbl <- data.frame(
  Dataset = factor(dataset_order, levels = dataset_order),
  Fold = c(
    543558 / 488497,
    17048468 / 3755046,
    2490864048 / 291533,
    2505641687 / 84476,
    2512377669 / 403,
    2682436250 / 28674511,
    3003063917 / 178268913,
    3056369047 / 111221570,
    3709552761 / 681691121
  )
)

bar_theme <- theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    axis.text.y        = element_text(size = 6, hjust = 1),
    axis.text.x        = element_text(size = 5.5),
    axis.title.x       = element_text(size = 6.5),
    axis.title.y       = element_blank(),
    plot.title         = element_text(hjust = 0, size = 7, face = "bold",
                                      margin = margin(b = 3)),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(4, 8, 4, 4)
  )

p_distinct <- ggplot(distinct_tbl,
                     aes(x = Distinct, y = Dataset, fill = Dataset)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = comma(Distinct)), hjust = -0.08, size = 1.9) +
  scale_x_continuous(
    trans = "log10",
    labels = label_number(accuracy = 1, scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.32))
  ) +
  scale_y_discrete(limits = rev(dataset_order)) +
  scale_fill_manual(values = pal) +
  labs(title = "a", x = expression("Distinct 31-mers ("*log[10]*")")) +
  bar_theme

p_unique <- ggplot(unique_tbl,
                   aes(x = Unique, y = Dataset, fill = Dataset)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = comma(Unique)), hjust = -0.08, size = 1.9) +
  scale_x_continuous(
    trans = "log10",
    labels = label_number(accuracy = 1, scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.32))
  ) +
  scale_y_discrete(limits = rev(dataset_order)) +
  scale_fill_manual(values = pal) +
  labs(title = "b", x = expression("Unique 31-mers ("*log[10]*")")) +
  bar_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_fold <- ggplot(fold_tbl,
                 aes(x = Fold, y = Dataset, fill = Dataset)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = ifelse(
      Fold >= 1000,
      paste0(format(round(Fold), big.mark = ","), "x"),
      paste0(round(Fold, 1), "x")
    )),
    hjust = -0.08,
    size = 1.9
  ) +
  scale_x_continuous(
    trans = "log10",
    labels = label_number(accuracy = 1, scale_cut = cut_short_scale(), suffix = "x"),
    expand = expansion(mult = c(0, 0.38))
  ) +
  scale_y_discrete(limits = rev(dataset_order)) +
  scale_fill_manual(values = pal) +
  labs(title = "c", x = expression("k-mer redundancy ("*log[10]*")")) +
  bar_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

top_row <- plot_grid(
  p_distinct,
  p_unique,
  p_fold,
  ncol = 3,
  align = "h",
  axis = "tb",
  rel_widths = c(1.35, 1, 1)
)

world <- ne_countries(scale = "medium", returnclass = "sf")

map_pal <- c(
  "African Pan-genome"  = "#009E73",
  "HPRC Pan-Genome"     = "#7B2D8B",
  "Chinese Pan-Genome"  = "#C45000",
  "Khoe-San Pan-Genome" = "#D4AF37",
  "Arab Pan-Genome"     = "#1F78B4"
)

cpc_df <- data.frame(
  Label = c(
    "Oroqen","Evenki","Hezhen","Man","Daur","Mongol","Chosen","Kazakh",
    "Kyrgyz","Uyghur","Yugur","Tu","Salar","Hui","Tibetan","Qiang","Tujia",
    "Bai","Drung","Lisu","Mosuo","Naxi","Yi","Achang","Jingpo","Deang",
    "Khatso","Jino","Blang","Wa","Han (North)","Han (Central)","Han (South)",
    "Miao","She","Yao","Bouyei","Dong","Zhuang","Kinh"
  ),
  Latitude  = c(50,50,48,42,49,47,39,45,42,43,40,37,37,34,31,32,30,26,27,27,
                27,27,26,25,25,23,25,23,22,22,39,34,27,26,26,25,26,26,24,21),
  Longitude = c(124,123,133,123,124,116,126,85,75,80,97,100,101,106,90,103,110,
                100,98,99,99,99,102,98,97,99,101,100,100,99,114,112,113,109,118,
                110,107,109,107,105),
  Dataset   = "Chinese Pan-Genome"
)

hprc_df <- data.frame(
  Label = c(
    "Puerto Rican","Gambian\n(Mandinka)","Caribbean\n(Barbados)","Colombian",
    "Peruvian","Mende\n(Sierra Leone)","Han Chinese\n(South)","Yoruba\n(Nigeria)",
    "Kinh\n(Vietnam)","Esan\n(Nigeria)","African Ancestry\n(SW USA)","Maasai\n(Kenya)",
    "Japanese","Bengali\n(Bangladesh)","British\n(England)","Pakistan",
    "Finnish","Toscani\n(Italy)","Gujarati\n(Houston)","Luhya\n(Kenya)",
    "Chinese Dai\n(Yunnan)","Sri Lankan Tamil\n(UK)","Han Chinese\n(Beijing)",
    "Telugu\n(UK)","Iberian\n(Spain)","Mexican Ancestry\n(Los Angeles)",
    "African American\n(St. Louis)"
  ),
  Latitude  = c(18.22,13.45,13.10,6.25,-12.04,7.87,23.13,7.38,10.82,6.75,35.00,-1.83,
                35.68,23.68,53.00,31.52,61.92,43.77,29.76,0.61,21.86,51.50,39.90,52.48,
                40.46,34.05,38.63),
  Longitude = c(-66.59,-16.57,-59.61,-75.56,-77.03,-11.92,113.26,3.93,106.63,6.25,
                -106.00,36.85,139.76,90.35,-1.50,74.36,25.75,11.25,-95.36,34.77,
                100.80,-0.12,116.40,-1.89,-3.75,-118.24,-90.20),
  Dataset   = "HPRC Pan-Genome"
)

african_df <- data.frame(
  Label = c(
    "African American\n(Atlanta)","African American\n(Baltimore/DC)",
    "African American\n(Chicago)","African American\n(Detroit)",
    "African American\n(Jackson, MS)","African American\n(Nashville)",
    "African American\n(New York)","African American\n(San Francisco)",
    "African American\n(Winston-Salem)","Caribbean\n(Barbados)",
    "Brazilian","Colombian","Dominican","Gabonese","Honduran",
    "Jamaican","Palenquero\n(Colombia)","Yoruba\n(Nigeria)","Puerto Rican"
  ),
  Latitude  = c(33.75,39.00,41.88,42.33,32.30,36.17,40.71,37.77,36.10,
                13.10,-14.23,6.25,18.74,-0.80,14.08,18.01,10.40,9.08,18.22),
  Longitude = c(-84.39,-76.61,-87.63,-83.05,-90.18,-86.78,-74.00,-122.42,-80.24,
                -59.61,-51.93,-75.56,-70.16,11.61,-87.21,-76.79,-74.72,8.68,-66.59),
  Dataset   = "African Pan-genome"
)

khoesan_df <- data.frame(
  Label     = "Ju|'hoansi\n(Nyae Nyae, Namibia)",
  Latitude  = -19.17,
  Longitude = 20.80,
  Dataset   = "Khoe-San Pan-Genome"
)

arab_df <- data.frame(
  Label = c("UAE","Saudi Arabia","Oman","Jordan","Egypt","Morocco","Syria","Yemen"),
  Latitude  = c(24.47, 24.71, 23.59, 31.95, 30.04, 34.02, 33.51, 15.37),
  Longitude = c(54.37, 46.68, 58.41, 35.93, 31.24, -6.84, 36.28, 44.19),
  Dataset   = "Arab Pan-Genome"
)

all_df   <- bind_rows(cpc_df, hprc_df, african_df, khoesan_df, arab_df)
label_df <- all_df |> distinct(Latitude, Longitude, .keep_all = TRUE)

us_east_labels <- c(
  "African American\n(New York)","African American\n(Baltimore/DC)",
  "African American\n(Detroit)","African American\n(Chicago)",
  "African American\n(Atlanta)","African American\n(Nashville)",
  "African American\n(Winston-Salem)","African American\n(Jackson, MS)",
  "African American\n(St. Louis)"
)

us_west_labels <- c(
  "African American\n(San Francisco)","African Ancestry\n(SW USA)",
  "Mexican Ancestry\n(Los Angeles)"
)

caribbean_labels <- c(
  "Puerto Rican","Dominican","Jamaican","Caribbean\n(Barbados)",
  "Honduran","Palenquero\n(Colombia)","Colombian","Gujarati\n(Houston)"
)

south_am_labels <- c("Peruvian","Brazilian")
mena_labels <- c("UAE","Saudi Arabia","Oman","Jordan","Egypt","Morocco","Syria","Yemen")
americas_all <- c(us_east_labels, us_west_labels, caribbean_labels, south_am_labels)

label_other   <- label_df |> filter(!Label %in% c(americas_all, mena_labels))
label_us_east <- label_df |> filter(Label %in% us_east_labels)
label_us_west <- label_df |> filter(Label %in% us_west_labels)
label_carib   <- label_df |> filter(Label %in% caribbean_labels)
label_south   <- label_df |> filter(Label %in% south_am_labels)
label_mena    <- label_df |> filter(Label %in% mena_labels)

repel_base <- list(
  size = 1.6,
  lineheight = 0.82,
  fontface = "bold",
  family = "sans",
  segment.size = 0.22,
  segment.color = "grey40",
  segment.alpha = 0.65,
  min.segment.length = 0.12,
  show.legend = FALSE
)

p_map <- ggplot() +
  geom_sf(data = world, fill = "#F7F7F7", colour = "grey65", linewidth = 0.18) +
  geom_point(
    data = all_df,
    aes(x = Longitude, y = Latitude, fill = Dataset),
    shape = 21,
    size = 2.2,
    colour = "white",
    stroke = 0.32,
    alpha = 0.88
  ) +
  do.call(geom_text_repel, c(list(
    data = label_other,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    max.overlaps = 200,
    force = 10,
    force_pull = 0.5,
    box.padding = 0.45,
    point.padding = 0.25,
    max.time = 8
  ), repel_base)) +
  do.call(geom_text_repel, c(list(
    data = label_us_east,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    nudge_x = -28,
    direction = "y",
    force = 8,
    force_pull = 0.05,
    box.padding = 0.7,
    point.padding = 0.25,
    max.overlaps = 100
  ), repel_base)) +
  do.call(geom_text_repel, c(list(
    data = label_us_west,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    nudge_x = -20,
    direction = "y",
    force = 4,
    force_pull = 0.1,
    box.padding = 0.6,
    point.padding = 0.25,
    max.overlaps = 50
  ), repel_base)) +
  do.call(geom_text_repel, c(list(
    data = label_carib,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    nudge_x = 10,
    nudge_y = -5,
    direction = "y",
    force = 5,
    force_pull = 0.2,
    box.padding = 0.55,
    point.padding = 0.25,
    max.overlaps = 80
  ), repel_base)) +
  do.call(geom_text_repel, c(list(
    data = label_south,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    nudge_x = -12,
    nudge_y = -5,
    force = 2
  ), repel_base)) +
  do.call(geom_text_repel, c(list(
    data = label_mena,
    mapping = aes(x = Longitude, y = Latitude, label = Label, colour = Dataset),
    nudge_x = 14,
    direction = "y",
    force = 6,
    force_pull = 0.15,
    box.padding = 0.6,
    point.padding = 0.25,
    max.overlaps = 80
  ), repel_base)) +
  scale_fill_manual(
    values = map_pal,
    name = "Reference dataset",
    breaks = c(
      "African Pan-genome",
      "HPRC Pan-Genome",
      "Chinese Pan-Genome",
      "Khoe-San Pan-Genome",
      "Arab Pan-Genome"
    ),
    labels = c(
      expression("African Pan-genome (" * italic(n) * " = 910)"),
      expression("HPRC Pan-Genome (" * italic(n) * " = 236)"),
      expression("Chinese Pan-Genome (" * italic(n) * " = 58)"),
      expression("Khoe-San Pan-Genome (" * italic(n) * " = 3)"),
      expression("Arab Pan-Genome (" * italic(n) * " = 53)")
    )
  ) +
  scale_colour_manual(values = map_pal, guide = "none") +
  coord_sf(xlim = c(-160, 160), ylim = c(-55, 72), expand = FALSE) +
  labs(title = "d") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.title      = element_text(hjust = 0, size = 7, face = "bold",
                                   margin = margin(b = 3, t = 2)),
    legend.position = "bottom",
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 6, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(2, 4, 2, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5), nrow = 2))

full_figure <- plot_grid(
  top_row,
  p_map,
  ncol = 1,
  rel_heights = c(0.36, 0.64)
)

print(full_figure)

ggsave(
  filename = out_path,
  plot = full_figure,
  width = 220,
  height = 170,
  units = "mm",
  dpi = 1200,
  bg = "white"
)

message("✅ Saved: ", out_path)