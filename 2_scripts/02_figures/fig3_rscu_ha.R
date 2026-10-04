suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

### =============== Panel A: RSCU heatmap ============

df <- read.csv(file = "data/RSCU/heatmap.csv", header = TRUE, row.names = 1)

annotation_col <- data.frame(
  Amino.acid = factor(rep(c("Phe-F", "Leu-L", "Ile-I", "Val-V", "Ser-S", "Pro-P",
                            "Thr-T", "Ala-A", "Tyr-Y", "His-H", "Gln-Q", "Asn-N",
                            "Lys-K", "Asp-D", "Glu-E", "Cys-C", "Arg-R", "Gly-G"),
                          times = c(2, 6, 3, 4, 6, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 6, 4)))
)
rownames(annotation_col) <- colnames(df)

ann_colors <- list(
  Amino.acid = c(
    "Phe-F" = "#F961C6", "Leu-L" = "#118AB2", "Ile-I" = "#5E548E",
    "Val-V" = "#7209B7", "Ser-S" = "#073B4C", "Pro-P" = "#9A6198",
    "Thr-T" = "#FFE6AA", "Ala-A" = "#FA7921", "Tyr-Y" = "#86A8E7",
    "His-H" = "#FF6F59", "Gln-Q" = "#06D6A0", "Asn-N" = "#E36414",
    "Lys-K" = "#F9C846", "Asp-D" = "#FF99C8", "Glu-E" = "#577590",
    "Cys-C" = "#3DA35D", "Arg-R" = "#EF476F", "Gly-G" = "#6A9488")
)

p1 <- pheatmap(df, color = colorRampPalette(brewer.pal(n = 7, name = "YlGnBu"))(50),
               cluster_rows = FALSE, cluster_cols = FALSE,
               show_rownames = TRUE,
               display_numbers = FALSE,
               gaps_row = c(5, 10, 15, 20, 25, 30, 35),
               gaps_col = c(2, 8, 11, 15, 21, 25, 29, 33, 35, 37, 39, 41, 43, 45, 47, 49, 55),
               annotation_col = annotation_col, angle_col = "45",
               annotation_colors = ann_colors,
               annotation_legend = TRUE, fontsize = 8, legend = TRUE)
print(p1)

ggsave(filename = "H5Nx_RSCU.jpg", plot = p1, width = 10, height = 5, dpi = 300)

### ============= Panel B: arginine codon usage across HA ====================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

df <- read.csv("data/RSCU/h5_ha_arg_by_subtype.csv")

desired_order <- c("H5N1", "H5N2", "H5N6", "H5N8", "Others")
df$Subtype <- factor(df$Subtype, levels = desired_order)

total_length <- 566

all_positions <- 1:total_length
all_subtypes  <- levels(df$Subtype)

full_grid <- expand.grid(
  AA_Position = all_positions,
  Subtype = all_subtypes,
  stringsAsFactors = FALSE
)

heatmap_data <- full_grid %>%
  left_join(df %>% dplyr::select(AA_Position, Subtype, AGA_AGG_Ratio),
            by = c("AA_Position", "Subtype"))
heatmap_data$Subtype <- factor(heatmap_data$Subtype, levels = desired_order)

domains <- data.frame(
  domain = c("Signal peptide", "HA1", "RBD", "Cleavage site", "HA2"),
  start  = c(1, 17, 124, 340, 347),
  end    = c(17, 340, 277, 347, total_length)
)

domain_colors <- c(
  "Signal peptide" = "#bdbdbd",
  "HA1"            = "#99d8c9",
  "RBD"            = "#2ca25f",
  "Cleavage site"  = "#fdbf6f",
  "HA2"            = "#fb9a99"
)

# Top panel: HA domain schematic
x_pad <- 20
x_limits <- c(0 - x_pad, total_length + x_pad)

p_schematic <- ggplot() +
  geom_rect(data = domains,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = domain),
            color = "black", alpha = 0.8, linewidth = 0.3) +
  scale_fill_manual(values = domain_colors) +
  geom_text(data = domains,
            aes(x = (start + end) / 2, y = 0.5, label = domain),
            size = 6, fontface = "bold") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0), name = "") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = NULL, name = "") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 10, 0, 10))

# Bottom panel: AGA/AGG ratio heatmap by subtype
subtype_labels <- data.frame(
  Subtype = factor(desired_order, levels = desired_order),
  x = total_length + 3,
  y = desired_order,
  label = desired_order
)

p_heatmap <- ggplot(heatmap_data, aes(x = AA_Position, y = Subtype, fill = AGA_AGG_Ratio)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(data = subtype_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, hjust = 0, vjust = .5) +
  scale_fill_gradient2(
    low = "#053061", mid = "#ffffbf", high = "#d73027",
    midpoint = 2/6,       
    limits = c(0, 1),
    na.value = "lightgrey",
    name = NULL,
    breaks = c(0, 0.33, 0.5, 0.75, 1),
    labels = c("0 (CGN only)", "0.33 (random)", "0.5", "0.75", "1 (AGA/AGG only)")
  ) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0), name = "Amino Acid Position") +
  scale_y_discrete(limits = rev(desired_order)) +
  theme_void(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.title.x = element_text(size = 15),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal",
        plot.margin = margin(0, 30, 15, 15))

combined <- p_schematic / p_heatmap + plot_layout(heights = c(1, 2.5))
print(combined)

ggsave("HA_schematic_heatmap_subtypes.jpg", combined, width = 18, height = 4, dpi = 300)