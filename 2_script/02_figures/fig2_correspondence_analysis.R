suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
})

### ============= Panel A: 2D CoA, updated with different files ========

# Repeat this block for each subtypes (H5N1,H5N2,H5N6,H5N8,Others)
df  <- read.csv(file = "data/RSCU/COA_CG.csv", header = TRUE,
                na.strings = c("", " ", "N/A", "NULL"))
df1 <- read.csv(file = "data/RSCU/EIGEN_CG.csv", header = FALSE)
df2 <- read.csv(file = "data/RSCU/CODON_CG.csv", header = TRUE)

df$Segment <- factor(df$Segment,
                     levels = c("HA","M1","NA","NP","NS1","PA","PB1","PB2"))
df$Subtype <- factor(df$Subtype, levels = c("N1","N2","N6","N8","Others"))

subtype_color <- c("N1"="#953638", "N2"="#1c516c", "N6"="#d8cd9d",
                   "N8"="#fda649", "Others"="#C5C6C7")
gene_color <- c("HA"="#c14f60","M1"="#ff94b2","NA"="#f0a529","NP"="#21c5c4",
                "NS1"="#C8C2D6","PA"="#8FBC8f","PB1"="#1d4e08","PB2"="#005cb8")

df2$vector_length <- sqrt(df2$Axis1^2 + df2$Axis2^2)
top_codons <- df2[order(-df2$vector_length), ][1:15, ]

p1 <- ggplot() +
  geom_point(data = df, aes(x = Axis1, y = Axis2, color = Subtype, fill = Subtype),
             size = 2, alpha = 1) +
  stat_ellipse(data = df, aes(x = Axis1, y = Axis2, color = Subtype, fill = Subtype),
               geom = "polygon", alpha = 0.1, level = 0.95) +
  geom_segment(data = top_codons,
               aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(.15, "cm"), type = "closed"),
               color = "grey", alpha = .6, size = 1) +
  geom_text_repel(data = top_codons,
                  aes(x = Axis1, y = Axis2, label = label),
                  inherit.aes = FALSE,
                  color = "black", size = 4,
                  max.overlaps = 20,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  box.padding = 0.3,
                  point.padding = 0.3) +
  theme_test(base_rect_size = 2) +
  geom_hline(yintercept = 0.0, color = "red", linewidth = 1.0, linetype = "dashed") +
  geom_vline(xintercept = 0.0, color = "red", linewidth = 1.0, linetype = "dashed") +
  xlab(paste0("Axis1 (", round(df1[1, 3], 4) * 100, "%)")) +
  ylab(paste0("Axis2 (", round(df1[2, 3], 4) * 100, "%)")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-.2, .2, 0.05), limits = c(-.2, .2)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(-.12, .12, 0.04), limits = c(-.12, .12)) +
  scale_color_manual(values = subtype_color) +
  scale_fill_manual(values = subtype_color) +
  ggtitle("H5Nx Complete Genome") +
  theme(plot.title = element_text(size = 26, hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text  = element_text(size = 18),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.15, "cm"),
        plot.margin = unit(c(.5, 1.2, .5, .5), "cm"),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.2, "cm"),
        legend.margin = margin(l = 20, unit = "pt"))

print(p2)

ggsave(filename = "COA_CG.jpg", plot = p1, width = 8, height = 6,
       device = "jpg", dpi = 300)

### ==================== Panel B: 3D CoA per segment (PB2 example) =================

suppressPackageStartupMessages({
  library(plotly)
})

# Repeat this block for each segment (HA, M1, NA, NP, NS1, PA, PB1, PB2),

df  <- read.csv(file = "data/RSCU/COA_PB2.csv", header = TRUE)
df1 <- read.csv(file = "data/RSCU/EIGEN_PB2.csv", header = FALSE)

df <- df %>% mutate(color_group = Subtype, shape_group = Host)

p3d <- plot_ly(df) %>%
  add_trace(
    x = ~Axis1, y = ~Axis2, z = ~Axis3,
    color = ~color_group,
    colors = subtype_color,
    symbol = ~shape_group,
    symbols = c("Human" = "square", "Avian" = "circle", "Nonhuman Mammal" = "diamond"),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 7, opacity = .8, line = list(width = .5, color = "black"))
  ) %>%
  layout(
    scene = list(
      # Axis percentages are read from EIGEN_PB2.csv (rows 1-3, inertia column)
      xaxis = list(title = list(text = "Axis 1 (28.25%)", font = list(size = 22))),
      yaxis = list(title = list(text = "Axis 2 (22.31%)", font = list(size = 22))),
      zaxis = list(title = list(text = "Axis 3 (10.21%)", font = list(size = 22))),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
      backgroundcolor = "rgb(230,230,230)",
      gridcolor = "rgb(255,255,255)",
      showbackground = TRUE
    ),
    legend = list(
      title = list(text = "<b>Subtype</b>", font = list(size = 26, family = "Arial")),
      font = list(size = 24, family = "Arial"),
      itemsizing = "constant"
    )
  )

print(p3d)