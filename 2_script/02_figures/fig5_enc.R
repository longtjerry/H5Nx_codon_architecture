suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

segments <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

# Expected ENC under GC mutation pressure alone
f <- function(x) 2 + x + 29 / (x^2 + (1 - x)^2)
x  <- seq(0, 1.0, 0.001)
df1 <- data.frame(x, y = f(x))

subtype_color <- c("H5N1" = "#953638", "H5N2" = "#1c516c", "H5N6" = "#d8cd9d",
                   "H5N8" = "#fda649", "Others" = "#C5C6C7")

host_shape <- c("Human" = 22, "Avian" = 21, "Nonhuman Mammal" = 23)

plot_list <- list()

for (segment in segments) {
  filename <- paste0("data/ENC_GC/ENC_", segment, ".csv")
  df <- read.csv(file = filename, header = TRUE)
  
  summary(lm(Nc ~ GC3s, data = df))
  
  df$Subtype <- factor(df$Subtype, levels = c("H5N1", "H5N2", "H5N6", "H5N8", "Others"))
  
  p <- ggplot(df, aes(x = GC3s, y = Nc, fill = Subtype, shape = Host)) +
    geom_point(color = "black", alpha = 0.8, size = 1.5, stroke = 0.3) +
    geom_line(data = df1, aes(x, y), color = "red", size = 0.8, inherit.aes = FALSE) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(30, 60, 10), limits = c(30, 61)) +
    scale_fill_manual(values = subtype_color) +
    scale_shape_manual(values = host_shape) +
    ggtitle(segment) +
    theme(plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title   = element_text(size = 9),
          axis.text    = element_text(size = 7),
          legend.position = "none") 
  
  plot_list[[segment]] <- p
}

combined_plot <- plot_list[["HA"]] + plot_list[["M1"]] + plot_list[["NA"]] + plot_list[["NP"]] +
  plot_list[["NS1"]] + plot_list[["PA"]] + plot_list[["PB1"]] + plot_list[["PB2"]] +
  plot_layout(ncol = 4, nrow = 2) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))

legend_fill_df <- data.frame(
  Subtype = factor(names(subtype_color), levels = names(subtype_color)), x = 1, y = 1)
legend_fill_plot <- ggplot(legend_fill_df, aes(x = x, y = y, fill = Subtype)) +
  geom_point(shape = 21, color = "black", size = 5, stroke = 0.6) +
  scale_fill_manual(values = subtype_color, name = "Subtype") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"))

legend_shape_df <- data.frame(
  Host = factor(names(host_shape), levels = names(host_shape)), x = 1, y = 1)
legend_shape_plot <- ggplot(legend_shape_df, aes(x = x, y = y, shape = Host)) +
  geom_point(fill = "gray50", color = "black", size = 5, stroke = 0.6) +
  scale_shape_manual(values = host_shape, name = "Host") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"))

get_legend_grob <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) return(tmp$grobs[[leg]]) else return(NULL)
}

legend_fill  <- get_legend_grob(legend_fill_plot)
legend_shape <- get_legend_grob(legend_shape_plot)

if (!is.null(legend_fill) && !is.null(legend_shape)) {
  legend_combined <- wrap_elements(full = legend_fill) + wrap_elements(full = legend_shape) +
    plot_layout(ncol = 2)
} else {
  legend_combined <- NULL
}

if (!is.null(legend_combined)) {
  final_plot <- combined_plot / legend_combined + plot_layout(heights = c(10, 1))
} else {
  final_plot <- combined_plot
}

print(final_plot)

### ============== Figure S2A: whole-genome ENC-GC3s ==========

df2 <- read.csv("data/ENC_GC/ENC_CG.csv", header = TRUE)

p1 <- ggplot(df2, aes(x = GC3s, y = Nc, fill = Subtype, shape = Host)) +
  geom_point(color = "black", alpha = 0.8, size = 3, stroke = 0.3) +
  geom_line(data = df1, aes(x, y), color = "red", size = 0.8, inherit.aes = FALSE) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0.0, 1, 0.2), limits = c(0.0, 1.0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(30, 60, 5), limits = c(30, 61)) +
  scale_fill_manual(values = subtype_color) +
  scale_shape_manual(values = host_shape) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  ggtitle("Complete Genome") +
  theme(plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        legend.text = element_text(size = 15))

print(p1)

p_mini <- ggplot(df2, aes(x = GC3s, y = Nc, fill = Subtype, shape = Host)) +
  geom_point(alpha = 0.8, size = 3) +
  theme_bw() +
  scale_fill_manual(values = subtype_color) +
  scale_shape_manual(values = host_shape) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  theme(legend.position = "none",
        panel.grid   = element_blank(),
        axis.title   = element_blank(),
        axis.ticks   = element_blank(),
        axis.text    = element_blank())

print(p_mini)


ggsave(filename = "ENC_segments.jpg", plot = final_plot,
       width = 14, height = 8, device = "jpg", dpi = 720)
ggsave(filename = "ENC_CG.jpg", plot = p1,
       width = 6, height = 4, device = "jpg", dpi = 720)
ggsave(filename = "ENC_CG_mini.jpg", plot = p_mini,
       width = 4, height = 3, device = "jpg", dpi = 720)