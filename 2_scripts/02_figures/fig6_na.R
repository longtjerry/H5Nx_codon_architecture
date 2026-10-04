suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(dplyr)
})

subtype_color <- c("H5N1" = "#953638", "H5N2" = "#1c516c", "H5N6" = "#d8cd9d",
                   "H5N8" = "#fda649", "Others" = "#C5C6C7")


### ============ Figure S2B: whole-genome neutrality plot =======

df <- read.csv("data/ENC_GC/Neutrality_CG.csv", header = TRUE)
df$Subtype <- factor(df$Subtype, levels = c("H5N1", "H5N2", "H5N6", "H5N8", "Others"))

n_subtypes <- length(levels(df$Subtype))

p_CG <- ggplot(df, aes(x = GC3s, y = GC12)) +
  geom_point(aes(fill = Subtype), color = "black", shape = 21,
             alpha = 0.8, size = 2, stroke = 0.3) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
              color = "red", linetype = "dashed", linewidth = 1.2, fullrange = TRUE) +
  geom_smooth(aes(color = Subtype), method = "lm", se = FALSE,
              formula = y ~ x, linewidth = 0.8, fullrange = TRUE) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(30, 60, 5), limits = c(30, 60)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(30, 60, 5), limits = c(30, 60)) +
  scale_fill_manual(values = subtype_color, name = "Subtype") +
  scale_color_manual(values = subtype_color, name = "Subtype") +
  stat_regline_equation(
    mapping = aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = y ~ x, color = "red", size = 3,
    label.x = 32, label.y = 58, show.legend = FALSE) +
  stat_regline_equation(
    aes(color = Subtype,
        label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = y ~ x, size = 3,
    label.x = 38,
    label.y = seq(40, 32, length.out = n_subtypes),
    show.legend = FALSE) +
  theme_classic() +
  labs(title = "Complete genome", x = "GC3s", y = "GC12") +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        legend.position = "right",
        legend.text = element_text(size = 10))

print(p_CG)

### ============= Figure 6: per-segment neutrality plots ==============

segments <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

process_gene <- function(segment) {
  filename <- paste0("data/ENC_GC/Neutrality_", segment, ".csv")
  df <- read.csv(file = filename, header = TRUE)
  df$Subtype <- factor(df$Subtype, levels = c("H5N1", "H5N2", "H5N6", "H5N8", "Others"))
  
  n_subtypes <- length(unique(df$Subtype))
  
  p <- ggplot(df, aes(x = GC3s, y = GC12)) +
    geom_point(aes(fill = Subtype), color = "black", shape = 21,
               alpha = 0.8, size = 2, stroke = 0.3) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
                color = "red", linetype = "dashed", size = 1.5, fullrange = TRUE) +
    geom_smooth(aes(color = Subtype), method = "lm", se = FALSE,
                formula = y ~ x, size = 0.8, fullrange = TRUE) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(30, 60, 5), limits = c(30, 60)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(30, 60, 5), limits = c(30, 60)) +
    scale_fill_manual(values = subtype_color, name = "Subtype") +
    scale_color_manual(values = subtype_color, name = "Subtype") +
    stat_regline_equation(
      mapping = aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
      formula = y ~ x, color = "red", size = 3,
      label.x = 32, label.y = 58, show.legend = FALSE) +
    stat_regline_equation(
      aes(color = Subtype,
          label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
      formula = y ~ x, size = 3,
      label.x = 40,
      label.y = seq(40, 32, length.out = n_subtypes),
      show.legend = FALSE) +
    theme_classic() +
    labs(title = segment, x = "GC3s", y = "GC12") +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title  = element_text(size = 9),
          axis.text   = element_text(size = 7),
          legend.position = "none")
  return(p)
}

segments_plot <- lapply(segments, process_gene)

combined_plot <- wrap_plots(segments_plot, ncol = 4) +
  theme(legend.position = "bottom",
        legend.justification = "right",
        legend.text = element_text(size = 11))

print(combined_plot)

regression_results <- data.frame()

for (segment in segments) {
  filename <- paste0("data/ENC_GC/Neutrality_", segment, ".csv")
  df <- read.csv(file = filename, header = TRUE)
  df$Subtype <- factor(df$Subtype, levels = c("H5N1", "H5N2", "H5N6", "H5N8", "Others"))
  
  segment_results <- df %>%
    group_by(Subtype) %>%
    do({
      if (nrow(.) >= 2) {
        model <- lm(GC12 ~ GC3s, data = .)
        data.frame(
          Segment   = segment,
          Slope     = coef(model)[2],
          Intercept = coef(model)[1],
          R_squared = summary(model)$r.squared,
          P_value   = summary(model)$coefficients[2, 4],
          Observations = nrow(.)
        )
      } else {
        data.frame(Segment = segment, Slope = NA, Intercept = NA,
                   R_squared = NA, P_value = NA, Observations = nrow(.))
      }
    }) %>%
    ungroup() %>%
    mutate(
      Mutation_pressure = round(Slope * 100, 1),
      Natural_selection = round((1 - Slope) * 100, 1),
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE ~ "NS"
      )
    )
  
  regression_results <- bind_rows(regression_results, segment_results)
}

print(regression_results)
write.csv(regression_results, file = "regression_results.csv", row.names = FALSE)

df_cg <- read.csv("data/ENC_GC/Neutrality_CG.csv", header = TRUE)
model_cg <- lm(GC12 ~ GC3s, data = df_cg)
cat("\nWhole-genome regression: slope =", round(coef(model_cg)[2], 3),
    ", R2 =", round(summary(model_cg)$r.squared, 4), "\n")

ggsave("Neutrality_CG.jpg", p_CG, width = 6, height = 4, dpi = 720, bg = "white")
ggsave("Neutrality_segments.jpg", combined_plot,
       width = 14, height = 8, device = "jpg", dpi = 300)