suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
})


section_levels <- c("Heart-Atrial Appendage", "Heart-Left Ventricle",
                    "Kidney-Cortex", "Kidney-Medulla",
                    "Brain-Amygdala", "Brain-Anterior cingulate cortex",
                    "Brain-Caudate", "Brain-Cerebellar Hemisphere",
                    "Brain-Cerebellum", "Brain-Cortex",
                    "Brain-Frontal Cortex", "Brain-Hippocampus",
                    "Brain-Hypothalamus", "Brain-Nucleus accumbens",
                    "Brain-Spinal cord", "Brain-Putamen",
                    "Brain-Substantia nigra", "Brain-Pituitary",
                    "Lung", "Small intestine", "Liver", "Muscle-Skeletal")

### ============== Panel A: strain highlight ============

df <- read.csv(file = "data/CAI/metadata_CAI.csv", header = TRUE,
               na.strings = c("", " ", "N/A", "NULL"))

df$Section <- factor(df$Section, levels = section_levels)
facet_order <- sort(unique(df$Gene))
df$Gene <- factor(df$Gene, levels = facet_order)

name_map <- c("A/Yangzhou/125/2022_H5N6|GCA_039266715.1" = "YZ125_H5N6",
              "A/Hong_Kong/483/1997"  = "HK483_H5N1",
              "A/Hong_Kong/486/97"    = "HK486_H5N1",
              "A/Indonesia/5/2005"    = "IDN5_H5N1")
highlight_strain <- names(name_map)

data <- df %>%
  filter(Section != "Genomic") %>%
  mutate(highlight = ifelse(Strain %in% highlight_strain,
                            name_map[Strain],
                            "Other H5Nx"))

data_highlight <- data %>% filter(highlight != "Other H5Nx")
stopifnot(nrow(data_highlight) > 0)

p1 <- ggplot(data, aes(x = Section, y = CAI)) +
  coord_flip() +
  geom_violin(fill = "grey90", color = "grey60") +
  geom_point(data = data_highlight,
             aes(color = highlight), size = 1, alpha = 0.9) +
  scale_color_manual(
    values = c("YZ125_H5N6" = "#E41A1C",
               "HK483_H5N1" = "#377EB8",
               "HK486_H5N1" = "#4DAF4A",
               "IDN5_H5N1"  = "#984EA3"),
    name = "Strain") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 16))) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8)) +
  facet_wrap(~ Gene, ncol = 3)

print(p1)

ggsave(filename = "CAI_highlight.jpg",
       plot = p1,
       device = "jpeg",
       width = 8, height = 8, units = "in",
       dpi = 300, bg = "white")

### ============== Panel B: host-grouped CAI + statistics ====================

suppressPackageStartupMessages({
  library(tidyr)
  library(ggtext)
  library(TOSTER)
})

df <- read.csv(file = "data/CAI/metadata_CAI.csv", header = TRUE,
               na.strings = c("", " ", "N/A", "NULL")) %>%
  filter(Gene == "CG") %>%        
  filter(Organ != "Genomic")      

kw_p <- df %>%
  group_by(Section) %>%
  summarise(p = kruskal.test(CAI ~ Host)$p.value, .groups = "drop")

effsize_df <- df %>%
  group_by(Section) %>%
  kruskal_effsize(CAI ~ Host) %>% 
  dplyr::select(Section, n, effsize, magnitude) %>%
  left_join(kw_p, by = "Section") %>%
  mutate(eta2H_label = paste0("eta2H = ", sprintf("%.4f", effsize)))

write.csv(effsize_df, "KW_effect_size_by_section.csv", row.names = FALSE)
print(effsize_df)

n_negligible <- sum(effsize_df$effsize < 0.01, na.rm = TRUE)
message(sprintf("%d/%d sections with eta2H < 0.01", n_negligible, nrow(effsize_df)))

dunn_df <- df %>%
  group_by(Section) %>%
  dunn_test(CAI ~ Host, p.adjust.method = "BH") %>%
  dplyr::select(Section, group1, group2, n1, n2, p, p.adj, p.adj.signif)

write.csv(dunn_df, "Dunn_posthoc_by_section.csv", row.names = FALSE)

sig_pairs <- dunn_df %>%
  filter(p.adj < 0.05) %>%
  mutate(pair = paste(group1, "vs", group2))
print(table(sig_pairs$pair))

bound <- 0.01

run_tost <- function(dat, g1, g2) {
  x <- dat$CAI[dat$Host == g1]
  y <- dat$CAI[dat$Host == g2]
  if (length(x) < 2 || length(y) < 2) {
    return(data.frame(comparison = paste(g1, "vs", g2),
                      mean_diff = NA, ci_low = NA, ci_high = NA,
                      p_tost = NA, equivalent = NA))
  }
  tt  <- t.test(y, x)
  res <- t_TOST(x = y, y = x,
                low_eqbound = -bound, high_eqbound = bound,
                eqbound_type = "raw", alpha = 0.05)
  
  pv <- res$TOST$p.value
  tost_rows <- grep("TOST", rownames(res$TOST)) 
  p2 <- if (length(tost_rows) == 2) max(pv[tost_rows]) else max(pv)
  
  data.frame(
    comparison = paste(g1, "vs", g2),
    mean_diff  = mean(y) - mean(x),
    ci_low     = tt$conf.int[1],
    ci_high    = tt$conf.int[2],
    p_tost     = p2,
    equivalent = p2 < 0.05
  )
}

tost_results <- df %>%
  group_by(Section) %>%
  group_modify(~ bind_rows(
    run_tost(.x, "Avian", "Human"),
    run_tost(.x, "Avian", "Nonhuman Mammal"),
    run_tost(.x, "Human", "Nonhuman Mammal")
  )) %>%
  ungroup()

write.csv(tost_results, "TOST_equivalence_by_section.csv", row.names = FALSE)

tost_summary <- tost_results %>%
  group_by(comparison) %>%
  summarise(n_equivalent = sum(equivalent), total = n(), .groups = "drop")
print(tost_summary)

n_table <- df %>%
  count(Host, Section) %>%
  pivot_wider(names_from = Host, values_from = n, values_fill = 0)
write.csv(n_table, "sample_size_by_section.csv", row.names = FALSE)

section_organ <- data.frame(
  Section = section_levels,
  Organ = c(rep("Heart", 2), rep("Kidney", 2), rep("Brain", 14),
            "Lung", "Small intestine", "Liver", "Muscle")
)

summary_data <- df %>%
  group_by(Host, Section) %>%
  summarise(
    mean_CAI = mean(CAI, na.rm = TRUE),
    sd_CAI   = sd(CAI, na.rm = TRUE),
    n        = n(),
    se_CAI   = sd_CAI / sqrt(n),
    ci_lower = mean_CAI - qt(0.975, n - 1) * se_CAI,
    ci_upper = mean_CAI + qt(0.975, n - 1) * se_CAI,
    .groups  = "drop"
  ) %>%
  filter(n >= 2)

summary_data$Section <- factor(summary_data$Section, levels = section_levels)

host_colors  <- c("Human" = "#E41A1C", "Avian" = "#377EB8",
                  "Nonhuman Mammal" = "#4DAF4A")
organ_colors <- c("Heart" = "#E74C3C", "Kidney" = "#3498DB", "Brain" = "#9B59B6",
                  "Lung" = "#1ABC9C", "Small intestine" = "#F39C12",
                  "Liver" = "#D35400", "Muscle" = "#7F8C8D")

colored_labels <- sapply(section_levels, function(s) {
  organ <- section_organ$Organ[section_organ$Section == s]
  paste0("<span style='color:", organ_colors[organ], "'>", s, "</span>")
})

label_df <- effsize_df %>%
  left_join(summary_data %>%
              group_by(Section) %>%
              summarise(y_max = max(ci_upper, na.rm = TRUE), .groups = "drop"),
            by = "Section") %>%
  mutate(y_pos   = y_max * 1.005,
         eta2H_lab = sprintf("%.3f", effsize))

p_final <- ggplot(summary_data, aes(x = Section, y = mean_CAI,
                                    color = Host, fill = Host, group = Host)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_text(data = label_df,
            aes(x = Section, y = y_pos, label = eta2H_lab),
            inherit.aes = FALSE, size = 2.6, vjust = 0) +
  annotate("text", x = 1, y = max(summary_data$ci_upper, na.rm = TRUE) * 1.035,
           label = "eta2H (Kruskal-Wallis effect size)", hjust = 0, size = 3) +
  scale_x_discrete(labels = colored_labels) +
  scale_color_manual(values = host_colors) +
  scale_fill_manual(values = host_colors) +
  labs(y = "Codon Adaptation Index (CAI)", x = "Organ/Tissue Section") +
  theme_bw() +
  theme(axis.text.x  = ggtext::element_markdown(angle = 45, hjust = 1, size = 11),
        axis.text.y  = element_text(size = 9),
        axis.title   = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.margin  = unit(c(1, 1, 1, 1), "cm")) +
  coord_cartesian(ylim = c(min(summary_data$ci_lower, na.rm = TRUE),
                           max(summary_data$ci_upper, na.rm = TRUE) * 1.05))

print(p_final)

ggsave("CAI_host.jpg", p_final, width = 12, height = 7, dpi = 300)