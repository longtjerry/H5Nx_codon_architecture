suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
  library(tidyr)
  library(ggtext)
})

df <- read.csv(file = "data/CAI/metadata_CAI.csv", header = TRUE,
               na.strings = c("", " ", "N/A", "NULL")) %>%
  filter(Gene != "CG") %>%
  filter(Organ != "Genomic")

subtype_color <- c("H5N1" = "#953638", "H5N2" = "#1c516c", "H5N6" = "#d8cd9d",
                   "H5N8" = "#fda649", "Others" = "#C5C6C7")

organ_colors <- c(
  "Heart"           = "#E74C3C",
  "Kidney"          = "#3498DB",
  "Brain"           = "#9B59B6",
  "Lung"            = "#1ABC9C",
  "Small intestine" = "#F39C12",
  "Liver"           = "#D35400",
  "Muscle"          = "#7F8C8D"
)

section_abbr <- c(
  "Heart-Atrial Appendage"          = "Hrt-AA",
  "Heart-Left Ventricle"             = "Hrt-LV",
  "Kidney-Cortex"                   = "Kid-Ctx",
  "Kidney-Medulla"                  = "Kid-Med",
  "Brain-Amygdala"                  = "Amy",
  "Brain-Anterior cingulate cortex" = "ACC",
  "Brain-Caudate"                   = "Cau",
  "Brain-Cerebellar Hemisphere"     = "Cer-H",
  "Brain-Cerebellum"                = "Cer",
  "Brain-Cortex"                    = "Ctx",
  "Brain-Frontal Cortex"            = "F-Ctx",
  "Brain-Hippocampus"               = "Hip",
  "Brain-Hypothalamus"              = "Hyp",
  "Brain-Nucleus accumbens"         = "NAc",
  "Brain-Spinal cord"               = "S-Cord",
  "Brain-Putamen"                   = "Put",
  "Brain-Substantia nigra"          = "SN",
  "Brain-Pituitary"                 = "Pit",
  "Lung"                            = "Lung",
  "Small intestine"                 = "S-Int",
  "Liver"                           = "Liv",
  "Muscle-Skeletal"                 = "Sk-M"
)


section_order <- c("Heart-Atrial Appendage", "Heart-Left Ventricle", "Kidney-Cortex",
                   "Kidney-Medulla", "Brain-Amygdala", "Brain-Anterior cingulate cortex",
                   "Brain-Caudate", "Brain-Cerebellar Hemisphere", "Brain-Cerebellum",
                   "Brain-Cortex", "Brain-Frontal Cortex", "Brain-Hippocampus",
                   "Brain-Hypothalamus", "Brain-Nucleus accumbens", "Brain-Spinal cord",
                   "Brain-Putamen", "Brain-Substantia nigra", "Brain-Pituitary",
                   "Lung", "Small intestine", "Liver", "Muscle-Skeletal")

gene_summary_data <- df %>%
  group_by(Subtype, Gene, Section) %>%
  summarise(
    mean_CAI = mean(CAI, na.rm = TRUE),
    sd_CAI   = sd(CAI, na.rm = TRUE),
    n        = n(),
    se_CAI   = sd_CAI / sqrt(n),
    ci_lower = mean_CAI - qt(0.975, n - 1) * se_CAI,
    ci_upper = mean_CAI + qt(0.975, n - 1) * se_CAI,
    .groups  = "drop"
  ) %>%
  mutate(Section = factor(Section, levels = section_order))

section_organ_map <- df %>%
  dplyr::select(Section, Organ) %>%
  distinct()

colored_labels <- sapply(section_order, function(section) {
  organ <- section_organ_map$Organ[section_organ_map$Section == section][1]
  color <- organ_colors[organ]
  if (is.na(color)) color <- "#000000"
  paste0("<span style='color:", color, "'>", section_abbr[section], "</span>")
})

### ========== Panel A: tissue-specific CAI =======

p1 <- ggplot(gene_summary_data, aes(x = Section, y = mean_CAI,
                                    color = Subtype, fill = Subtype,
                                    group = Subtype)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.15, color = NA) +
  geom_line(size = 0.8) +
  geom_point(size = 1.2, shape = 21, fill = "white", stroke = 0.8) +
  facet_wrap(~ Gene, ncol = 4) +
  theme_bw() +
  theme(axis.text.x  = ggtext::element_markdown(angle = 45, hjust = 1, size = 9),
        axis.text.y  = element_text(size = 9),
        axis.title   = element_text(face = "bold", size = 12),
        strip.text   = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold"),
        legend.position  = "bottom",
        legend.text  = element_text(size = 14),
        panel.grid.minor = element_blank(),
        plot.margin  = unit(c(1, 1, 1, 1), "cm")) +
  labs(y = "Codon Adaptation Index (CAI)",
       x = "Organ/Tissue Section",
       color = "H5Nx Subtype:",
       fill  = "H5Nx Subtype:") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = colored_labels) +
  scale_color_manual(values = subtype_color) +
  scale_fill_manual(values = subtype_color)

print(p1)

### ==================== Panel B: human genome reference CAI ====================

genomic_df <- read.csv(file = "data/CAI/metadata_genomic.csv", header = TRUE,
                       na.strings = c("", " ", "N/A", "NULL"))

genomic_summary <- genomic_df %>%
  group_by(Subtype, Gene) %>%
  summarise(
    mean_CAI = mean(CAI, na.rm = TRUE),
    sd_CAI   = sd(CAI, na.rm = TRUE),
    n        = n(),
    se_CAI   = sd_CAI / sqrt(n),
    ci_lower = mean_CAI - qt(0.975, n - 1) * se_CAI,
    ci_upper = mean_CAI + qt(0.975, n - 1) * se_CAI,
    .groups  = "drop"
  )

p_genomic <- ggplot(genomic_summary, aes(x = Gene, y = mean_CAI,
                                         color = Subtype, fill = Subtype,
                                         group = Subtype)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.2, color = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title  = element_text(face = "bold", size = 10),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title  = element_text(hjust = .5, face = "bold", size = 14),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  labs(title = "Human genome",
       y = "Codon Adaptation Index (CAI)",
       x = "Viral Gene",
       color = "Virus Subtype",
       fill  = "Virus Subtype") +
  scale_color_manual(values = subtype_color) +
  scale_fill_manual(values = subtype_color)

print(p_genomic)

genomic_ranking <- genomic_summary %>%
  group_by(Subtype) %>%
  summarise(Overall_Genomic_CAI = mean(mean_CAI), .groups = "drop") %>%
  arrange(desc(Overall_Genomic_CAI))
cat("\nOverall genomic CAI by subtype:\n")
print(genomic_ranking)

gene_ranking <- genomic_summary %>%
  group_by(Gene) %>%
  summarise(Mean_Genomic_CAI = mean(mean_CAI), .groups = "drop") %>%
  arrange(desc(Mean_Genomic_CAI))
cat("\nGenomic CAI by gene:\n")
print(gene_ranking)

# Export rankings for Table S2
write.csv(genomic_ranking, "genomic_CAI_ranking_by_subtype.csv", row.names = FALSE)
write.csv(gene_ranking, "genomic_CAI_ranking_by_gene.csv", row.names = FALSE)


ggsave(filename = "CAI_tissue.jpg",
       plot = p1,
       device = "jpeg",
       width = 14, height = 8, units = "in",
       dpi = 300, bg = "white")

ggsave(filename = "CAI_genomic.jpg",
       plot = p_genomic,
       device = "jpeg",
       width = 6, height = 4, units = "in",
       dpi = 300, bg = "white")