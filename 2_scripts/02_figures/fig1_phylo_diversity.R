suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ape)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(phangorn)
})

### ============== Panel A: phylogenetic trees ====================

tree_HA <- treeio::read.tree("data/trees/gene_HA_aln.fasta.treefile")
tree_NA <- treeio::read.tree("data/trees/gene_NA_aln.fasta.treefile")
data    <- read.csv("data/subtype.csv", header = TRUE,
                    na.strings = c("", " ", "N/A", "NULL"))

headers <- readLines("data/trees/list.txt")
headers <- sub("^>", "", headers)
stopifnot(length(headers) > 0)

# UFBoot support summary (informational)
bs <- as.numeric(tree_HA$node.label)
summary(bs)
table(cut(bs, c(0, 70, 80, 95, 100)))

bs <- as.numeric(tree_NA$node.label)
summary(bs)
table(cut(bs, c(0, 70, 80, 95, 100)))

make_key <- function(x) {
  x <- sub("\\|.*$", "", x)                     
  x <- gsub("[()]", "_", x)                       
  x <- sub("[ _]*H[0-9]N[0-9]_*$", "", x)        
  x <- gsub(" ", "_", x)                         
  x <- gsub("_+", "_", x)                        
  sub("_$", "", x)                               
}

data$Key    <- make_key(data$Seq_name)
headers_key <- make_key(headers)
tree_HA$Key <- make_key(tree_HA$tip.label)
tree_NA$Key <- make_key(tree_NA$tip.label)

a_ha <- sum(!tree_HA$Key %in% headers_key)
a_na <- sum(!tree_NA$Key %in% headers_key)
cat("(a) Tree tips not present in fasta - HA:", a_ha, "| NA:", a_na, "\n")
stopifnot(a_ha == 0, a_na == 0)

miss_csv <- setdiff(headers_key, data$Key)
cat("(b) Fasta entries missing from CSV:", length(miss_csv), "\n")
if (length(miss_csv) > 0) print(head(miss_csv, 10))

dup_conflict <- data %>%
  group_by(Key) %>%
  summarise(n_subtype = n_distinct(Subtype), .groups = "drop") %>%
  filter(n_subtype > 1)
cat("(c) Keys with conflicting Subtype assignments:", nrow(dup_conflict), "\n")
if (nrow(dup_conflict) > 0) {
  print(data %>% filter(Key %in% dup_conflict$Key))
  stop("Conflicting subtype assignments found - resolve manually before continuing")
}

key2subtype <- data$Subtype
names(key2subtype) <- data$Key

make_plot_data <- function(tree) {
  st <- key2subtype[tree$Key]
  data.frame(Seq_name = tree$tip.label,
             Subtype  = unname(st),
             stringsAsFactors = FALSE)
}
plot_data_HA <- make_plot_data(tree_HA)
plot_data_NA <- make_plot_data(tree_NA)
n_unmatched  <- sum(is.na(plot_data_HA$Subtype)) + sum(is.na(plot_data_NA$Subtype))
cat("Total unmatched tips:", n_unmatched, "\n")
stopifnot(n_unmatched == 0)

subtype_color <- c("H5N1" = "#953638", "H5N2" = "#1c516c", "H5N6" = "#d8cd9d",
                   "H5N8" = "#fda649", "Others" = "#C5C6C7")

plot_tree <- function(tree, plot_data, caption) {
  n_tip <- length(tree$tip.label)
  ggtree(tree) %<+% plot_data +
    geom_tree(size = 1) +
    geom_tippoint(aes(fill = Subtype), stroke = .5, shape = 21, size = 2) +
    scale_fill_manual(name = "Subtype", values = subtype_color,
                      guide = guide_legend(keywidth = .3, keyheight = .3, ncol = 2,
                                           override.aes = list(size = 6, alpha = 1))) +
    geom_treescale(fontsize = 3) +
    coord_cartesian(clip = "off",
                    ylim = c(1 - n_tip * 0.01, n_tip + n_tip * 0.01)) +
    theme_tree() +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.margin  = margin(10, 20, 10, 10))
}

p1 <- plot_tree(tree_HA, plot_data_HA, "HA")
p2 <- plot_tree(tree_NA, plot_data_NA, "NA")

calculate_multi_state_association <- function(tree, trait, n_permutations = 1000) {
  trait_factor <- as.factor(trait[tree$tip.label])
  stopifnot(!any(is.na(trait_factor)))
  
  phy <- phangorn::as.phyDat(trait_factor)
  observed_score <- phangorn::parsimony(tree, phy)
  
  perm_once <- function(i) {
    rt <- sample(trait_factor)
    names(rt) <- names(trait_factor)
    phangorn::parsimony(tree, phangorn::as.phyDat(rt))
  }
  
  if (.Platform$OS.type == "windows") {
    n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("trait_factor", "tree"), envir = environment())
    parallel::clusterEvalQ(cl, library(phangorn))
    random_scores <- parallel::parLapply(cl, 1:n_permutations, perm_once) |> unlist()
    parallel::stopCluster(cl)
  } else {
    random_scores <- parallel::mclapply(1:n_permutations, perm_once,
                                        mc.cores = max(1, parallel::detectCores() - 1)) |> unlist()
  }
  
  p_value     <- (sum(random_scores <= observed_score) + 1) / (n_permutations + 1)
  effect_size <- (mean(random_scores) - observed_score) / sd(random_scores)
  
  list(observed_parsimony    = observed_score,
       mean_random_parsimony = mean(random_scores),
       sd_random_parsimony   = sd(random_scores),
       p_value               = p_value,
       effect_size           = effect_size,
       n_states              = length(levels(trait_factor)),
       n_permutations        = n_permutations)
}

trait_vector <- setNames(key2subtype[tree_HA$Key], tree_HA$tip.label)
result <- calculate_multi_state_association(tree_HA, trait_vector,
                                            n_permutations = 1000)

cat("\n===== Multi-state Association Test Results =====\n")
cat("Observed parsimony:", result$observed_parsimony, "\n")
cat("Random mean +/- SD:", round(result$mean_random_parsimony, 1), "+/-",
    round(result$sd_random_parsimony, 1), "\n")
cat("p =", result$p_value, "\n")
cat("Effect size (z):", round(result$effect_size, 3), "\n")

sig <- result$p_value < 0.05
cat(ifelse(sig,
           "CONCLUSION: significant phylogenetic association - subtypes are non-randomly distributed on the HA tree\n",
           "CONCLUSION: no significant phylogenetic association - subtypes appear randomly distributed on the HA tree\n"))

write.csv(data.frame(t(unlist(result))), "HA_tree_subtype_association_test.csv",
          row.names = FALSE)


significance_message <- ifelse(
  result$p_value < 0.001, "Significant phylogenetic association (p < 0.001)",
  sprintf("Significant phylogenetic association (p = %.4f)", result$p_value)
)

combined_plot <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    subtitle = significance_message,
    theme = theme(plot.subtitle = element_text(size = 14, hjust = 0.5,
                                               color = ifelse(sig, "red", "black")))
  )

print(combined_plot)

ggsave(
  filename = "combined_tree.png",
  plot     = combined_plot,
  device   = "png",
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

### ============ Panel B: nucleotide diversity (pi) ====================

suppressPackageStartupMessages({
  library(pegas)
  library(Biostrings)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})

files <- list.files(pattern = "gene_.*_aln\\.fasta$")
cat(length(files), "alignment files found\n\n")

length_report <- data.frame(
  file = files,
  n_seq = NA_integer_,
  min_len = NA_integer_,
  max_len = NA_integer_,
  uniform = NA
)

for (i in seq_along(files)) {
  f <- files[i]
  seqs <- try(readDNAStringSet(f), silent = TRUE)
  if (inherits(seqs, "try-error")) {
    cat("[read failed]", f, "\n")
    next
  }
  w <- width(seqs)
  length_report[i, ] <- data.frame(
    file = f,
    n_seq = length(seqs),
    min_len = min(w),
    max_len = max(w),
    uniform = var(w) == 0
  )
  if (var(w) != 0) {
    cat("[inconsistent length]", f, ": n =", length(seqs),
        ", min =", min(w), ", max =", max(w),
        " -> re-check in MAFFT or remove abnormal sequences\n")
  } else {
    cat("[OK]", f, ": n =", length(seqs), ", length =", max(w), "\n")
  }
}

write.csv(length_report, "alignment_length_check.csv", row.names = FALSE)
cat("\nFiles with inconsistent lengths:", sum(!length_report$uniform, na.rm = TRUE), "\n\n")


safe_nuc_div <- function(fasta_file) {
  if (!file.exists(fasta_file)) return(NA_real_)
  seqs <- try(readDNAStringSet(fasta_file), silent = TRUE)
  if (inherits(seqs, "try-error")) return(NA_real_)
  
  w <- width(seqs)
  if (var(w) != 0) {
    warning(basename(fasta_file), " inconsistent sequence lengths (min=", min(w),
            ", max=", max(w), "); returning NA. Re-check in MAFFT first.")
    return(NA_real_)
  }
  nuc.div(as.DNAbin(seqs), pairwise.deletion = TRUE)
}

get_n <- function(fasta_file) {
  if (!file.exists(fasta_file)) return(NA_integer_)
  length(try(readDNAStringSet(fasta_file), silent = TRUE))
}


segments <- c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2")
subtypes <- c("H5N1", "H5N2", "H5N6", "H5N8", "Others")

pi_matrix <- matrix(NA_real_, nrow = length(segments), ncol = length(subtypes),
                    dimnames = list(segments, subtypes))
n_matrix <- pi_matrix

for (seg in segments) {
  for (sub in subtypes) {
    f <- paste0("gene_", seg, "_", sub, "_aln.fasta")
    pi_matrix[seg, sub] <- safe_nuc_div(f)
    n_matrix[seg, sub]  <- get_n(f)
    cat(sprintf("%-4s %-6s pi = %.5f (n = %s)\n", seg, sub,
                ifelse(is.na(pi_matrix[seg, sub]), NA, pi_matrix[seg, sub]),
                ifelse(is.na(n_matrix[seg, sub]), "NA", n_matrix[seg, sub])))
  }
}


write.csv(pi_matrix, "pi_diversity_matrix.csv")
write.csv(n_matrix,  "pi_diversity_n_matrix.csv")

pi_df <- as.data.frame(pi_matrix) %>%
  rownames_to_column("Segment") %>%
  pivot_longer(-Segment, names_to = "Subtype", values_to = "Pi") %>%
  left_join(
    as.data.frame(n_matrix) %>%
      rownames_to_column("Segment") %>%
      pivot_longer(-Segment, names_to = "Subtype", values_to = "n"),
    by = c("Segment", "Subtype")
  )
write.csv(pi_df, "pi_diversity_long.csv", row.names = FALSE)

cell_labels <- matrix(
  ifelse(is.na(pi_matrix), "n.a.",
         paste0(sprintf("%.4f", pi_matrix), "\n(n=", n_matrix, ")")),
  nrow = nrow(pi_matrix), dimnames = dimnames(pi_matrix)
)

pheatmap(pi_matrix[, c("H5N1","H5N2","H5N6","H5N8")],
         main = "Nucleotide diversity (pi) of major H5Nx subtypes",
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "#c0392b"))(100),
         display_numbers = cell_labels[, c("H5N1","H5N2","H5N6","H5N8")],
         fontsize = 11,
         na_col = "grey90",
         filename = "pi_diversity_heatmap.jpg",
         width = 8, height = 8)

cat("\nDone. Outputs: combined_tree.png, pi_diversity_heatmap.jpg,\n")
cat("HA_tree_subtype_association_test.csv, alignment_length_check.csv, pi_diversity_*.csv\n")