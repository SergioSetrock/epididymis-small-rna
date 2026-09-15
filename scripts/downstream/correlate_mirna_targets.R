# --- correlate_mirna_targets.R ---
suppressMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
small_counts_file <- args[1]
large_counts_file <- args[2]
meta_file         <- args[3]
target_dir        <- args[4]
out_dir           <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting miRNA-mRNA Integration and Correlation Analysis...")

# --- 1. Load and Clean Metadata ---
meta <- read.csv(meta_file, row.names=1)

# --- 2. Load and Normalize Small RNA Counts ---
message("Loading and normalizing Small RNA counts...")
small_counts <- read.table(small_counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
# Find common samples
common_small <- intersect(colnames(small_counts), rownames(meta))
small_counts <- small_counts[, common_small]

dds_small <- DESeqDataSetFromMatrix(countData = round(small_counts), 
                                    colData = meta[common_small, , drop=FALSE], 
                                    design = ~ 1)
dds_small <- estimateSizeFactors(dds_small)
# Log2 transform normalized counts (adding 1 to avoid log(0))
small_norm <- log2(counts(dds_small, normalized=TRUE) + 1)

# --- 3. Load and Normalize Large RNA Counts ---
message("Loading, mapping, and normalizing Large RNA counts...")
large_counts <- read.table(large_counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
common_large <- intersect(colnames(large_counts), rownames(meta))
large_counts <- large_counts[, common_large]

# Map Large RNA Ensembl IDs to Gene Symbols
clean_ids <- gsub("\\..*$", "", rownames(large_counts))
symbols <- mapIds(org.Mm.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace rownames with symbols (aggregate duplicates by summing counts)
large_counts$Symbol <- symbols
large_counts <- large_counts %>% 
  filter(!is.na(Symbol) & Symbol != "") %>%
  group_by(Symbol) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()
rownames(large_counts) <- large_counts$Symbol
large_counts$Symbol <- NULL

dds_large <- DESeqDataSetFromMatrix(countData = round(large_counts), 
                                    colData = meta[common_large, , drop=FALSE], 
                                    design = ~ 1)
dds_large <- estimateSizeFactors(dds_large)
large_norm <- log2(counts(dds_large, normalized=TRUE) + 1)

# Ensure we only compare the exact same samples between both datasets
final_samples <- intersect(colnames(small_norm), colnames(large_norm))
small_norm <- small_norm[, final_samples]
large_norm <- large_norm[, final_samples]

message(paste("Integration ready: Found", length(final_samples), "common samples."))

# --- 4. Correlate miRNA-Target Pairs ---
csv_files <- list.files(target_dir, pattern = "01_miRNA_Targets_.*\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No target tables found.")
}

for (file in csv_files) {
  comp_name <- gsub("01_miRNA_Targets_", "", basename(file))
  comp_name <- gsub("\\.csv$", "", comp_name)
  
  df_targets <- read.csv(file)
  if (nrow(df_targets) == 0) next
  
  message(paste("\nCalculating correlations for:", comp_name))
  
  results_list <- list()
  
  for (i in 1:nrow(df_targets)) {
    mature_id <- df_targets$mature_mirna_id[i]
    target_gene <- df_targets$target_symbol[i]
    
    # Check if target is in our Large RNA matrix
    if (!target_gene %in% rownames(large_norm)) next
    
    # Reverse-engineer the mature ID back to the Ensembl Symbol (e.g., mmu-miR-329-3p -> Mir329)
    core_mir <- gsub("mmu-mir-|mmu-miR-|-5p|-3p", "", mature_id)
    # Find the corresponding row in the small RNA matrix
    small_idx <- grep(paste0("(?i)^mir", core_mir, "$"), rownames(small_norm), perl=TRUE)
    
    if (length(small_idx) == 0) next
    
    mir_gene <- rownames(small_norm)[small_idx[1]] # Take the first match
    
    # Extract the expression vectors
    expr_mir <- as.numeric(small_norm[mir_gene, final_samples])
    expr_tar <- as.numeric(large_norm[target_gene, final_samples])
    
    # Calculate Spearman Correlation
    cor_res <- suppressWarnings(cor.test(expr_mir, expr_tar, method = "spearman"))
    
    results_list[[i]] <- data.frame(
      miRNA_Symbol = mir_gene,
      Mature_ID = mature_id,
      Target_Gene = target_gene,
      Database_Hits = df_targets$Hit_Count[i],
      Spearman_Rho = cor_res$estimate,
      P_Value = cor_res$p.value
    )
  }
  
  if (length(results_list) == 0) {
    message("  -> No matching pairs found in expression matrices.")
    next
  }
  
  df_cor <- bind_rows(results_list) %>%
    filter(!is.na(Spearman_Rho)) %>%
    # Filter for NEGATIVE correlations and standard significance
    filter(Spearman_Rho < 0 & P_Value < 0.05) %>%
    arrange(Spearman_Rho) # Most negative at the top
  
  if (nrow(df_cor) == 0) {
    message("  -> No significant negative correlations found.")
    next
  }
  
  # Save Table
  out_csv <- file.path(out_dir, paste0("01_Negative_Correlations_", comp_name, ".csv"))
  write.csv(df_cor, out_csv, row.names = FALSE)
  message(paste("  -> Found", nrow(df_cor), "significant negative correlations. Saved table."))
  
  # --- 5. Plot Top 6 Correlated Pairs ---
  top_pairs <- head(df_cor, 6)
  
  for (j in 1:nrow(top_pairs)) {
    mir <- top_pairs$miRNA_Symbol[j]
    tar <- top_pairs$Target_Gene[j]
    rho <- round(top_pairs$Spearman_Rho[j], 3)
    pval <- signif(top_pairs$P_Value[j], 3)
    
    plot_data <- data.frame(
      Sample = final_samples,
      miRNA = as.numeric(small_norm[mir, final_samples]),
      mRNA = as.numeric(large_norm[tar, final_samples])
    )
    
    p <- ggplot(plot_data, aes(x = miRNA, y = mRNA)) +
      geom_point(size = 3, color = "darkblue", alpha = 0.7) +
      geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
      theme_bw() +
      labs(title = paste(mir, "vs", tar),
           subtitle = paste("Spearman Rho:", rho, "| P-value:", pval),
           x = paste(mir, "Expression (log2 Norm Counts)"),
           y = paste(tar, "Expression (log2 Norm Counts)")) +
      theme(plot.title = element_text(size = 14, face = "bold", hjust=0.5),
            plot.subtitle = element_text(size = 11, hjust=0.5))
    
    ggsave(file.path(out_dir, paste0("02_Scatter_", comp_name, "_", mir, "_vs_", tar, ".png")), plot = p, width = 6, height = 5)
  }
  message("  -> Generated scatter plots for top pairs.")
}

message("\nIntegration complete!")
