# --- enrich_mirna_targets.R ---
suppressMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
target_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting GO pathway enrichment for miRNA targets...")

# Find the target files we just generated
csv_files <- list.files(target_dir, pattern = "01_miRNA_Targets_.*\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No target tables found in the specified directory.")
}

for (file in csv_files) {
  comp_name <- gsub("01_miRNA_Targets_", "", basename(file))
  comp_name <- gsub("\\.csv$", "", comp_name)
  
  df <- read.csv(file)
  
  if (nrow(df) == 0) {
    message(paste("  -> The table for", comp_name, "is empty. Skipping."))
    next
  }
  
  message(paste("\nProcessing enrichment for:", comp_name))
  
  # Extract unique target genes
  # Optional: We could filter here by Hit_Count > 1 to use only highly robust targets,
  # but we will use all unique ones for a global view.
  target_genes <- unique(df$target_symbol)
  message(paste("  -> Evaluating", length(target_genes), "unique target genes..."))
  
  # Run Gene Ontology (Biological Process)
  go_enrich <- enrichGO(gene          = target_genes,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "SYMBOL", # Our tables already use gene symbols
                        ont           = "BP",     # Biological Process
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.1,
                        readable      = FALSE)
  
  if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
    message("  -> No significantly enriched GO pathways found.")
    writeLines("No significant GO pathways found.", file.path(out_dir, paste0("02_GO_Enrichment_", comp_name, "_NO_SIG.txt")))
    next
  }
  
  # Save the results table
  res_table <- as.data.frame(go_enrich)
  write.csv(res_table, file.path(out_dir, paste0("02_GO_Enrichment_", comp_name, ".csv")), row.names = FALSE)
  
  # Generate a Dotplot with the top 15 results
  num_paths <- min(15, nrow(res_table))
  p_dot <- dotplot(go_enrich, showCategory = num_paths) +
    theme_bw() +
    labs(title = paste("Biological Pathways Targeted by DE miRNAs\n", gsub("_", " ", comp_name))) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10))
  
  ggsave(file.path(out_dir, paste0("03_GO_Dotplot_", comp_name, ".png")), plot = p_dot, width = 9, height = 7)
  
  message("  -> Success! Pathway table and Dotplot generated.")
}

message("\nPathway analysis completed!")
