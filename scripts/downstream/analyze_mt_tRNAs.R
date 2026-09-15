# --- analyze_mt_tRNAs.R ---
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
deg_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting mitochondrial tRNA (mt_tRNA) analysis...")

csv_files <- list.files(deg_dir, pattern = "01_Significant_DEGs_.*_annotated\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No annotated significant DEG tables found.")
}

for (file in csv_files) {
  comp_name <- gsub("01_Significant_DEGs_", "", basename(file))
  comp_name <- gsub("_annotated\\.csv$", "", comp_name)
  
  df <- read.csv(file)
  
  # 1. Filter strictly for mt_tRNAs (even if they are not significant, to see the trend, 
  # but let's stick to your significant tables if that's what was passed)
  df_mt <- df %>%
    filter(Biotype == "mt_tRNA")
  
  if (nrow(df_mt) == 0) {
    message(paste("  -> No mt_tRNAs found in", comp_name, "- Skipping."))
    next
  }
  
  message(paste("\nProcessing", nrow(df_mt), "mt_tRNAs for:", comp_name))
  
  # 2. Clean up the data for plotting
  df_mt <- df_mt %>%
    # Ensure we drop any NA symbols
    filter(!is.na(Symbol) & Symbol != "") %>%
    arrange(log2FoldChange) %>%
    mutate(
      Direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
      # Extract the specific amino acid letter from the Ensembl symbol if formatted like 'mt-Ta'
      Amino_Acid = str_extract(Symbol, "(?<=mt-T)[a-zA-Z]") 
    )
  
  # Lock factor levels for ordered plotting
  df_mt$Symbol <- factor(df_mt$Symbol, levels = df_mt$Symbol)
  
  # Save the table
  out_csv <- file.path(out_dir, paste0("01_mt_tRNAs_", comp_name, ".csv"))
  write.csv(df_mt, out_csv, row.names = FALSE)
  
  # 3. Create a Lollipop Chart
  p_lolli <- ggplot(df_mt, aes(x = Symbol, y = log2FoldChange, color = Direction)) +
    geom_segment(aes(x = Symbol, xend = Symbol, y = 0, yend = log2FoldChange), color = "gray50", size = 1) +
    geom_point(size = 4) +
    coord_flip() +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue")) +
    theme_bw() +
    labs(title = paste("Mitochondrial tRNA Regulation\n", gsub("_", " ", comp_name)),
         subtitle = "Impact on Mitochondrial Translation Machinery",
         x = "Mitochondrial tRNA",
         y = "Log2 Fold Change") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.position = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  
  plot_height <- max(4, nrow(df_mt) * 0.3)
  
  ggsave(file.path(out_dir, paste0("02_Lollipop_mt_tRNAs_", comp_name, ".png")), plot = p_lolli, width = 7, height = plot_height)
  
  message("  -> Success! Generated mt_tRNA table and Lollipop plot.")
}

message("\nmt_tRNA analysis completed!")
