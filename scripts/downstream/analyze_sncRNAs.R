# --- analyze_sncRNAs.R ---
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
deg_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting snRNA and snoRNA functional classification...")

csv_files <- list.files(deg_dir, pattern = "01_Significant_DEGs_.*_annotated\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No annotated significant DEG tables found.")
}

for (file in csv_files) {
  comp_name <- gsub("01_Significant_DEGs_", "", basename(file))
  comp_name <- gsub("_annotated\\.csv$", "", comp_name)
  
  df <- read.csv(file)
  
  # 1. Filter strictly for Differentially Expressed snRNAs and snoRNAs
  df_snc <- df %>%
    filter(Regulation != "Not_Significant" & Biotype %in% c("snRNA", "snoRNA"))
  
  if (nrow(df_snc) == 0) {
    message(paste("  -> No significant snRNAs/snoRNAs found in", comp_name, "- Skipping."))
    next
  }
  
  message(paste("\nProcessing", nrow(df_snc), "snRNAs/snoRNAs for:", comp_name))
  
  # 2. Classify Subtypes based on Gene Symbols
  # snoRNAs: SNORD (C/D box - Methylation), SNORA (H/ACA box - Pseudouridylation)
  # snRNAs: usually Rnu or U-RNAs (Spliceosome)
  df_snc <- df_snc %>%
    mutate(
      Subtype = case_when(
        str_detect(Symbol, "(?i)snord") ~ "snoRNA: C/D Box (rRNA Methylation)",
        str_detect(Symbol, "(?i)snora") ~ "snoRNA: H/ACA Box (rRNA Pseudouridylation)",
        Biotype == "snoRNA" ~ "snoRNA: Other/Unclassified",
        Biotype == "snRNA" ~ "snRNA: Spliceosome Component",
        TRUE ~ "Unknown"
      )
    ) %>%
    arrange(desc(log2FoldChange))
  
  # Lock the factor levels so the waterfall plot orders them perfectly by fold change
  df_snc$Symbol <- factor(df_snc$Symbol, levels = df_snc$Symbol)
  
  # Save the classified table
  out_csv <- file.path(out_dir, paste0("01_Classified_sncRNAs_", comp_name, ".csv"))
  write.csv(df_snc, out_csv, row.names = FALSE)
  
  # 3. Create a Waterfall Plot of Log2FoldChanges
  p_waterfall <- ggplot(df_snc, aes(x = Symbol, y = log2FoldChange, fill = Subtype)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    coord_flip() + # Flip to make gene names readable
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    labs(title = paste("snRNA & snoRNA Regulation\n", gsub("_", " ", comp_name)),
         subtitle = "Impact on Splicing and Ribosomal Machinery",
         x = "Gene Symbol",
         y = "Log2 Fold Change",
         fill = "Functional Subtype") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.text.y = element_text(size = 9, face = "bold"),
          legend.position = "bottom",
          legend.direction = "vertical") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  
  # Dynamic height based on number of genes
  plot_height <- max(5, nrow(df_snc) * 0.25)
  
  ggsave(file.path(out_dir, paste0("02_Waterfall_sncRNAs_", comp_name, ".png")), plot = p_waterfall, width = 8, height = plot_height)
  
  message("  -> Success! Generated classification table and waterfall plot.")
}

message("\nsnRNA and snoRNA analysis completed!")
