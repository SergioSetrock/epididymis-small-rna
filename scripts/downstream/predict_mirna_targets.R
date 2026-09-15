# --- predict_mirna_targets.R ---
suppressMessages({
  library(dplyr)
  library(multiMiR)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
deg_dir <- args[1]     
out_dir <- args[2]     

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting miRNA target prediction using multiMiR...")

csv_files <- list.files(deg_dir, pattern = "01_Significant_DEGs_.*_annotated\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No annotated significant DEG tables found.")
}

for (file in csv_files) {
  comp_name <- gsub("01_Significant_DEGs_", "", basename(file))
  comp_name <- gsub("_annotated\\.csv$", "", comp_name)
  
  df <- read.csv(file)
  
  # 1. Filter out NA values and grab only miRNAs
  df_mir <- df %>%
    filter(Regulation != "Not_Significant" & Biotype == "miRNA" & !is.na(Symbol) & Symbol != "")
  
  if (nrow(df_mir) == 0) {
    message(paste("  -> No valid significant miRNAs found in", comp_name, "- Skipping."))
    next
  }
  
  message(paste("\nProcessing", nrow(df_mir), "miRNAs for:", comp_name))
  
  # 2. Generate an expanded net of miRBase variants for each symbol
  base_names <- tolower(df_mir$Symbol)
  base_names <- sub("^mir", "", base_names) 
  
  mir_queries <- c()
  for (base in base_names) {
    mir_queries <- c(mir_queries,
                     paste0("mmu-mir-", base),       
                     paste0("mmu-miR-", base),       
                     paste0("mmu-miR-", base, "-5p"), 
                     paste0("mmu-miR-", base, "-3p"), 
                     paste0("mmu-miR-", base, "a"),   
                     paste0("mmu-miR-", base, "b"))
  }
  mir_queries <- unique(mir_queries)
  
  message("  -> Querying databases with expanded miRBase variants...")
  
  target_results <- tryCatch({
    # We use suppressWarnings to silence multiMiR's internal matrix math bugs
    suppressWarnings(
      get_multimir(org     = "mmu",
                   mirna   = mir_queries,
                   table   = "all", 
                   summary = TRUE,
                   predicted.cutoff.type = "p", 
                   predicted.cutoff      = 10)  
    )
  }, error = function(e) {
    message("  -> Database query error or no targets found.")
    return(NULL)
  })
  
  if (is.null(target_results) || nrow(target_results@data) == 0) {
    message("  -> No targets found in databases for these miRNAs.")
    next
  }
  
  # 3. Clean and organize the results (CRITICAL FIX: dplyr::select)
  res_df <- target_results@data %>%
    dplyr::select(mature_mirna_id, target_symbol, target_entrez, type, database) %>%
    dplyr::group_by(mature_mirna_id, target_symbol) %>%
    dplyr::summarise(
      Evidence_Type = paste(unique(type), collapse = ", "),
      Databases_Found = paste(unique(database), collapse = ", "),
      Hit_Count = n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(Hit_Count)) 
  
  # Save the master target table
  out_file <- file.path(out_dir, paste0("01_miRNA_Targets_", comp_name, ".csv"))
  write.csv(res_df, out_file, row.names = FALSE)
  
  message(paste("  -> Success! Found", nrow(res_df), "unique miRNA-Target interactions."))
  message(paste("  -> Saved to:", out_file))
}

message("\nAll miRNA target predictions completed!")
