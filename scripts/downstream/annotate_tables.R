# --- annotate_tables.R ---
# Script independiente para anotar tablas de DESeq2 con Symbol y Biotype

suppressMessages({
  library(dplyr)
  library(org.Mm.eg.db)
})

fractions <- c("large", "small")

for (frac in fractions) {
  message(paste("\n--- Anotando fracción:", toupper(frac), "---"))
  
  res_dir <- file.path("results", paste0(frac, "_fraction"), "deseq2_apeglm_all_samples")
  biotype_file <- file.path(paste0(frac, "_fraction"), "data", "gene_biotype_map.tsv")
  
  if (!file.exists(biotype_file)) {
    warning(paste("No se encontró el archivo de biotipos:", biotype_file))
    next
  }
  biotype_map <- read.table(biotype_file, header=FALSE, sep="\t", 
                            col.names=c("GeneID", "Biotype"), fill=TRUE)
  
  csv_files <- list.files(res_dir, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("_annotated\\.csv$", csv_files)]
  
  for (file in csv_files) {
    message(paste("  Procesando:", basename(file)))
    
    df <- read.csv(file)
    
    if (!"GeneID" %in% colnames(df)) {
      next
    }
    
    # --- CORRECCIÓN: ¿Qué pasa si la tabla está vacía (0 genes DE)? ---
    if (nrow(df) == 0) {
      message("    -> Tabla vacía (0 genes). Guardando estructura vacía...")
      df$Symbol <- character()
      df$Biotype <- character()
      out_file <- sub("\\.csv$", "_annotated.csv", file)
      write.csv(df, out_file, row.names = FALSE)
      next # Saltar al siguiente archivo
    }
    # -------------------------------------------------------------------

    clean_ids <- gsub("\\..*$", "", df$GeneID)
    
    symbols <- mapIds(org.Mm.eg.db,
                      keys = clean_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
    
    df$Symbol <- symbols
    
    # Unir con el mapa de biotipos
    df <- left_join(df, biotype_map, by = "GeneID")
    
    # Reordenar columnas (ID, Symbol, Biotype al frente)
    cols <- colnames(df)
    front_cols <- c("GeneID", "Symbol", "Biotype")
    remaining_cols <- setdiff(cols, front_cols)
    df <- df[, c(front_cols, remaining_cols)]
    
    out_file <- sub("\\.csv$", "_annotated.csv", file)
    write.csv(df, out_file, row.names = FALSE)
  }
}

message("\n¡Anotación completada con éxito!")
