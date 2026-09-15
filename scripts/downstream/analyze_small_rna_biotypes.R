# --- analyze_small_rna_biotypes.R ---
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(broom)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file <- args[2]
biotype_file <- args[3]
deg_dir <- args[4]
out_dir <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
message("Loading data for small RNA biotype analysis...")

counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)
biotype_map <- read.table(biotype_file, header=FALSE, sep="\t", col.names=c("GeneID", "Biotype"), fill=TRUE)

common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]
meta$Sample <- rownames(meta)
meta$Group <- factor(paste(meta$Treatment, meta$Generation, sep="_"))

# --- 2. Process Raw Counts for Global Biotypes ---
counts$GeneID <- rownames(counts)

long_counts <- counts %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Count") %>%
  left_join(biotype_map, by = "GeneID") %>%
  left_join(meta, by = "Sample") %>%
  mutate(Biotype = ifelse(is.na(Biotype) | Biotype == "", "Unknown", Biotype))

sample_biotypes <- long_counts %>%
  group_by(Sample, Treatment, Generation, Group, Biotype) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percentage = (TotalCount / sum(TotalCount)) * 100) %>%
  ungroup()

group_biotypes <- sample_biotypes %>%
  group_by(Group, Biotype) %>%
  summarise(Mean_Percentage = mean(Percentage), .groups = "drop") %>%
  mutate(Biotype_Clean = ifelse(Mean_Percentage < 1, "Other (<1%)", Biotype)) %>%
  group_by(Group, Biotype_Clean) %>%
  summarise(Mean_Percentage = sum(Mean_Percentage), .groups = "drop")

# --- TASK 1A: Plot Global Biotype Pie Charts (Individual Plots with Percentages in Legend) ---
message("Generating individual global biotype pie charts...")

# Create a consistent color palette across all plots so "miRNA" is always the same color
unique_biotypes <- unique(group_biotypes$Biotype_Clean)
getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
base_colors <- getPalette(length(unique_biotypes))
names(base_colors) <- unique_biotypes

for (grp in unique(group_biotypes$Group)) {
  df_grp <- group_biotypes %>% filter(Group == grp)
  
  # Create a custom label that includes the percentage
  df_grp <- df_grp %>%
    mutate(LegendLabel = paste0(Biotype_Clean, " (", round(Mean_Percentage, 1), "%)"))
  
  # Map the dynamic legend labels to the base colors
  grp_colors <- base_colors[df_grp$Biotype_Clean]
  names(grp_colors) <- df_grp$LegendLabel
  
  p_pie <- ggplot(df_grp, aes(x = "", y = Mean_Percentage, fill = LegendLabel)) +
    geom_bar(stat = "identity", width = 1, color="white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = grp_colors) +
    theme_void() +
    labs(title = paste("Average Biotype Distribution:", grp), fill = "Biotype") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust=0.5),
          legend.title = element_text(face="bold", size=12),
          legend.text = element_text(size=11))
  
  ggsave(file.path(out_dir, paste0("01_Global_Biotypes_Pie_", grp, ".png")), plot = p_pie, width = 8, height = 6)
}

# --- TASK 1B: Statistics on Global Biotypes ---
message("Calculating Two-Way ANOVA for Biotype composition shifts...")

top_biotypes <- group_biotypes %>% group_by(Biotype_Clean) %>% summarise(Total = sum(Mean_Percentage)) %>% arrange(desc(Total)) %>% head(5) %>% pull(Biotype_Clean)
stats_results <- data.frame()

for (bio in top_biotypes) {
  if (bio == "Other (<1%)") next
  df_stat <- sample_biotypes %>% filter(Biotype == bio)
  model <- aov(Percentage ~ Treatment * Generation, data = df_stat)
  tidy_model <- tidy(model) %>% filter(term != "Residuals") %>% mutate(Biotype = bio)
  stats_results <- bind_rows(stats_results, tidy_model)
}

stats_results <- stats_results %>%
  select(Biotype, term, statistic, p.value) %>%
  mutate(Significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  arrange(Biotype, p.value)

write.csv(stats_results, file.path(out_dir, "02_Biotype_ANOVA_Statistics.csv"), row.names = FALSE)
# --- TASK 2: Biotypes of Differentially Expressed Genes (DEGs) ---
message("Generating individual Pie Charts for Differentially Expressed Genes with consistent colors...")

# Locate all annotated significant DEG files generated in previous steps
csv_files <- list.files(deg_dir, pattern = "01_Significant_DEGs_.*_annotated\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  message("No annotated significant DEG tables found in the specified directory.")
} else {
  
  # STEP 2A: Pre-scan all files to extract ALL unique DEG biotypes
  # This ensures we build a master palette so colors remain consistent across Up/Down and different generations
  all_deg_biotypes <- c()
  deg_data_list <- list()
  
  for (file in csv_files) {
    df_deg <- read.csv(file)
    if (nrow(df_deg) > 0) {
      df_deg <- df_deg %>%
        mutate(Biotype = ifelse(is.na(Biotype) | Biotype == "", "Unknown", Biotype))
      
      all_deg_biotypes <- c(all_deg_biotypes, unique(df_deg$Biotype))
      deg_data_list[[file]] <- df_deg # Store the dataframe in memory to prevent re-reading
    }
  }
  
  # STEP 2B: Generate the master color palette for DEGs
  unique_deg_biotypes <- sort(unique(all_deg_biotypes))
  getPalette_deg <- colorRampPalette(brewer.pal(8, "Set2"))
  deg_base_colors <- getPalette_deg(length(unique_deg_biotypes))
  names(deg_base_colors) <- unique_deg_biotypes
  
  # STEP 2C: Generate the pie charts applying the master palette
  for (file in names(deg_data_list)) {
    # Extract the comparison name from the file name (e.g., "iAs_vs_Control_in_F2")
    comp_name <- gsub("01_Significant_DEGs_", "", basename(file))
    comp_name <- gsub("_annotated\\.csv$", "", comp_name)
    
    df_deg <- deg_data_list[[file]]
    
    # Calculate counts and percentages of biotypes within Upregulated and Downregulated subsets
    deg_biotypes <- df_deg %>%
      group_by(Regulation, Biotype) %>%
      summarise(GeneCount = n(), .groups = "drop") %>%
      group_by(Regulation) %>%
      mutate(Percentage = (GeneCount / sum(GeneCount)) * 100) %>%
      ungroup()
    
    # Generate separate pie charts for Upregulated and Downregulated DEGs
    for (reg in unique(deg_biotypes$Regulation)) {
      df_reg <- deg_biotypes %>% filter(Regulation == reg)
      
      # Create a complex label containing Biotype, Count, and Percentage
      df_reg <- df_reg %>%
        mutate(LegendLabel = paste0(Biotype, " (Count: ", GeneCount, ", ", round(Percentage, 1), "%)"))
      
      # Map the master colors to the newly generated dynamic labels
      grp_colors_deg <- deg_base_colors[df_reg$Biotype]
      names(grp_colors_deg) <- df_reg$LegendLabel
      
      p_deg <- ggplot(df_reg, aes(x = "", y = GeneCount, fill = LegendLabel)) +
        geom_bar(stat = "identity", width = 1, color="white") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = grp_colors_deg) + # Apply the consistent master palette
        theme_void() +
        labs(title = paste(reg, "DEG Biotypes\n", gsub("_", " ", comp_name)), fill = "Biotype") +
        theme(plot.title = element_text(size = 14, face = "bold", hjust=0.5),
              legend.title = element_text(face="bold", size=12),
              legend.text = element_text(size=11))
      
      ggsave(file.path(out_dir, paste0("03_DEG_Biotypes_Pie_", comp_name, "_", reg, ".png")), plot = p_deg, width = 8, height = 6)
    }
  }
}

message("\nSmall RNA Biotype analysis completed successfully!")