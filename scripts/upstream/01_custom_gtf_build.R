# --- Load Libraries ---
# We need 'rtracklayer' to read/export GFF/GTF files
# and 'dplyr' for combining the data frames.
library(rtracklayer)
library(dplyr) # <-- Changed from 'plyr'

# --- Define File Names ---
gencode_file <- "gencode.vM25.annotation.gtf.gz"
trna_file <- "mm10.tRNAscan.gtf.gz"
pirna_file <- "pirnadb.v1_7_6.mm10.gtf.gz"

output_file <- "mm10.RNA_full.custom.gtf"

# --- Step 1: Load the files ---
print(paste("Loading Gencode (base)...", gencode_file))
genome <- readGFF(gencode_file)
print("Gencode loaded.")

print(paste("Loading tRNAs...", trna_file))
tRNAs <- readGFF(trna_file)
print("tRNAs loaded.")

print(paste("Loading piRNAs...", pirna_file))
piRNAs <- readGFF(pirna_file)
print("piRNAs loaded.")

# --- Step 2: Combine the files ---
print("Combining files...")

# 
# bind_rows() is the dplyr equivalent of rbind.fill()
# It's also cleaner, as we can list all data frames in one call.
genome_complete <- bind_rows(genome, tRNAs, piRNAs)
# This one line replaces the two rbind.fill() lines

print("Files combined.")

# --- Step 3: Export the final GTF ---
print(paste("Exporting final GTF to:", output_file))
export(genome_complete, output_file, format = 'gtf')

print("SUCCESS!")
print(paste("Your custom GTF has been created:", output_file))