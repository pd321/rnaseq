log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(goseq))

# Read in the deg file
rnk_df <- readr::read_tsv(snakemake@input[["degxls"]])

# Get mean fpkm
rnk_df$meanfpkm <- rowMeans(dplyr::select(rnk_df, dplyr::ends_with("_FPKM")), na.rm = TRUE)

# Filter for mean fpkm
rnk_df <- dplyr::filter(rnk_df, meanfpkm >= snakemake@params[["mean_fpkm_cutoff"]])

# Keep only relevant cols
rnk_df <- dplyr::select(rnk_df, GeneName, log2FoldChange)

# Make gene names uppercase
rnk_df <- dplyr::mutate(rnk_df, GeneName = toupper(GeneName))

# Sort by log2fc
rnk_df <- dplyr::arrange(rnk_df, dplyr::desc(log2FoldChange))

# Write output
write.table(x=rnk_df, file=snakemake@output[["gsearnk"]], sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Reproducibility
sessioninfo::session_info()
