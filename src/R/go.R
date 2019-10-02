log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(goseq))


# Read deg output
deg_xls <- readr::read_tsv(snakemake@input[["degxls"]])

# Remove lowly expressed genes
deg_xls <- dplyr::filter(deg_xls, baseMean > snakemake@params[["basemean_cutoff"]])

# goseq has problems with ensids with version so remove them if present
deg_xls$Ens_id <- gsub("\\.\\d+$", "", deg_xls$Ens_id)


# Subset to significant genes
sig_genes <- dplyr::filter(deg_xls,  padj < snakemake@params[["padj_cutoff"]])

# Create a vector of deg yes/no required by goseq
de_binary_vect <- as.integer(deg_xls$Ens_id %in% sig_genes$Ens_id)
names(de_binary_vect) <- deg_xls$Ens_id

# Go seq analysis
pwf <- goseq::nullp(DEgenes = de_binary_vect, genome = snakemake@params[["genome"]], id = "ensGene")
go <- goseq::goseq(pwf, genome = snakemake@params[["genome"]], id = "ensGene", test.cats = c("GO:BP", "GO:MF"))

# Do fdr adjustment
go$fdr <- p.adjust(go$over_represented_pvalue)

# Write output
write.table(go, file = snakemake@input[["gores"]], sep = "\t", row.names = FALSE)

# Reproducibility
sessioninfo::session_info()
