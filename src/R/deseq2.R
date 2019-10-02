set.seed(1024)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

# -----------------------functions------------------------------
#' Plot PCA
#'
#' @param rlog_dds dds An rlog transformed deseq2 object
#' @param plot_out_path FILEPATH Path where to save pca plot
#'
#' @return A dataframe giving the positions of each sample on PC1 and PC2
#' @export
#'
#' @examples
#' plot_pca(rlog_dds = rlog_dds_obj, plot_out_path = /path/to/pca.pdf)
plot_pca <- function(rlog_dds, plot_out_path) {
  # Vars
  intgroup <- "condition"
  plot_title <- "PCA"
  
  # Getting the pca df from deseq2
  pca_dat <- DESeq2::plotPCA(object = rlog_dds, intgroup = intgroup, returnData = TRUE)
  
  # Setting up the pca colours
  col <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
           RColorBrewer::brewer.pal(n = 4, name = "Set1"))
  
  # Plotting the PCA
  pca_gg <- ggplot(data = pca_dat, aes(x = PC1, y = PC2, color = group, label = name)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = col) +
    ggrepel::geom_text_repel(size = 3) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "bottom") +
    cowplot::background_grid(major = "xy") +
    ggtitle(plot_title)
  
  # Save the plot
  cowplot::save_plot(filename = plot_out_path, base_height = 7, base_width = 7, plot = pca_gg)
  
  # Return the PCA df
  return(pca_dat)
}
#'
#'
#' Deseq2 constrasts
#'
#' @param case Name of the case condition
#' @param control Name of the control condition
#' @param deseq2_obj A deseq2 obj
#' @param annot Dataframe having addtional info to add to deseq2 res e.g. geneinfo, fpkm
#' @param outdir Directory to output xls files to
#' @param lfc_threshold The threshold of log2fc for which walds p-values will be calculated
#'
#' @return None
#' @export
#'
#' @examples
deseq2_res <- function(case, control, deseq2_obj, annot, outdir, lfc_threshold = 0.1){
  # Setup contrast
  contrast <- c("condition", case, control)
  
  # Get deseq results for current contrast
  dds_res <- DESeq2::results(object = deseq2_obj, contrast = contrast, lfcThreshold = lfc_threshold)
  
  # Shrink the logfoldchanges for more robust estimates
  dds_res <- DESeq2::lfcShrink(dds = deseq2_obj, contrast = contrast, res = dds_res)
  dds_res <- tibble::as_tibble(dds_res, rownames = "Ens_id")
  
  # Setup annot df to merge
  # Remove FPKM samples which are not relevant to current contrast
  samples_to_exclude <- dds$sample_name[!dds$condition %in% c(case, control)]
  if(length(samples_to_exclude) > 0){
    samples_to_exclude <- paste0(samples_to_exclude, "_FPKM")
  }
  annot <- dplyr::select(annot, -samples_to_exclude)
  
  # Merge annot with deseq2 res
  dds_res_merged <- dplyr::left_join(dds_res, annot, by = "Ens_id")
  dds_res_merged <- dplyr::arrange(dds_res_merged, padj)
  
  # Write output
  write.table(x=dds_res_merged, file=file.path(outdir, paste0(case, "_vs_", control, ".xls")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Read in the kallsito out files into a names vector
kal_out_files = snakemake@input[["counts"]]
kal_out_name = basename(dirname(kal_out_files))
kal_out_df <- data.frame(kallisto_out = kal_out_files, sample_name = kal_out_name, stringsAsFactors = FALSE)
#
# Read in the metadata df
metadata <- read.csv("metadata.tsv", sep = "\t", stringsAsFactors = FALSE)
metadata <- merge(metadata, kal_out_df, by = "sample_name")
metadata$condition <- as.factor(metadata$condition)
rownames(metadata) <- metadata$sample_name
#
# Create tximport object
txi <- tximport::tximport(
  files = metadata$kallisto_out, type = "kallisto", txIn = TRUE,
  tx2gene = readr::read_tsv(snakemake@input[["tx2gene"]])
)
#
# Create deseq2 object
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = metadata, design = ~condition)
#
# Estimate size factors , estimate dispersions and perform nbinomWaldTest
dds <- DESeq2::DESeq(object = dds, parallel = TRUE)
#
#
# Get information about genes
gene_info <- read.csv(file = snakemake@input[["geneinfo"]], sep = "\t", stringsAsFactors = FALSE)
#
# Get fpkm values for downstream processing and tag them with _FPKM
fpkm <-  as.data.frame(DESeq2::fpkm(dds, robust = TRUE), stringsAsFactors = FALSE)
colnames(fpkm) <- paste0(colnames(fpkm), "_FPKM")
fpkm$Ens_id <- rownames(fpkm)
#
annot_df <- merge(gene_info, fpkm, by = "Ens_id")
#
# Perform contrast wise ops
# Transform all contrasts into a list with each element being a vector of 2 elements [1] case [2] control
contrasts <- unlist(strsplit(snakemake@params[["contrasts"]], ",")) %>%
  purrr::map(~unlist(strsplit(.x, "_vs_")))
#
purrr::map(contrasts, ~deseq2_res(case = .x[1],
                                  control = .x[2],
                                  deseq2_obj = dds,
                                  annot = annot_df,
                                  outdir = "results/deg/"))
#
# Get rlog transformed deseq2 object,  used for pca
rld <- DESeq2::rlog(object = dds, blind = FALSE)
#
# Plot the pca
pca_dat <- plot_pca(rlog_dds = rld, plot_out_path = "results/deg/pca.pdf")


# Reproducibility
sessionInfo()

