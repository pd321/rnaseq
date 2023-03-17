set.seed(1024)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(PCAtools))

# -----------------------functions------------------------------
#' Plot PCA
#'
#' @param rlog_dds dds An rlog transformed deseq2 object
#' @param deseq2_obj A deseq2 obj
#' @param plot_out_path FILEPATH Path where to save pca plot
#'
#' @return A dataframe giving the positions of each sample on PC1 and PC2
#' @export
#'
#' @examples
#' plot_pca(rlog_dds = rlog_dds_obj, plot_out_path = /path/to/pca.pdf)
plot_pca <- function(rlog_dds, deseq2_obj, plot_out_path) {

  # Getting the pca obj using pcatools
  p <- PCAtools::pca(assay(rlog_dds), metadata = colData(deseq2_obj))

  pdf(file = plot_out_path, width = 10, height = 10)
  print(PCAtools::biplot(p, colby = 'condition', legendPosition = 'right'))
  dev.off()

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
deseq2_res <- function(case, control, deseq2_obj, annot, outdir){
  # Setup contrast
  contrast <- c("condition", case, control)
  
  # Get deseq results for current contrast
  dds_res <- DESeq2::results(object = deseq2_obj, contrast = contrast)
  
  # Shrink the logfoldchanges for more robust estimates
  dds_res <- DESeq2::lfcShrink(dds = deseq2_obj, coef = paste0("condition_", case, "_vs_", control), res = dds_res)
  dds_res <- tibble::as_tibble(dds_res, rownames = "GENEID")
  
  # Setup annot df to merge
  # Remove FPKM samples which are not relevant to current contrast
  samples_to_exclude <- dds$sample_name[!dds$condition %in% c(case, control)]
  if(length(samples_to_exclude) > 0){
    samples_to_exclude <- paste0(samples_to_exclude, "_FPKM")
  }
  annot <- dplyr::select(annot, -samples_to_exclude)
  
  # Merge annot with deseq2 res
  dds_res_merged <- dplyr::left_join(dds_res, annot, by = "GENEID")
  dds_res_merged <- dplyr::arrange(dds_res_merged, padj)
  
  # Write output
  write.table(x=dds_res_merged, file=file.path(outdir, paste0(case, "_vs_", control, ".xls")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Read in the kallsito out files into a names vector
salmon_out_files = snakemake@input[["counts"]]
salmon_out_names = basename(dirname(salmon_out_files))
salmon_out_df <- data.frame(salmon_out = salmon_out_files, sample_name = salmon_out_names, stringsAsFactors = FALSE)
#
# Read in the metadata df
metadata <- read.csv(snakemake@input[["metadata_file"]], sep = "\t", stringsAsFactors = FALSE)
metadata <- merge(metadata, salmon_out_df, by = "sample_name")
metadata$condition <- as.factor(metadata$condition)
rownames(metadata) <- metadata$sample_name
#
# Create tximport object
txi <- tximport::tximport(
  files = metadata$salmon_out, type = "salmon", txIn = TRUE,
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
fpkm$GENEID <- rownames(fpkm)
#
annot_df <- dplyr::left_join(gene_info, fpkm, by = "GENEID")
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
pca_dat <- plot_pca(rlog_dds = rld, deseq2_obj = dds, plot_out_path = "results/deg/pca.pdf")


# Reproducibility
sessionInfo()
