# Create a Summarized Experiment object a user supplied annotation file
# and the quant.sf files for each sample processed by the REdiscoverTE pipeline.
#
# NOTE: path to REdiscoverTE RDS and rmsk Annotation is hard-coded
# ------------------------------------------------------------------------------
library(optparse)
library(data.table)
library(SummarizedExperiment)
library(tidyverse)


# get command line args ---------------------------------------------------
option_list = list(
  make_option(c("-d", "--project_dir"), 
              type = "character", 
              default = NULL, 
              help = "Path to PROJECT directory containing quants dir and annotation file", 
              metavar = "project_dir"),
  make_option(c("--cores"),
              type = "integer",
              default = 4,
              help = "Number of cores to use when processing quant files",
              metavar = "cores"),
  make_option(c("--rediscover_rds"),
              type = "character",
              default = "/mnt/data/gdata/human/REdiscoverTE_hg38/REdiscoverTE_hg38_GFP/GENCODE.V26.Basic_Gene_Annotation_md5_GFP.rds",
              help = "path to REdiscoverTE annotation rds file",
              metavar = "rediscover_rds"),
  make_option(c("--rmsk_rds"),
              type = "character",
              default = "/mnt/data/gdata/human/REdiscoverTE_hg38/rmsk_annotation.RDS",
              help = "path to repeat masker rds annotation file",
              metavar = "rmsk_rds"),
  make_option(c("--gtf"),
              type = "character",
              default = "/home/gcalendo/data/projects/yb5-combined-analysis/doc/gencode.v26.annotation.gtf.gz",
              help = "path to gtf gene annotation file",
              metavar = "gtf")
  )

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)


# Function used to process quant files ------------------------------------
process_quant_file <- function(quant_file, rediscover_annot, rmsk_annot) {
  gencove_annot <- as.data.table(readRDS(rediscover_annot))
  rmsk_annot <- as.data.table(readRDS(rmsk_annot))
 
  DT <- data.table::fread(quant_file, sep = "\t")
  DT <- DT[NumReads > 0]
  
  transcript_counts <- data.table::merge.data.table(
    x = DT,
    y = gencove_annot,
    by.x = "Name",
    by.y = "md5"
  )
  
  repElem_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_annot,
    by.x = "Name",
    by.y = "md5"
  )
 
  gene_counts <- transcript_counts[,.(count = sum(NumReads, na.rm = TRUE)), by = .(symbol)]
  repElem_counts <- repElem_counts[, .(repElem = paste(repClass, repFamily, repName, sep = "."), NumReads)][, .(count = sum(NumReads, na.rm = TRUE)), by = .(repElem)]
  
  list("gene" = gene_counts, "re" = repElem_counts)
}


# list and name all quant files
quant_files <- list.files(opt$project_dir,
                          pattern = "quant.sf",
                          recursive = TRUE,
                          full.names = TRUE)

names(quant_files) <- str_extract(quant_files, "SRR[0-9]+")

# process the quant files
processed <- parallel::mclapply(
  X = quant_files, 
  FUN = process_quant_file,
  opt$rediscover_rds,
  opt$rmsk_rds,
  mc.cores = opt$cores)

# extract gene and RE counts as dataframes
gene_counts <- map(processed, purrr::pluck, "gene")
re_counts <- map(processed, purrr::pluck, "re")

# collapse into single dataframe and round the counts for edgeR
gene_df <- bind_rows(gene_counts, .id = "sample_name") %>% mutate(count = round(count))
re_df <- bind_rows(re_counts, .id = "sample_name") %>% mutate(count = round(count))

# coerce count dataframes to count matrices
gene_mat <- gene_df %>% 
  pivot_wider(names_from = sample_name, 
              values_from = count, 
              values_fill = 0) %>% 
  column_to_rownames(var = "symbol")

re_mat <- re_df %>% 
  pivot_wider(names_from = sample_name,
              values_from = count,
              values_fill = 0) %>% 
  column_to_rownames(var = "repElem")

# read in GTF file for gene annotations
gtf <- as.data.frame(rtracklayer::import(here("doc", "gencode.v26.annotation.gtf.gz")))

# filter for protein coding genes in gene matrix
cds <- subset(gtf, 
              type == "gene" & gene_type == "protein_coding", 
              c("gene_id", "gene_name", "seqnames", "start", "end", "width", "strand")
)

# filter counts matrix for only protein coding genes
gene_mat <- subset(gene_mat, rownames(gene_mat) %in% cds$gene_name)
