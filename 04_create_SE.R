# Create a Summarized Experiment object a user supplied annotation file
# and the quant.sf files for each sample processed by the REdiscoverTE pipeline.
#
# NOTE: There are hard-coded paths to annotation files
# ------------------------------------------------------------------------------
suppressMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(SummarizedExperiment, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))

# get commandline arguments
option_list = list(
  make_option(c("-d", "--quants_dir"), 
              type = "character", 
              default = NULL, 
              help = "Path to quants/ directory containing sub-directories for each SAMPLE", 
              metavar = "quants_dir"),
  make_option(c("-f", "--annotation"),
              type = "character",
              default = NULL,
              help = "Path to annotation file. Must at least contain a column called 'sample_name' and 'group'. Other metadata is optional",
              metavar = "annotation"),
  make_option(c("-o", "--out_dir"),
              type = "character",
              default = ".",
              help = "Location to save the final SE object RDS file",
              metavar = "out_dir"),
  make_option(c("-c", "--cores"),
              type = "integer",
              default = 4,
              help = "Number of cores to use when processing quant files",
              metavar = "cores")
  )

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

# Function used to process a single quant file ----------------------------------
# annotation files used in process function
message("Reading in annotation files for processing steps...")
gencove_annot <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/REdiscoverTE_hg38_GFP/GENCODE.V26.Basic_Gene_Annotation_md5_GFP.rds"))
rmsk_annot <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/rmsk_annotation.RDS"))

process_quant_file <- function(quant_file) {
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

# main processing steps --------------------------------------------------------
# list all quant files in quants dir
message(paste0("Reading in list of quant files from ", opt$quants_dir, "..."))
quant_files <- list.files(opt$quants_dir,
                          pattern = "quant.sf",
                          recursive = TRUE,
                          full.names = TRUE)

# name each quant file with its SRA accession
names(quant_files) <- regmatches(quant_files, regexpr("SRR[0-9]+", quant_files))

# process the quant files in parallel
message(paste("Processing quant files using", opt$cores, "cores..."))
processed <- parallel::mclapply(
  X = quant_files,
  FUN = process_quant_file,
  mc.cores = opt$cores)

# extract gene and RE counts as dataframes
gene_counts <- map(processed, purrr::pluck, "gene")
re_counts <- map(processed, purrr::pluck, "re")

# collapse into single dataframe and round the counts
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
message("Reading in GTF file for gene annotations...")
cds <- readRDS("/home/gcalendo/data/projects/geo-rna-seq/doc/gencode.v26.annotation.cds.rds")

# filter counts matrix for only protein coding genes
message("Filtering for protein coding genes...")
gene_mat <- subset(gene_mat, rownames(gene_mat) %in% cds$gene_name)

# Check that the same samples exist in the RE matrix and gene matrix
message("Checking colnames align in gene and RE matrices...")
if (all(colnames(gene_mat) == colnames(re_mat))) {
  counts <- rbind(gene_mat, re_mat)
} else {
  stop("colnames in gene matrix do not match colnames in RE matrix.")
}

# Read in the annotation information
message("Reading in sample annotation information...")
annotation_path <- opt$annotation
ext <- tools::file_ext(annotation_path)
annot_df <- switch(ext,
                   csv = vroom::vroom(annotation_path, delim = ','),
                   tsv = vroom::vroom(annotation_path, delim = '\t'),
                   txt = vroom::vroom(annotation_path, delim = '\t')
                   )

# check that 'sample_name' and 'group' are present in annotation
stopifnot('sample_name' %in% colnames(annot_df))
stopifnot('group' %in% colnames(annot_df))

# convert the sample_name column to rownames
message("Converting annotation from tibble to data.frame...")
col_data <- annot_df %>%
  column_to_rownames(var = "sample_name")

# coerce the group column to a factor
message("Coercing 'group' to be a factor variable...")
col_data$group <- factor(col_data$group)

# reorder colnames of counts to match rownames of annotation
message("Reordering the count matrix by sample metadata...")
counts <- counts[, rownames(col_data)]

# check that order of colnames in counts is same as rownames in col_data
if (!(all(rownames(col_data) == colnames(counts)))) {
  stop("colnames of counts matrix and rownames of sample annotation do not match.")
}

# create SummarizedExperiment object
message("Creating Summarized Experiment object...")
se <- SummarizedExperiment(assays = list(counts = as.matrix(counts)), 
                           colData = col_data)

# write se RDS file to out_dir
message(paste("Writing SummarizedExperiment object to", file.path(opt$out_dir, "summarized_experiment.rds")))
saveRDS(se, file = file.path(opt$out_dir, "summarized_experiment.rds"))
message("Done.")
