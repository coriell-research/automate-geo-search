# Create a MultiAssayExperiment object a user supplied annotation file
# and the quant.sf files for each sample processed by the REdiscoverTE pipeline.
#
# NOTE: There are hard-coded paths to annotation files
# ------------------------------------------------------------------------------
suppressMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(MultiAssayExperiment, warn.conflicts = FALSE, quietly = TRUE))

# get commandline arguments
option_list <- list(
  make_option(c("-d", "--quants_dir"),
    type = "character",
    default = NULL,
    help = "Path to quants/ directory containing sub-directories for each SAMPLE",
    metavar = "quants_dir"
  ),
  make_option(c("-f", "--annotation"),
    type = "character",
    default = NULL,
    help = "Path to annotation file. Must at least contain a column called 'sample_name' and 'group'. Other metadata is optional",
    metavar = "annotation"
  ),
  make_option(c("-o", "--out_dir"),
    type = "character",
    default = ".",
    help = "Location to save the final MultiAssayExperiment object RDS file",
    metavar = "out_dir"
  ),
  make_option(c("-c", "--cores"),
    type = "integer",
    default = 4,
    help = "Number of cores to use when processing quant files",
    metavar = "cores"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

message("Reading in annotation files for processing steps...")
gencode_dt <- readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/gencode_dt.rds")
exon_dt <- readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/exon_dt.rds")
intron_dt <- readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/intron_dt.rds")
intergenic_dt <- readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/intergenic_dt.rds")

# Function for processing a single quant file
process_quant_file <- function(quant_file) {
  DT <- fread(quant_file, sep = "\t")
  DT <- DT[NumReads > 0]

  # Inner join the gene annotation to the quants
  transcript_counts <-
    merge.data.table(
      x = DT,
      y = gencode_dt,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  # summarize to the gene name level and round to integers
  gene_counts <- transcript_counts[, .(count = round(sum(NumReads, na.rm = TRUE), digits = 0)), by = .(symbol)]

  # Inner join counts with exonic RE annotation ---
  sample_exon_dt <-
    merge.data.table(
      x = DT,
      y = exon_dt,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  exon_counts <- sample_exon_dt[
    ,
    repElem := paste(repClass, repFamily, repName, sep = ".")
  ][,
    .(count = sum(NumReads, na.rm = TRUE)),
    by = "repElem"
  ][
    ,
    count := round(count, digits = 0)
  ]

  # Inner join counts with Intronic annotation ---
  sample_intron_dt <-
    merge.data.table(
      x = DT,
      y = intron_dt,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  intron_counts <- sample_intron_dt[
    ,
    repElem := paste(repClass, repFamily, repName, sep = ".")
  ][,
    .(count = sum(NumReads, na.rm = TRUE)),
    by = "repElem"
  ][
    ,
    count := round(count, digits = 0)
  ]

  # Inner join counts with Intronic annotation ---
  sample_intergenic_dt <-
    merge.data.table(
      x = DT,
      y = intergenic_dt,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  intergenic_counts <- sample_intergenic_dt[
    ,
    repElem := paste(repClass, repFamily, repName, sep = ".")
  ][,
    .(count = sum(NumReads, na.rm = TRUE)),
    by = "repElem"
  ][
    ,
    count := round(count, digits = 0)
  ]

  # combine the intronic and intergenic counts into a single data.table ---
  re_counts <- rbindlist(
    list(intron_counts, intergenic_counts)
  )[,
    .(count = sum(count)),
    by = "repElem"
  ]

  data_list <- list(
    "gene" = gene_counts,
    "exon" = exon_counts,
    "intron" = intron_counts,
    "intergenic" = intergenic_counts,
    "re" = re_counts
  )

  return(data_list)
}

# main processing steps --------------------------------------------------------
message(paste0("Reading in list of quant files from ", opt$quants_dir, "..."))
quant_files <- list.files(opt$quants_dir,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

# name quant files by SRA accession from filename
names(quant_files) <- regmatches(quant_files, regexpr("SRR[0-9]+", quant_files))

# Run process_quant_file func over all quant files in parallel
message(paste("Processing quant files using", opt$cores, "cores..."))
processed <- parallel::mclapply(
  X = quant_files,
  FUN = process_quant_file,
  mc.cores = opt$cores
)

# extract individual data.tables of counts
gene_dts <- lapply(processed, function(x) x[["gene"]])
re_dts <- lapply(processed, function(x) x[["re"]])
exon_dts <- lapply(processed, function(x) x[["exon"]])
intron_dts <- lapply(processed, function(x) x[["intron"]])
intergenic_dts <- lapply(processed, function(x) x[["intergenic"]])

# collapse each into single data.tables
gene_dt <- rbindlist(gene_dts, idcol = "sample_name")
re_dt <- rbindlist(re_dts, idcol = "sample_name")
exon_dt <- rbindlist(exon_dts, idcol = "sample_name")
intron_dt <- rbindlist(intron_dts, idcol = "sample_name")
intergenic_dt <- rbindlist(intergenic_dts, idcol = "sample_name")

# cast each dt wider and coerce to matrices
gene_mat <- as.matrix(dcast(gene_dt, symbol ~ sample_name, value.var = "count", fill = 0), rownames = "symbol")
re_mat <- as.matrix(dcast(re_dt, repElem ~ sample_name, value.var = "count", fill = 0), rownames = "repElem")
exon_mat <- as.matrix(dcast(exon_dt, repElem ~ sample_name, value.var = "count", fill = 0), rownames = "repElem")
intron_mat <- as.matrix(dcast(intron_dt, repElem ~ sample_name, value.var = "count", fill = 0), rownames = "repElem")
intergenic_mat <- as.matrix(dcast(intergenic_dt, repElem ~ sample_name, value.var = "count", fill = 0), rownames = "repElem")

# check that all quant files were processed
stopifnot(all(names(quant_files) %in% colnames(gene_mat)))
stopifnot(all(names(quant_files) %in% colnames(re_mat)))
stopifnot(all(names(quant_files) %in% colnames(exon_mat)))
stopifnot(all(names(quant_files) %in% colnames(intron_mat)))
stopifnot(all(names(quant_files) %in% colnames(intergenic_mat)))

# reorder the matrices by the quant file names
gene_mat <- gene_mat[, names(quant_files)]
re_mat <- re_mat[, names(quant_files)]
exon_mat <- exon_mat[, names(quant_files)]
intron_mat <- intron_mat[, names(quant_files)]
intergenic_mat <- intergenic_mat[, names(quant_files)]

# Read in the sample metadata
message("Reading in sample annotation information...")
meta_df <- fread(opt$annotation)
stopifnot("sample_name" %in% colnames(meta_df))
stopifnot("group" %in% colnames(meta_df))

# convert to data.frame for colData
setDF(meta_df, rownames = meta_df$sample_name)
meta_df <- subset(meta_df, select = -sample_name)
stopifnot(all(names(quant_files) %in% rownames(meta_df)))

# reorder the metadata to match colnames of matrices
meta_df <- meta_df[names(quant_files), ]

# create MultiAssayExperiment object from each assay and metadata
exp_list <- list(
  "gene counts" = gene_mat,
  "RE counts" = re_mat,
  "exon RE counts" = exon_mat,
  "intron RE counts" = intron_mat,
  "intergenic RE counts" = intergenic_mat
)

message("Creating MultiAssayExperiment object...")
mae <- MultiAssayExperiment(experiments = exp_list, colData = meta_df)

# write mae .rds file to out_dir
message(paste("Writing MultiAssayExperiment object to", file.path(opt$out_dir, "MultiAssayExperiment.rds")))
saveRDS(mae, file = file.path(opt$out_dir, "MultiAssayExperiment.rds"))
message("Done.")
