# GEO Search

Automate downloading and mapping of RNA-Seq reads from cancer datasets present in GEO. 

## Step 1. Query GeoDataSets for studies of interest

The first argument to the script is the query string, the second is the desired out directory for the results file. 
The query string used should be the same formatting used when searching NCBI. Be sure to escape special characters
in the query. 

The script returns a file called 'search-results.tsv' written in the specified `$OUT_DIRECTORY`. The search-results.tsv 
file contains separate lines with the BIOPROJECT ID, the number of samples, and the title of the project. 
**Skim through this file and delete the lines that you do not want to process before proceeding**

```bash
./01_search.sh \
    "(\"expression profiling by high throughput sequencing\"[DataSet Type]) AND decitabine AND \"Homo sapiens\"[porgn]" \
    $OUT_DIRECTORY
```

## Step 2. Get the run information for each project

**NOTE: If you have already identified a dataset of interest the get_run_info_by_PRJNA.sh 
script takes a BioProject ID as input and returns run-info and metadata**. This script
is generally easier to use than the script below.

The next step parses the search-results.tsv file and pulls the standard run-info.csv 
file and an additional metadata.txt file from SRA. The script will output each 
file in a new, separate directory named with the BIOPROJECT ID for each record present in the 
search results. **To protect against large BIOPROJECTS being imported, the script limits processing to BIOPROJECTS with fewer than (100) samples** 
This can be changed by modifying the variables in the `02_get_run_info.sh` script.

```bash
./02_get_run_info.sh path/to/search-results.tsv $OUT_DIRECTORY
```

## Step 2b. Manually clean sample metadata

Since the run-info.csv file often does not contain useful sample metadata
information an additional metadata.txt file is downloaded into the same directory.
This file will contain SRA accessions along with associated metadata but is 
poorly formatted and needs to be cleaned before any automated processing. 

I typically create a copy of the metadata.txt file and name it `annotation.txt` (csv or tsv is also acceptable).
Within this annotation file I create new columns for the `sample_name` (SRA accession)
and the `group` (the grouping factor used in the experiment). **These two columns are neccessary for further processing**.
Additional metadata columns can also be present but `sample_name` and `group` are required 
if you want to automatically create `SummarizedExperiment` objects as described below (Step 4).

## Step 3. Download fastqs and map with Salmon

This step downloads all SRA accessions for each BIOPROJECT, maps each sample, 
then removes each fastq file after successful completion. 

To run downloading and mapping on all BIOPROJECTS downloaded:

```bash
for PROJECT in $(find . -name "PRJNA*" -type d); do
    python 03_download_and_process.py $PROJECT/run-info.csv;
done
```

By default, the script will download the fastq files from SRA into the PROJECT folder and then perform quantification
with Salmon, placing the individual quant results in `${PROJECT}/quants/${SRA_ACC}_quants/quant.sf`. The default values 
for the salmon_idx, out_directory, and number of threads used in `03_download_and_process.py` script can be changed. Use
`python 03_download_and_process.py --help` for more information.

In some cases, `fasterq-dump` will fail and the script will silently move onto the
next project. **Be sure to check that the number of samples processed (i.e. {SRA_acc}_quants/) is equal to the number of samples listed in `search-results.tsv for that project**

## Step 4. Create SummarizedExperiment objects for each PROJECT.

This step will run an Rscript that combines the count data for each sample in 
a PROJECT along with the sample annotations in order to create a final 
`SummarizedExperiment` object to be used in downstream DE analyses. 

To run the script on all projects:

```bash
for PROJECT in $(find . -name "PRJNA*" -type d); do
    Rscript 04_create_SE.R -d $PROJECT/quants/ -f $PROJECT/annotation.txt -o $PROJECT -c 8;
done
```

where:

- `-d` is the directory containing all `${SRA_accession}_quants/` directories for a project.
- `-f` is the sample annotation file with `sample_name` and `group` columns.
- `-o` is the output directory
- `-c` is the number of cores to use in processing

Use `Rscript 04_create_SE.R --help` for more information.

Running this script *should* produce a file called `summarized_experiment.rds` located
in the out_directory. 

## Step 4b. Create MultiAssayExperiment objects

The above Rscript that returns a `SummarizedExperiment` object pre-filters some
of the data and returns all counts as one single count matrix with colData. Instead
of returning this pre-filtered structure the `04b_create_MAE.R` script will return
a `MultiAssayExperiment` object with one assay for gene counts (all features, 
not just protein coding), RE counts (intergenic and intronic counts summed), 
exonic RE counts, intronic RE counts, and intergenic RE counts. The commandline
arguments are the same as above.
