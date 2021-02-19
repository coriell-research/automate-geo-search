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
    "$OUT_DIRECTORY
```

## Step 2. Get the run information for each project

The next step parses the search-results.tsv file and pulls the standard run-info.csv file from SRA. The script will
output each run-info.csv file in a new, separate directory named with the BIOPROJECT ID for each record present in the 
search results. **To protect against large BIOPROJECTS being imported, the script limits processing to BIOPROJECTS with fewer than (40) samples** 
This can be changed by modifying the variables in the `02_get_run_info.sh` script.

```bash
./02_get_run_info.sh path/to/search-results.tsv $OUT_DIRECTORY
```

## Step 3. Download fastqs and map with Salmon

The last step downloads all SRA accessions for each BIOPROJECT, maps each sample, then removes each fastq file after successful
completion. 

To run downloading and mapping on all BIOPROJECTS downloaded:

```bash
for PROJECT in $(find . -name "PRJNA*" -type d); do
    python 03_download_and_process.py $PROJECT/run-info.csv;
done
```

By default, the script will download the fastq files from SRA into the PROJECT folder and then perform quantification
with Salmon, placing the individual quant results in `${PROJECT}/quants/${SRA_ACC}_quants/quant.sf`. The default values 
for the salmon_idx, out_directory, and number of threads used in `03_download_and_process.py` script be changed. 