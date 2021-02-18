# E. Directory structure

> This segment provides an overview of the structure of the directories in the main scripts `plants-pipeline` directory (Section E1) and the data directory `pipeline-data` (Section E2).

***

## E1. Main scripts directory (`plants-pipeline`)

```
plants-pipeline
- config
  - setup_lin.sh *
  - setup_mac.py *
  - setup_data.py
  - query_data_path.py
  - constants.py *
  - requirements.txt
- main
  - get_species_list.py
  - get_support_files.py
  - run_job.py
  - update_all_tpm_matrices.py
- preprocess
  - po_parser
  - despatch.py
  - ena.py
  - sra.py
  - match.py (deprecated)
  - match2.py
  - sra_species.py
  - iohelper.py
- download
  - despatch.py
  - download_functions.py
  - bash_download_template.py
  - min_bash_download_tempalte.py
  - extract_runinfo.py
  - checkfiles.py
  - helpers.py
  - validate.py (deprecated)
- postprocess
  - pull_tpm_matrix.py
  - get_qc_stats.py
  - calc_pcc_matrix.py (deprecated)
  - coexpression1.py
  - calc_f1_scores.py
```

The files marked with an `*` in the directory tree above can be edited by the user. The other files should be left as it is for the pipeline to run.

***

### plants-pipeline/main

This directory contain all the scripts called for running a bulk batch of download. They generally rely on the list of species given in the latest file of `pipeline-data/preprocess/job-list`, filering our species without any CDS link provided.

This segment describes the purpose of each file. For details on how to implement the scripts, refer to [**C & D. Download guide**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/download.md) for how to specify the arguments.

### plants-pipeline/main/get_species_list.py

This script obtains the list of species belonging to the specified clade and the corresponding species taxanomic IDs.

### plants-pipeline/main/get_support_files.py

This script:
- ensures that runtable for each species is downloaded
- ensures CDS for the species is downloaded
- each downloaded CDS is processed into kallisto index
- remove CDS files

These are performed in preparation for the download.

### plants-pipeline/main/run_job.py

This is the main script which runs the download jobs for species. For each species,
- Run IDs are downloaded in parallel processes
- Multiple threads are used for kallisto pseudoalignment
- 3 attempts are made to download previously failed run IDs.

The status of the download can be viewed from `pipeline-data/download/logs/status`.

### plants-pipeline/main/update_all_tpm_matrices.py

This script helps to compile the tpm values from all stored files for each Run ID, and return tpm matrices per species in `pipeline-data/postprocess/tpm-matrices`.

***

### plants-pipeline/preprocess

This directory contains the scripts needed to obtain files before the download commence. The files obtained are the species list and runtables. The annotation scripts are also found here.

This segment describes the purpose of each file. For details on how to implement the scripts, refer to [**C & D. Download guide**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/download.md) for how to specify the arguments.

### plants-pipeline/preprocess/po_parser

This directory contains the Ontology object that is used to parsed the .obo format Plant Ontology file. An Ontology instance for Plant Ontology was saved as the variable `po` in the file `po.py` here, which can be imported to other scripts in the pipeline. If other ontologies are to be used, the `po.py` script can be edited to include load an Ontology instance with other .obo files that is not Plant Ontology.

### plants-pipeline/preprocess/despatch.py

The script to download runtables for one species, with options to chose either or both of SRA and ENA, and whether to annotate or not.

### plants-pipeline/preprocess/sra.py

The script to download runtable from SRA for a specific species

### plants-pipeline/preprocess/ena.py

The script to download runtable (metadata) from ENA for a specific species. Note that the metadata provided from ENA API is different from the runtables provided by NCBI SRA. In particular, the number of columns here are fewer and has less information than SRA runtables.

### plants-pipeline/preprocess/match2.py

The script to annotate the po term for each run ID provided in the runtables of each species. This is a second edition of the script, therefore having the number 2.

### plants-pipeline/preprocess/match.py

This is our legacy annotation script at `plants-pipeline/preprocess/match.py`. This method is deprecated because it requires a deterministic column names of the runtables while SRA provide inconsistent naming in the column headers, which makes this method less robust.

However, in future, if SRA standardise the column headers of all runtables, then this method might become useable. Therefore we decide to leave the script and components here in the pipeline.

### plants-pipeline/preprocess/iohelper.py

Helper functions are defined here and can be used across the pipeline if imported into the local script.

***

### plants-pipeline/download/

This directory contains the scripts needed to perform download, kallisto pseudoalignment and file organisation.

This segment describes the purpose of each file. For details on how to implement the scripts, refer to [**C & D. Download guide**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/download.md) for how to specify the arguments.

### plants-pipeline/download/despatch.py

This script despatch a job to download one species.

### plants-pipeline/download/download_functions.py

This script contains all the functions needed to run the download in parallel or linearly.

### plants-pipeline/download/bash_download_template.py
### plants-pipeline/download/min_bash_download_template.py

These are bash template scripts used to run the download for each experiment/sample. The variables for each Run ID will be modified and the scripts compiled in a file output stored in `pipeline-data/download/bash-jobfiles`.

### plants-pipeline/download/extract_runinfo.py

This script is used to extract runinfo files from kallisto output and delete the individual runinfo files for each sample/run ID, compiling the runinfo into a combined file in `pipeline-data/download/runinfo-main` and `pipeline-data/download/logs/runinfo`.

### plants-pipeline/download/checkfiles.py

This script helps check if run IDs are successfully downloaded.

### plants-pipeline/download/helpers.py

Helper functions are defined here and can be used across the pipeline if imported into the local script.

***

### plants-pipeline/postprocess/pull_tpm_matrix.py

Script to create a tpm matrix for one specified species, based on downloaded run IDs for that species.

### plants-pipeline/postprocess/get_qc_stats.py

This script generates quality control (QC) statistics with recommended thresholds. Using the generated visualisations and summary statistics, the user can modify the tresholds and run this script again to define the new thresholds.

### plants-pipeline/postprocess/calc_pcc_matix.py

Deprecated as the file to be stored will be too big. It will be more efficent to calculate the PCC values on the fly and store only the final output of interest.

### plants-pipeline/postprocess/coexpression1.py

Calculates the percentage ribosomal gene neighbours (by coexpression) each gene have. This prepares the data needed to summarise the F1 scores for each combination of cutoffs in the next step.

### plants-pipeline/postprocess/calc_f1_scores.py

Summarise the F1 scores for each combination of PCC cutoff and percentage ribosomal gene cutoff.

***

***

## E2. Data directory (`pipeline-data`)

The directory tree is as such:

```
pipeline-data
- preprocess
  - ena-annotations
  - ena-runtables
  - job-list
  - species-list
  - sra-annotation2s
  - sra-annotations
  - sra-runtables
- download
  - bash-jobfiles
  - bash-tmp
  - cds
  - fastq-tmp
  - idx
  - kallisto-out
  - logs
    - initiation
    - progress
    - runinfo
    - runtime
    - status
  - runinfo-main
- postprocess
  - f1-stats
  - gene-classifications
  - pcc-matrices
  - percentage-matrices
  - qc-jointplots
  - qc-matrices
  - qc-out
  - qc-summary
  - tpm-matrices
```

***

### pipeline-data/preprocess/sra-runtables

This is where the runtables (metadata) for each species, downloaded from NCBI SRA RunSelector, is stored at. The runtables are in the unmodified comma-seperated values as provided by SRA. The file names are in the format of YYYYMMDD-HHMMSS-taxidXXXX_sra_runtables.txt. It is important to maintain this convention and not modify it as the scripts look for the latest runtable for each queried species using this format.

Older runtables can be manually deleted, or it can be left in the directory for archive purposes.

### pipeline-data/preprocess/sra-annotation2s

This is a tab-delimited file with the first column being the Run ID of the experiment and the second columns being the PO term annotation assigned by LSTrAP-Kingdom.

If the field in the second column is blank, it means that LSTrAP-Kingdom cannot find the information from the SRA runtable for this species.

The annotation is made by the script in `plants-pipeline/preprocess/match2.py`. This is a second edition of the script, therefore having the number 2.

### pipeline-data/preprocess/sra-annotations

Output from our legacy annotation script at `plants-pipeline/preprocess/match.py`. This method is deprecated because it requires a deterministic column names of the runtables while SRA provide inconsistent naming in the column headers, which makes this method less robust.

However, in future, if SRA standardise the column headers of all runtables, then this method might become useable. Therefore we decide to leave the script and components here in the pipeline.

### pipeline-data/preprocess/ena-runtables

This is where the runtables (metadata) for each species, downloaded from ENA API, is stored at. The runtables are in the unmodified tab-seperated values as provided by ENA. The file names are in the format of YYYYMMDD-HHMMSS-taxidXXXX_sra_runtables.txt. It is important to maintain this convention and not modify it as the scripts look for the latest runtable for each queried species using this format.

Older runtables can be manually deleted, or it can be left in the directory for archive purposes.

### pipeline-data/preprocess/ena-annotations

Similar to `sra-annotations`.

### pipeline-data/preprocess/species-list

When calling `plants-pipeline/main/get_species_list.py`, the generated list of species within the specified clade will be stored here. The list here is to be left as it is.

### pipeline-data/preprocess/job-list

Same as species-list but this list is to be editted to fill in the cds links for each species, if available, before running the remaining steps in the pipeline.

***

### pipeline-data/download/bash-jobfiles

A compilation of all the bash scripts needed to download the desired run IDs. The script for each run ID is on one line and this file is then used to run each line of bash commands in parallel.

The template is found in `plants-pipeline/download/bash_download_template.py`.

### pipeline-data/download/bash-tmp

Deprecated.

### pipeline-data/download/cds

This is a directory to temporarily store downloaded CDS files. Once the kallisto index is built, the associated CDS file will be deleted.

### pipeline-data/download/fastq-tmp

This is used to temporarily store downloaded fastq files, which will then be used for pseudoalignement and quantification by kallisto. The fastq file is then deleted once the kallisto output is produced.

### pipeline-data/download/idx

Stores kallisto index generated from the CDS of each species.

### pipeline-data/download/kallisto-out

Stores kallisto output. Mirrors the directory structure of ENA. For example, SRR123456 will be stored as `SRR123/SRR123456.zip`. SRR3334567 will be stores as `SRR333/007/SRR3334567.zip`.

### pipeline-data/download/logs

- initiation: logs the time of the run IDs whose downloaded were initiated/attempted.
- progress: logs whether each run ID has been successfully downloaded. This is used to keep track and redownload only undownloaded run IDs.
- runtime: Logs the time taken to download and process each run ID.
- runinfo: Logs the compiled runinfo from each download attempt.
- status: Logs the success rates for each attempt of each species and the total time taken for each attempt.

### pipeline-data/download/runinfo-main

The most updated runinfo compiled for each species. Will be automatically generated by combining the runinfo files for the species from `pipeline-data/download/logs`.

***

### pipeline-data/postprocess/tpm-matrices

All the tpm values for the genes in all experiments for a species are compiled in TPM matrcies, which are stored here.

### pipeline-data/postprocess/qc-matrices

Stores the matrices of the number of samples passing QC, given a set of log processed or % pseudoaligned thresholds.

Sample for a species:

```
log_processed|p_pseudoaligned 10  20  30  40  50  60  70  80  90
6.0 227 212 210 207 152 26  0 0 0
6.2 227 212 210 207 152 26  0 0 0
6.4 223 208 206 203 152 26  0 0 0
6.6000000000000005  221 206 204 201 152 26  0 0 0
6.800000000000001 216 201 199 197 148 26  0 0 0
7.000000000000001 175 162 160 159 144 26  0 0 0
7.200000000000001 74  73  73  73  67  1 0 0 0
7.400000000000001 9 9 9 9 8 0 0 0 0
7.600000000000001 4 4 4 4 4 0 0 0 0
7.800000000000002 0 0 0 0 0 0 0 0 0
8.000000000000002 0 0 0 0 0 0 0 0 0
```

### pipeline-data/postprocess/qc-summary

A summary of the recommended thresholds to set and the number of experiments passing the QC.

### pipeline-data/postprocess/qc-jointplots

Stores the jointplots of the samples for visualisation on how to set thresholds.

### pipeline-data/postprocess/qc-out

Stores a json file of the run IDs passing the QC, for downstream use.

### pipeline-data/postprocess/gene-classifications

For annotation benchmarking purposes. This is a tsv file obtained from [Mercator 4.2](https://www.plabipd.de/portal/mercator4), using the same CDS used in this pipeline as the input. The file for each species here define the Mapman Bin assigned to each gene in the CDS.

### pipeline-data/postprocess/f1-stats

For F1 scores benchmarking purposes

### pipeline-data/postprocess/percentage-matrices

For F1 scores benchmarking purposes

### pipeline-data/postprocess/pcc-matrices

Deprecated as the file to be stored will be too big. It will be more efficent to calculate the PCC values on the fly and store only the final output of interest.

***
