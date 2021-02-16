# LSTrAP-Kingdom: an automated pipeline to generate annotated gene expression atlases for kingdoms of life

By [**NTU Plants Systems Biology and Evolution Laboratory**](https://www.plant.tools)

This repository is found in [this Github Repository](https://github.com/wirriamm/plants-pipeline). Do create pull requests for issues/bugs and feature requests.

[Contact me](mailto:will0046@e.ntu.edu.sg) for feedback or reporting bugs.

[A. First Setup]()

## A. First setup

Clone this repository to your machine. Ensure the working directory is in the main plants-pipeline directory.

### A1. Setting up Python environment

To setup this pipeline for the first time, create a virtual environment. Ensure you have installed 'Python 3.8' and 'pip 20.1.1'. Run the following command below. If your local 'Python 3.8' installation is aliased as python3, replace `python3.8` in the command below with `python3`. `virtualenv` can be replaced with any other Python environment management tools you prefer like `venv` which comes with the standard Python library.
```
virtualenv -p python3.8 proj_env
pip install --upgrade pip
pip install -r config/requirements.txt
```

### A2. Setting data path

To setup the directories for this project's data repository, run this command, replacing `path/to/data/repository/` with your desired path. The pipeline's main data repository named `pipeline-data` will be stored within this directory you specified. Note that `path/to/data/repository/` should be separate from this repository's `plants-pipeline` directory and should not be tracked by this repository's version control.
```
python3 config/setup_data.py -p /path/to/data/repository/
```
A data directory will be created at `/path/to/data/repository/pipeline-data` if it does not already exists. If it exists, sub-directories that are missing will be created. The variable for `DATA_PATH` will then be updated as `/path/to/data/repository/pipeline-data` in `config/constants.py` file.

Then, you can open the `config/constants.py` file. Check that `DATA_PATH` is correct.

### A3. Aspera download

To download Aspera ascp program and kallisto pseudoaligner, run our setup script:

For MacOS,
```
source config/setup_mac.sh
```
For Linux,
```
source config/setup_lin.sh
```

In `config/constants.py` file, also edit the `ASPERA_SSH_KEY` variable accordingly based on your local machine installation of Aspera ascp.

For linux machines, it is usually at:
`ASPERA_SSH_KEY = "/home/user/.aspera/cli/etc/asperaweb_id_dsa.openssh"`

For Macintosh, it is usally at:
`ASPERA_SSH_KEY = "/Users/user/Applications/Aspera CLI/etc/asperaweb_id_dsa.openssh"`

### A4. Other dependencies for web scraping

Then, ensure that your machine has the following programs installed
- chromium
- chromedriver (compatible with your version of chromium)

### A5. Dotenv file

Create a `.env` file in the `plants-pipeline` directory, speciying your own NCBI API key as:
```
NCBI_API_KEY = "{enter_the_key}"
```
A personal NCBI API Key can be obtained by signing up for an NCBI account. Guides can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK25497/#top)

## B. For each subsequent run of the pipeline

### Set up environment and dependencies

Begin by entering the main directory of this pipeline, which is `plants-pipeline` if you cloned from this Github repository.

```
cd /path/to/plants-pipeline
```

Run these commands to set up the environment for each new session.

To activate this python environment and all its packages:
```
source proj_env/bin/activate
```

Then, to setup kallisto, ascp commands and other dependencies,

For MacOS,
```
source config/setup_mac.sh
```
For Linux,
```
source config/setup_lin.sh
```

## C. Bulk Run

### C1. Get species list
```
python main/get_species_list.py -n viridiplantae -t 33090
```
- `-n` is for the taxanomic name of the family or any other ranks above species.
- `-t` is for the taxanomic id of that name.

The species list will be found in `pipeline-data/preprocess/species-list/` directory.

The same list will also be foind in `pipeline-data/preprocess/job-list/` directory. This tsv file is the one to be editted by the user in an Excel or equivalent interface. For the species to be downloaded, user can input the ftp/http url to the cds file for the species under the `cds_link` column. For species that the user does not want to download (number of experiments too low), the user can delete the row from this file or just leave the field under `cds_link` blank. Only species rows with a filled CDS link will be processed by the pipeline.

### C2. Download cds and runtables, create kallisto index
```
python main/get_support_files.py
```
This script will run based on the latest timestamped file in `pipeline-data/preprocess/job-list/`. (It is the file in which you have added the urls under `cds_link` column.)

Note that occassionally, NCBI website does face a lagging server, therefore web scraping may fail. In that case, rerun this segment, or try again another day.

### C3. Running the download job
```
nohup python main/run_job.py -m ascp-bash -w 20 -t 2 &
```
The compulsory arguments are:
- `-m` is for one of the three download methods: 'ascp-bash', 'ascp-python' or 'curl'. `ascp-bash` is the desired method of this pipeline.
- `-w` is to set the number of workers for multiprocessing.
- `-t` is to set the number of threads when performing kallisto quant.

`nohup` is used to run the job in the background. Any other equivalent programs may be used.

### C4. Checking the download status

Download status can be checked from `pipeline-data/download/logs/status/YYYYMMDD-HHMMSS_job.log`. This file will be updated as the download proceeds.

### C5. Reinitiating download

If upon checking the status logs for a batch of mass download, you notice that many species have a high number of failed downloads, you can re-run the command of section C3. This will download only the experiments/Run IDs that have previously failed.

### C6. Extracting run info tables and tpm matrices

After bulk download and/or re-initiating bulk downloads to a satisfactory coverage of experiments/Run IDs, a TPM matrices per species can be extracted from the kallisto output as such:
```
python main/update_all_tpm_matrices.py
```

The TPM matrix, with file name timestamped and labelled by taxid, can be found in `pipeline-data/postprocess/tpm-matrices`.

## D. Small jobs: Running for just one species

### D1. Get runtables

Runtables can be obtained from NCBI's SRA or ENA. The runtables will be saved in `pipeline-data/preprocess/sra-runtables` and `pipeline-data/preprocess/ena-runtables`. Files are labelled by their taxanomic id.

To download the runtables, call the following command, replacing 'Arabidopsis thaliana' with the species name of interest. Take note to wrap the species name in quotation marks.

```
python preprocess/main.py -s 'Arabidopsis thaliana' -d sra
```
or
```
python preprocess/main.py -s 'Arabidopsis thaliana' -d ena
```

The available options are:
- `-s` for specifying the full species name. Mandatory.
- `-d` for specifying the database source, accepting either `sra` or `ena`. Mandatory.
- `-b` is an optional tag, which can be specified without any arguments, to extract organ annotations from the runtables.

Organ annotations (if activated) will be stored in `data-pipeline/preprocess/sra-annotations` and `data-pipeline/preprocess/ena-annotations` respectively.


### D2. Download one species

Call the `despatch.py` script with the following arguments.
- `-s` is for the 3 letter alias for the species name. For example, Arabidopsis thaliana should have the alias 'Ath'.
- `-c` is for the name (not full path) of the CDS fasta file found in `pipeline-data/download/cds/`. It will be good to set up a convention such as 'Ath.cds.fasta' for the file naming.
- `-m` is for one of the three download methods: 'ascp-bash', 'ascp-python' or 'curl'.
- `-l` is an optional tag to indicate if download is to be done linearly. By default, download will be in parallel processes. This only applies if --method chosen was 'ascp-python' or 'curl'
- `-w` is an optional tag to set the number of workers for miltiprocessing. If download is to be done linearly, this argument will be ignored. By default, number of workers is set to 8.
- `-t` is an optional tag to set the number of threads when performing kallisto quant.

For example, some possible commands are:

**Main use case:**
```
python download/despatch.py -s Ath -c Ath.cds.fasta -m ascp-bash -w 20
```

**Other use cases:**
```
python download/despatch.py -s Ath -c Ath.cds.fasta -m ascp-python
python download/despatch.py -s Ath -c Ath.cds.fasta -m ascp-python -w 4
python download/despatch.py -s Ath -c Ath.cds.fasta -m curl -l
```

To run the download in the background, instead, run this:
```
nohup python download/despatch.py -s Ath -c Ath.cds.fasta -m ascp-bash -w 20 &
```
The stdout will be found in the file 'nohup.out'. However, due to multi-processing, the stdout might not make much sense. Instead, the `pipeline-data/download/log` files can be more informative.

### D3. Checking the logs

Once completed, the time for download could be checked from log file in pipeline-data/download/logs/runtime, labelled by the time the download was initiated.

Once the download have completed, two logfiles can be checked. Logfiles' names begin with a timestamp in the format 'YYYYMMDD-HHMMSS', indicating the time at which the download batch was initiated.

`pipeline-data/download/logs/` contain several directories:
- `initiation` directory stores the files containing the timestamp at which each Run ID in the batch was initiated.
- `runtime` directory stores the files containing the timestamp, library layout and download error type (if any) of each Run ID. However, if the child process for downloading the RunID was broken, it will not be logged in this file. The Run ID should still be found in the logfile from `initiation`.
- `progress` directory stores files indicating whether or not each Run ID for the species has been successfully downloaded.
- `status` directory stores the logs for the download status of a bulk download run.
- `runinfo` directory stores the history of runinfo from different batches of downloads. However for the compiled runinfo table per species, refer instead to `pipeline-data/download/runinfo-main` directory.

_To be updated_

### D4. Checking the kallisto output files and reinitiating new batch for missing files

After downloading a batch, to check the successful downloads against the Run IDs in the runtable, run the following file.
- `-s` specifies the species three-letter alias

```
python download/checkfiles.py -t 3702
```

This script will:
- Check `pipline-data/download/kallisto-out` for successful downloads of kallisto output files.
- Update a new log file in `pipline-data/download/logs/progress/[timestamp]-[spe]-progress.log`.
- Update `pipline-data/download/runinfo-main` log file for the species' runinfo table.

### D5. Redownloading failed Run IDs

Step D2 can be repeated.
