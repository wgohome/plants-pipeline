# Plants Transcriptomics Pipeline

---

**NTU Plants Systems Biology and Evolution Laboratory** (Mutwil's Lab)

This repository is found in [Github Repository](https://github.com/wirriamm/plants-pipeline). Create pull requests for issues/bugs.

[Contact me](mailto:will0046@e.ntu.edu.sg)

---

## First setup

Clone the repository to your machine. Ensure you are in the main plants-pipeline directory.

To setup this pipeline for the first time, create a virtual environment. Ensure you have installed 'Python 3.8' and 'pip 20.1.1'. Run the following command below. If your local 'Python 3.8' installation is aliased as python3, replace `python3.8` in the command below with `python3`. `virtualenv` can be replaced with any other Python environment management tools you prefer like `venv` which comes with the standard Python library.
```
virtualenv -p python3.8 proj_env
pip install --upgrade pip
pip install -r requirements.txt
```

To setup the directories for this project's data repository, run this command, replacing `path/to/data/repository/` with your desired path.
```
python3 cofig/setup_data.py -p /path/to/data/repository/
```

Then, open `config/constants.py` file and edit the `ASPERA_SSH_KEY` and `DATA_PATH` variables accordingly based on your local machine. In general, `ASPERA_SSH_KEY` is located in a standard install path for Mac and another path for Linux.

```
DATA_PATH = "path/to/data/repository/"
ASPERA_SSH_KEY = "/Users/[your-username]/Applications/Aspera CLI/etc/asperaweb_id_dsa.openssh" # For MacOS
ASPERA_SSH_KEY = "/home/.aspera/cli/Aspera CLI/etc/asperaweb_id_dsa.openssh" # For Linux
```

## For subsequent runs

Begin by entering the main directory of this pipeline, which is `plants-pipeline` if your cloned from this Github repository. Run these commands to set up the environment for each new session.

To activate this python environment and all its packages:
```
source proj_env/bin/activate
```

Then, to setup kallisto and ascp commands,

For MacOS,
```
source config/setup_mac.sh
```
For Linux,
```
source config/setup_lin.sh
```

### Download one species

Call the `despatch.py` script with the following arguments.
- `-i` is for the 3 letter alias for the species name. For example, Arabidopsis thaliana should have the alias 'Ath'.
- `-c` is for the name (not full path) of the CDS fasta file found in `pipeline-data/download/cds/`. It will be good to set up a convention such as 'Ath.cds.fasta' for the file naming.

```
python download/despatch.py -i Ath -c Ath.cds.fasta
```

To run the download in the background, instead, run this:
```
nohup python download/despatch.py -i Ath -c Ath.cds.fasta &
```
The stdout will be found in the file 'nohup.out'. However, due to multi-processing, the stdout might not make much sense. Instead, the log files can be more informative.

### Checking the logs

Once completed, the time for download could be checked from log file in pipeline-data/download/logs/time, labelled by the timestamp the download was initiated.

Once the download have completed, two logfiles can be checked. Logfiles' names begin with a timestamp in the format 'YYYYMMDD-HHMMSS', indicating the time at which the download batch was initiated.

`pipeline-data/download/logs/` contain two directories:
- `initiation` directory stores the files containing the timestamp at which each Run ID in the batch was initiated.
- `runtime` directory stores the files containing the timestamp, library layout and download error type (if any) of each Run ID. However, if the child process for downloading the RunID was broken, it will not be logged in this file. The Run ID should still be found in the logfile from `initiation`.

### Checking the kallisto output files and reinitiating new batch for missing files

_To be updated_















