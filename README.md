# Plants Transcriptomics Pipeline

---

**NTU Plants Systems Biology and Evolution Laboratory**

This repository is found in [Github Repository](https://github.com/wirriamm/plants-pipeline). Create pull requests for issues/bugs.

[Contact me](mailto:will0046@e.ntu.edu.sg)

---

### First setup

Clone the repository to you rmachine. Ensure you are in the plants-pipeline directory.

To setup for the first time, create a virtual environment. Ensure you have installed Python 3.8. Run the following command below. If it is aliased as python3, replace `python3.8` with `python3`. `virtualenv` can be replaced with any other Python environment management tools you prefer.
```
virtualenv -p python3.8 proj_env
pip install --upgrade pip
pip install -r requirements.txt
```

To setup directories tree for this project's data repository, run this command, replacing `path/to/data/repository/` with your desired path.
```
python3 ./cofig/setup_data.py -p /path/to/data/repository/
```

Then, open `config/constants.py` file and edit the `ASPERA_SSH_KEY` and `DATA_PATH` variables accordingly.

```
DATA_PATH = "path/to/data/repository/"
ASPERA_SSH_KEY = "/Users/wirriamm/Applications/Aspera CLI/etc/asperaweb_id_dsa.openssh" # For MacOS
ASPERA_SSH_KEY = "/home/.aspera/cli/Aspera CLI/etc/asperaweb_id_dsa.openssh" # For Linux
```

### For subsequent runs

Run these commands to set up the environment for each use.

To activate this python environment and all its packages:
```
source proj_env/bin/activate
```

Then, to setup kallisto and ascp commands,

For MacOS,
```
source setup_mac.sh
```
For Linux,
```
source setup_lin.sh
```

### Download one species

Call the `despatch.py` script with the following arguments. `-i` is for the 3 letter alias for the species name and `-c` is for the name of the CDS fasta file found in pipeline_data/download/cds/ .

```
python download/despatch.py -i Ath -c Ath.cds.fasta
```

Once completed, the time for download could be checked from log file in /pipeline-data/download/logs/time, labelled by the timestamp the download was initiated.
