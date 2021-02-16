# A. First Setup of the Pipeline

> This segment only needs to be implemented at the first setup of this repository on your local machine/server.

Begin by cloning the repository to a desired directory on your local machine/server using either one or two of these commands:

By http:
```
git clone https://github.com/wirriamm/plants-pipeline.git
```

Or alternatively, by ssh:
```
git clone git@github.com:wirriamm/plants-pipeline.git
```

Ensure the working directory is in the main plants-pipeline directory.
```
cd plants-pipeline
```

## A1. Setting up Python environment

To setup this pipeline for the first time, create a virtual environment. Ensure you have installed 'Python 3.8' and 'pip 20.1.1'. Run the following command below. If your local 'Python 3.8' installation is aliased as python3, replace `python3.8` in the command below with `python3`. `virtualenv` can be replaced with any other Python environment management tools you prefer like `venv` which comes with the standard Python library.

```
virtualenv -p python3.8 proj_env
pip install --upgrade pip
pip install -r config/requirements.txt
```

## A2. Setting data path

To setup the directories for this project's data repository, run this command, replacing `path/to/data/repository/` with your desired path. The pipeline's main data repository named `pipeline-data` will be stored within this directory you specified. Note that `path/to/data/repository/` should be separate from this repository's `plants-pipeline` directory and should not be tracked by this repository's version control.
```
python3 config/setup_data.py -p /path/to/data/repository/
```
A data directory will be created at `/path/to/data/repository/pipeline-data` if it does not already exists. If it exists, sub-directories that are missing will be created. The variable for `DATA_PATH` will then be updated as `/path/to/data/repository/pipeline-data` in `config/constants.py` file.

Then, you can open the `config/constants.py` file. Check that `DATA_PATH` is correct.

To query the data path from the command line, run:

```
python config/query_data_path.py
```

## A3. Aspera download

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

## A4. Other dependencies for web scraping

Then, ensure that your machine has the following programs installed
- chromium
- chromedriver (compatible with your version of chromium)

## A5. Dotenv file

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
