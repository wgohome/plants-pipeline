### First setup

Ensure you are in the plants-pipeline directory.

To setup for the first time, create a virtual environment. Ensure you have installed Python 3.8. Run the following command. If it is aliased as python3, replace `python3.8` with `python3`. `virtualenv` can be replaced with any other Python environment management tools you prefer.
```
virtualenv -p python3.8 proj_env
pip install --upgrade pip
pip install -r requirements.txt
```

To setup directories tree for this project's data repository, run this command, replacing `path/to/data/repository/` with your desired path.
```
python3 ./cofig/setup_data.py -p path/to/data/repository/
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
