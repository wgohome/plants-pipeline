### First setup
To setup for the first time, create a virtual environment. Ensure you have installed Python 3.8. Run the following command. If it is aliased as python3, replace `python3.8` with `python3`.
```
python3.8 -m venv proj_env
pip install --upgrade pip
pip install -r requirements.txt
```

### For subsequent runs
```
source proj_env/bin/activate
```
For MacOS,
```
source setup_mac.sh
```
For Linux,
```
source setup_lin.sh
```
