# B. For each subsequent run of the pipeline

> These commands need to be run everytime the pipeline is accessed from a new terminal session. They will load the python environment with the installed packages, and add ascp and kallisto commands to the global environment $PATH. If kallisto or ascp(Aspera CLI) is not downloaded, they will also be downloaded.

## Set up environment and dependencies

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
