# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import os
# Relative imports of CONSTANTS in config/constants.py
from config.constants import *

out_path = f"{DATA_PATH}/download/fastq_tmp/{'/'.join(unpaired.split('/')[:-1])}"

# Pull out latest log for species and find out completed/incomplete Run IDs
# Check if kallisto output present
# If kallisto output present, move RunID kallisto file to kallisto-tmp, arranged in ENA dir structure
# If kallisto output present, delete fastq from fastq-tmp
# Create log to update completed RunIDs and incomplete RunIDs

def check_kallisto(runid):
    kal_path = f"{DATA_PATH}/kallisto_tmp/{runid}/run_info.json"
