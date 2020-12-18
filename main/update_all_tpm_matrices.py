# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import argparse
import pandas as pd
import os
# Relative imports
from config.constants import DATA_PATH
from download import helpers
from postprocess import pull_tpm_matrix
from main import run_job

taxids = run_job.get_valid_jobs()
for taxid in taxids:
    pull_tpm_matrix.write_tpm_matrix(taxid)
