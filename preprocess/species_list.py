# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################
"""
Script to get an updated species list file with:
- taxid
- species name
- number of sra experiments
"""

from preprocess.sra_species import make_species_report

master_df = make_species_report()
