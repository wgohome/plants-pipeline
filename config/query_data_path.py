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
import argparse
from config.constants import DATA_PATH, ASPERA_SSH_KEY

print(f"Current DATA_PATH is:\n{DATA_PATH}")
if not os.path.exists(DATA_PATH):
    print("This path is invalid.")
else:
    print()
print(f"Current ASPERA_SSH_KEY path is:\n{ASPERA_SSH_KEY}\n")
