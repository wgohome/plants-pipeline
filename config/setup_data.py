# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
    PARENT_MODULE = re.search('^(.*plants-pipeline)', realpath(dirname(__file__))).group()
################################################################################

import os
import argparse
import pdb

parser = argparse.ArgumentParser(description = 'This script sets up the directory structure of the data repository for the plants-pipeline scripts.', epilog = 'By Mutwil Lab')
parser.add_argument('-p', '--path', nargs=1, metavar='data_path',
                    help='Enter the path of the directory where the data directories for this project is to be created',
                    dest='data_path', type=str, required=True)
args = parser.parse_args()
DATA_PATH = args.data_path[0] + '/pipeline-data/'
assert os.path.exists(f"{args.data_path[0]}") == True, "Target directory invalid."

# Default directories, minimally needed to run the download pipeline. Additional files for post-processing can be added
dests = {
    'download': [
        'cds',
        'idx',
        'bash-jobfiles',
        'bash-tmp',
        'fastq-tmp',
        'kallisto-out',
        'runinfo-main',
        'logs/initiation',
        'logs/runtime',
        'logs/progress',
        'logs/runinfo',
        'logs/status'
    ],
    'preprocess': [
        'ena-runtables',
        'sra-runtables',
        'ena-annotations',
        'sra-annotations',
        'species-list',
        'job-list'
    ],
    'postprocess': [
        ''
    ]
}

def makedirs(dest):
    os.makedirs(f"{DATA_PATH}{dest}", exist_ok=True)

if __name__ == "__main__":
    if os.path.exists(DATA_PATH):
        print("WARNING: This DATA_PATH already exists. But subdirectories that were not existing will be created.")
    for p, c in dests.items():
        for dest in c:
            if not os.path.exists(f"{DATA_PATH}/{p}/{dest}"):
                makedirs(f"{p}/{dest}")
                print(f"Created: {p}/{dest}")
    print("🚀 Created directories for pipeline data repository")
    # Update constants file
    f = open(f"{PARENT_MODULE}/config/constants.py").readlines()
    to_write = []
    for line in f:
        if "DATA_PATH =" in line:
            to_write.append(f"DATA_PATH = \"{DATA_PATH}\"\n")
        else:
            to_write.append(line);
    with open(f"{PARENT_MODULE}/config/constants.py", 'w') as f:
        f.writelines(to_write)
    print(f"./config/constants.py has been edited with:\n\tDATA_PATH = \"{DATA_PATH}\"\n")
