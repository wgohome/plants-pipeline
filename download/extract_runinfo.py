# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import json
import argparse
# relative imports
from config.constants import DATA_PATH
import helpers

parser = argparse.ArgumentParser(description = 'This script is called from the shell scipt to read the runinfo string and compile it into a runinfo table.', epilog = 'By Mutwil Lab')
parser.add_argument('-r', '--runinfo-log', metavar='runinfo_log',
                    help='Enter the species runifo log to append to. For instance, \'[data-path]/download/logs/runinfo/20200730-220000-taxid3702_runinfo.txt\' for Arabidopsis thaliana.',
                    dest='runinfo_log', type=str, required=True)
parser.add_argument('-k', '--kallisto-path', metavar='kal_out',
                    help='Enter the full kallisto output path. For instance, \'[data-path]/download/kallisto-out/SRR123/007/SRR1234567\'.',
                    dest='kal_out', type=str, required=True)

args = parser.parse_args()
runinfo_log = args.runinfo_log
runinfo_path = args.kal_out
runinfo_string = open(runinfo_path, 'r').read()
runinfo = json.loads(runinfo_string)
values = [str(val) for val in runinfo.values()]
to_write = '\t'.join(values) + '\n'
with open(runinfo_log, 'a') as f:
    f.write(to_write)
