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
import pdb
# Relative imports
from preprocess import sra_species

parser = argparse.ArgumentParser(description = 'This script gets a list of species with their corresponding number of RNA-seq samples available.', epilog = 'By Mutwil Lab')
parser.add_argument('-n', '--name', nargs=1, metavar='name',
                    help="Enter the name of the family, eg: 'Viridiplantae'.",
                    dest='name', type=str, required=True)
parser.add_argument('-t', '--taxid', nargs=1, metavar='taxid', required=True, dest='taxid', help="Enter the taxid of the given family name.")
args = parser.parse_args()
name = args.name[0].lower()
taxid = int(args.taxid[0])

master_df = sra_species.make_species_report(name, taxid)
