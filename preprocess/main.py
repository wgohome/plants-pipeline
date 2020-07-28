# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

import numpy as np
import pandas as pd
import argparse
import pdb
# Relative imports
from config.constants import DATA_PATH
from preprocess import match, sra, ena, iohelper

parser = argparse.ArgumentParser(description = 'This script downloads the runtable from SRA or ENA to pipeline-data/preprocess directory.', epilog = 'By Mutwil Lab')
parser.add_argument('-s', '--species-name', metavar='species_name',
                    help='Enter the name of the species to be downloaded, in quotes. For example, "Arabidopsis thaliana".',
                    dest='species', type=str, required=True)
parser.add_argument('-d', '--database', metavar='database',
                    help='The choice of database source for the runtable, either \'sra\' (from NCBI) or \'ena\' from ENA.',
                    dest='database', type=str, required=True)
parser.add_argument('-a', '--annotate', action='store_true', default=False, required=False, dest='annotate', help="Include this optional tag if the runtable is supposed to be annotated with organ labels.")

args = parser.parse_args()
species = args.species.capitalize()
db = args.database.lower()
annotate = args.annotate
spe_id = iohelper.species_label(species)

sra_runtable_path = f"{DATA_PATH}/preprocess/sra-runtables/{spe_id}_sra_runtable.txt"
ena_runtable_path = f"{DATA_PATH}/preprocess/ena-runtables/{spe_id}_ena_runtable.txt"
tax_id = sra.get_taxanomic_id(species)

if db == 'sra':
    sra_df = sra.process_sra_runtable(tax_id, species) # Downloads SRA runtable, returns df
    if annotate:
        match.annotate(sra_df, species, db='sra') # Saves annotation file
elif db == 'ena':
    ena_df = ena.get_ena_runtable(tax_id)
    ena_df.to_csv(ena_runtable_path, sep='\t', index=False)
    if annotate:
        match.annotate(ena_df, species, db='ena')
print(f"Downloaded {db} runtable and written annotation for {species} ({spe_id}) 😄")
