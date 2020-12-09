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
import datetime as dt
import pandas as pd
import pdb
# Local imports
from config.constants import DATA_PATH
from preprocess.ena import get_taxanomic_id

def get_timestamp():
    return dt.datetime.now().strftime('%Y%m%d-%H%M%S')

def make_filename(db, filetype='runtable', spe_id=None, species=None):
    db = db.lower()
    filetype = filetype.lower()
    if spe_id == None:
        spe_id = species_label(species)
    return f"{DATA_PATH}/preprocess/{db}-{filetype}s/{get_timestamp()}-{spe_id}_{db}_{filetype}.txt"

def create_annotation_file(species, db):
    db = db.lower()
    spe_id = species_label(species)
    out_path = make_filename(db, filetype='annotation', spe_id=spe_id)
    with open(out_path, 'w') as f:
        f.write("")
    return out_path

def write_annotation(line, out_path):
    with open(out_path, 'a') as f:
        f.write(line)

# def get_species_list(path):
#     with open(path, 'r') as f:
#         lines = [line.strip() for line in f.readlines()]
#     return lines

def species_shortform(species):
    '''Arabidopsis thaliana => Ath'''
    assert len(species.split() == 2), "species name must comprise of 2 words - the genus and species"
    return species.split()[0][0] + species.split()[1][0:2]

def species_label(species):
    '''3702 => taxid3702'''
    directory = f"{DATA_PATH}/preprocess/job-list"
    latest_list = max(os.listdir(directory))
    df = pd.read_csv(f"{directory}/{latest_list}", sep='\t')
    taxid = df.loc[df['species'].str.lower() == species.lower()]['taxid'].values[0]
    pdb.set_trace()
    return f"taxid{taxid}"

def species_name(taxid):
    directory = f"{DATA_PATH}/preprocess/job-list"
    latest_list = max(os.listdir(directory))
    df = pd.read_csv(f"{directory}/{latest_list}", sep='\t')
    species_name = df.loc[df['taxid'] == taxid]['species'].values[0]
    return species_name

__all__ = ['create_annotation_file', 'write_annotation', 'get_species_list', 'species_shortform']
