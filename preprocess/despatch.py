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
from preprocess import match, sra, ena, iohelper, match2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script downloads the runtable from SRA or ENA to pipeline-data/preprocess directory.', epilog = 'By Mutwil Lab')
    parser.add_argument('-s', '--species-name', metavar='species_name',
                        help='Enter the name of the species to be downloaded, in quotes. For example, "Arabidopsis thaliana". If taxid is provided, species name is not neccessary.',
                        dest='species', type=str, required=False)
    parser.add_argument('-t', '--taxid', metavar='taxid',
                        help='Enter the taxanomic id of the species to be downloaded, in quotes. For example, "3702" for Arabidopsis thaliana. If taxid is not known, providing the species scientific name would suffice',
                        dest='taxid', type=int, required=False)
    parser.add_argument('-d', '--database', metavar='database',
                        help='The choice of database source for the runtable, either \'sra\' (from NCBI) or \'ena\' from ENA.',
                        dest='database', type=str, required=True)
    parser.add_argument('-a', '--annotate', action='store_true', default=False, required=False, dest='annotate', help="Include this optional tag if the runtable is supposed to be annotated with organ labels.")
    parser.add_argument('-b', '--annotate2', action='store_true', default=False, required=False, dest='annotate2', help="Include this optional tag if the runtable is supposed to be annotated with organ labels by Method 2.")
    args = parser.parse_args()
    species = args.species
    taxid = args.taxid
    db = args.database.lower()
    annotate = args.annotate
    annotate2 = args.annotate2
    if (taxid == None) and (species == None):
        raise Exception("Minimum requirement: enter either species or taxid argument.")

def fetch_runtable(species=None, taxid=None, db='sra', annotate=False, annotate2=False):
    if (taxid == None) and (species == None):
        raise Exception("‼️ Minimum requirement: enter either species or taxid argument.")
    elif taxid == None:
        taxid = sra.get_taxanomic_id(species)
    elif species == None:
        species = iohelper.species_name(taxid)
    spe_id = f"taxid{taxid}"
    if db == 'sra':
        # sra_runtable_path = iohelper.make_filename(db='sra', filetype='runtable', spe_id=spe_id)
        sra_df = sra.process_sra_runtable(taxid, species) # Downloads SRA runtable, returns df
        if annotate:
            match.annotate(sra_df, species, db='sra') # Saves annotation file
        if annotate2:
            match2.annotate2(sra_df, species, db='sra')
    elif db == 'ena':
        ena_df = ena.get_ena_runtable(taxid)
        ena_runtable_path = iohelper.make_filename(db='ena', filetype='runtable', spe_id=spe_id)
        ena_df.to_csv(ena_runtable_path, sep='\t', index=False)
        if annotate:
            match.annotate(ena_df, species, db='ena')
    print(f"Downloaded {db} runtable for {species} ({spe_id}) 📥")

if __name__ == '__main__':
    fetch_runtable(taxid=taxid, species=species, db=db, annotate=annotate)

__all__ = ['fetch_runtable']
