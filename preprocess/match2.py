# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

# Packages
import re
import argparse
import pandas as pd
import numpy as np
from nltk.stem.snowball import SnowballStemmer
import pdb
# Relative imports
from preprocess.po_parser.po import po # po is an Ontology class instance
from preprocess import iohelper
from download import helpers

stemmer = SnowballStemmer(language='english')

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description = 'This script annotates the specified sra runtable.', epilog = 'By Mutwil Lab')
#     parser.add_argument('-s', '--species-name', metavar='species_name',
#                         help='Enter the name of the species to be downloaded, in quotes. For example, "Arabidopsis thaliana". If taxid is provided, species name is not neccessary.',
#                         dest='species', type=str, required=False)
#     parser.add_argument('-t', '--taxid', metavar='taxid',
#                         help='Enter the taxanomic id of the species to be downloaded, in quotes. For example, "3702" for Arabidopsis thaliana. If taxid is not known, providing the species scientific name would suffice',
#                         dest='taxid', type=int, required=True)
#     parser.add_argument('-d', '--database', metavar='database',
#                         help='The choice of database source for the runtable, either \'sra\' (from NCBI) or \'ena\' from ENA.',
#                         dest='database', type=str, required=True)
#     args = parser.parse_args()
#     # species = args.species
#     taxid = args.taxid
#     db = args.database.lower()
#     annotate = args.annotate
#     if (taxid == None) and (species == None):
#         raise Exception("Minimum requirement: enter either species or taxid argument.")

def annotate2(df, species, db='sra'):
    out_path = iohelper.create_annotation_file(species, db, filetype='annotation2')
    po_to_stems = po.po_to_stems()
    for row in df.itertuples():
        words = []
        for col in row[1:]:
            if pd.isna(col):
                continue
            for word in col.split():
                words.append(stemmer.stem(word))
        # Search words in row against po terms
        matches = []
        for po_term, po_stems in po_to_stems.items():
            common_stems = set(words) & set(po_stems)
            if len(common_stems)>0:
                matches.append([len(common_stems)/len(po_stems), len(po_stems), po_term])
        if matches:
            matches = sorted(matches, key=lambda x: x[0:2], reverse=True)
            to_write = f"{row.Run}\t{matches[0][2]}\n"
        else:
            to_write = f"{row.Run}\t\n"
        iohelper.write_annotation(to_write, out_path)
    return None
