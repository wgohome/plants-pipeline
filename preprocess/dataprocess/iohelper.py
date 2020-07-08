# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

def create_annotation_file(species, study, db):
    db = db.lower()
    Spe = species_shortform(species)
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*pipeline)', abspath).group()
    out_path = f"{dirname(parent_module)}/data/{study}/out/{db}_annotations/{Spe}_{db}_annotation.txt"
    with open(out_path, 'w') as f:
        f.write("")
    return out_path

def write_annotation(line, out_path):
    with open(out_path, 'a') as f:
        f.write(line)

def get_species_list(path):
    with open(path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    return lines

def species_shortform(species):
    # Arabidopsis thaliana => Ath
    return species.split()[0][0] + species.split()[1][0:2]

__all__ = ['write_annotation', 'get_species_list', 'species_shortform']
