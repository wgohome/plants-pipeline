# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

from config.constants import DATA_PATH
from preprocess import match, sra

spe = 'Egu'
species = 'Elaeis guineensis'
sra_runtable_path = f"{DATA_PATH}/preprocess/sra-runtables/{spe}_sra_runtable.txt"
spe_df = sra.read_sra_runtable(sra_runtable_path)
match.annotate(spe_df, species, db='sra')
print("Done 😄")
