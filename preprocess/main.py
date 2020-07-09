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

sra_runtable_path = f"{DATA_PATH}/preprocess/sra-runtables/Ath_sra_runtable.txt"
ath_df = sra.read_sra_runtable(sra_runtable_path)
match.annotate(ath_df, 'Arabidopsis thaliana', db='sra')
print("Done 😄")
