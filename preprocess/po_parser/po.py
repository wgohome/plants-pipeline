# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

abspath = realpath(dirname(__file__))
parent_module = re.search('^(.*plants-pipeline)', abspath).group()
if __name__ == '__main__':
    sys.path.insert(0, parent_module)

################################################################################

import os
import re
import requests
# local imports
from .ontology import Ontology
from .term import Term

url = "https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo"
po_obo_path = f"{parent_module}/dependencies/po.obo"
if not os.path.exists(po_obo_path):
    response = requests.get(url)
    with open(po_obo_path, 'w') as f:
        f.write(response.text)

obo = open(po_obo_path).read()
po = Ontology(obo)

__all__ = ['po']
