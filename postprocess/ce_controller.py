import CE

# import argparse
# import numpy as np
# import pandas as pd
# import random
# from random import shuffle
# from statsmodels.stats.multitest import multipletests
# import pdb

CE.FDR = 0.05
CE.METHOD = 2

# species_ls = ["Cpa", "Cre", "Pab", "Osa", "Sly", "Smo", "Vvi", "Zma", "Ath", "Nta", "Egu"]
# species_ls = ["Cpa", "Cre", "Pab", "Osa", "Sly", "Smo", "Vvi", "Zma"]
species_ls = ["Ath"]
for species in species_ls:
    irt = f"../plants/{species}/"
    # ort = f"CE_test/more_tests/{species}/"
    CE.mapman_path = irt + f"{species}.mercator.txt"
    CE.hcca_path = irt + f"{species}.hcca"
    CE.tpm_path = irt + f"{species}_matrix.txt"
    CE.overrep_path = irt + f"{species}_overrep_bins.txt"
    CE.underrep_path = irt + f"{species}_underrep_bins.txt"
    print(species, "Method", CE.METHOD)
    CE.main()
    print()

print("Done")
