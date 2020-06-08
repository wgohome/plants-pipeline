from os.path import realpath, dirname
import sys

if sys.platform == 'linux':
    DATA_PATH = "/home/data1/will_data/pipeline-data"
    ASPERA_SSH_KEY = "/home/workstation/.aspera/cli/etc/asperaweb_id_dsa.openssh" # For Linux
elif sys.platform == 'darwin':
    DATA_PATH = "/Users/wirriamm/marek_lab/plants-pipeline/data-local"
    ASPERA_SSH_KEY = "/Users/wirriamm/Applications/Aspera CLI/etc/asperaweb_id_dsa.openssh" # For MacOS
else:
    print("Please manually set your constants")
    DATA_PATH = input("DATA_PATH: ")
    ASPERA_SSH_KEY = input("ASPERA_SSH_KEY: ")
