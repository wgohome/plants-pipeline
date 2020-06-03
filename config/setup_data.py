import os
import argparse

parser = argparse.ArgumentParser(description = 'This script sets up the directory structure of the data repository for the plants-pipeline scripts.', epilog = 'By Mutwil Lab')
parser.add_argument('-p', '--path', nargs=1, metavar='data_path',
                    help='Enter the path of the directory where the data directories for this project is to be created',
                    dest='data_path', type=str, required=True)
args = parser.parse_args()
DATA_PATH = args.data_path[0]

os.makedirs(f"{DATA_PATH}/pipeline-data/download/cds", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/download/idx", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/download/fastq", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/download/kallsito_out", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/download/logs/status", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/download/logs/time", exist_ok=True)

os.makedirs(f"{DATA_PATH}/pipeline-data/preprocess/out/ena_runtables", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/preprocess/out/sra_runtables", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/preprocess/out/ena_annotations", exist_ok=True)
os.makedirs(f"{DATA_PATH}/pipeline-data/preprocess/out/sra_annotations", exist_ok=True)

os.makedirs(f"{DATA_PATH}/pipeline-data/postprocess", exist_ok=True)

print("Created directories for pipeline data repository")
