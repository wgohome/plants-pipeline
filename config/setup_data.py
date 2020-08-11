import os
import argparse

parser = argparse.ArgumentParser(description = 'This script sets up the directory structure of the data repository for the plants-pipeline scripts.', epilog = 'By Mutwil Lab')
parser.add_argument('-p', '--path', nargs=1, metavar='data_path',
                    help='Enter the path of the directory where the data directories for this project is to be created',
                    dest='data_path', type=str, required=True)
args = parser.parse_args()
DATA_PATH = args.data_path[0] + '/pipeline-data/'

# @HERE: TO CONFIGURE
dests = {
    'download': [
        'cds',
        'idx',
        'bash-jobfiles',
        'bash-tmp',
        'fastq-tmp',
        'kallisto-out',
        'logs/initiation',
        'logs/runtime',
        'logs/progress',
        'logs/runinfo'
    ],
    'preprocess': [
        'ena-runtables',
        'sra-runtables',
        'ena-annotations',
        'sra-annotations',
        'species-list'
    ],
    'postprocess': [
        ''
    ]
}

def makedirs(dest):
    os.makedirs(f"{DATA_PATH}{dest}", exist_ok=True)

for p, c in dests.items():
    for dest in c:
        makedirs(f"{p}/{dest}")
        print(f"Created: {p}/{dest}")

print("🚀 Created directories for pipeline data repository")
