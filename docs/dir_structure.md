# E. Directory structure

> This segment provides an overview of the structure of the directories in the main scripts `plants-pipeline` directory and the data directory `pipeline-data`.

## Main scripts directory (`plants-pipeline`)



## Data directory (`pipeline-data`)

The directory tree is as such:

pipeline-data
├── download
│   ├── bash-jobfiles
│   ├── bash-tmp
│   ├── cds
│   ├── fastq-tmp
│   ├── idx
│   ├── kallisto-out
│   ├── logs
│   │   ├── initiation
│   │   ├── progress
│   │   ├── runinfo
│   │   ├── runtime
│   │   └── status
│   └── runinfo-main
├── postprocess
│   ├── f1-stats
│   ├── gene-classifications
│   ├── pcc-matrices
│   ├── percentage-matrices
│   ├── qc-jointplots
│   ├── qc-matrices
│   ├── qc-out
│   ├── qc-summary
│   └── tpm-matrices
└── preprocess
    ├── ena-annotations
    ├── ena-runtables
    ├── job-list
    ├── species-list
    ├── sra-annotation2s
    ├── sra-annotations
    └── sra-runtables
