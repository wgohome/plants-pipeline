# LSTrAP-Kingdom: an automated pipeline to generate annotated gene expression atlases for kingdoms of life

By [**NTU Plants Systems Biology and Evolution Laboratory**](https://www.plant.tools)

This repository is found in [this Github Repository](https://github.com/wirriamm/plants-pipeline), with an accompanying paper found [here (preprint version)](https://www.biorxiv.org/content/10.1101/2021.01.23.427930v1.full). Do create pull requests for issues/bugs and feature requests. [Contact me](mailto:will0046@e.ntu.edu.sg) for feedback or reporting bugs.

# Guides

[**A. First local setup of the pipeline**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/setup.md)

- This segment only needs to be implemented at the first setup of this repository on your local machine/server.

[**B. Initialization for each session**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/initialize.md)

- These commands need to be run everytime the pipeline is accessed from a new terminal session. They will load the python environment with the installed packages, and add ascp and kallisto commands to the global environment $PATH. If kallisto or ascp(Aspera CLI) is not downloaded, they will also be downloaded.

[**C & D. Download guide**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/download.md)

- **C. Bulk Download**

  - The steps to run a download job for multiple species are outlined here.

- **D. Small download job**

  - The steps to run a download job for a single species are outlined here.

[**E. Directory structure**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/dir_structure.md)

- This segment provides an overview of the structure of the directories in the main scripts `plants-pipeline` directory and the data directory `pipeline-data`.

[**F. Postprocessing**](https://github.com/wirriamm/plants-pipeline/blob/master/docs/postprocess.md)

- After the download job is completed, these are the steps needed to generate the TPM matrices and perform quality control, which includes:
  - Generating TPM matrices
  - Quality control
  - Performing coexpression to count number of ribosomal gene neighbours for every gene
- The F1 scores for the benchmark in the paper are generated using the scripts here.

[**G. Annotation Benchmark**]()
