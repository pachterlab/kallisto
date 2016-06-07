---
layout: page
title: "Getting started with 10X Chromium 3' digital expression analysis"
description: ""
---

The short tutorial below explains how to process a small (example) 10X Chromium 3' digital expression data. 

- Begin by making sure you have the (correct versions) of the necessary programs installed:
  - [Install __kallisto__](download.html). 
  - Install [python](https://www.python.org). It is convenient to work on a linux server for the processing of single-cell data (due to large files and the utility of multithreading) and an easy way to install python is via [Anaconda](https://www.continuum.io/downloads). Python 2.7 suffices and can be installed with `bash Anaconda2-4.0.0-Linux-x86_64.sh`. 
  - Install [Jupyter Notebook](http://jupyter.readthedocs.io/en/latest/install.html) (this will not be necessary if you've installed python with Anaconda as Jupyter Notebook is bundled along with other packages and libraries you will need). To run Jupyter Notebook remotely on a server type `ipython notebook --no-browser --port=8889`. Start an SSH tunnel on the local machine with `ssh -N -f -L localhost:8888:localhost:8889 remote_user@remote_host`. You can then run Juypter Notebook remotely via a local web-browser at `http://localhost:8888`.
  - Make sure that the installed version of `scikit-learn` is __at most__ 0.16.1. The latest version 0.17.1 has problems with the tSNE algorithm that is called in some of the Jupyter Notebooks you will use.  
<br />
- Download the "human-mouse" transcriptome from the [kallisto transcriptome website](http://bio.math.berkeley.edu/kallisto/transcriptomes/).

- Build the __kallisto__ index for the transcriptome with the command `kallisto index -i transcripts.idx Mus_musculus.GRCm38.rel79_Homo_sapiens.GRCh38.rel79.cdna.mix.fa.gz`.

- Copy the `config.json` file from [here](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/example_dataset) to a working directory of choice and edit the paths in the file as needed. If necessary also change the NUM_THREADS parameter to a value appropriate for the machine on which the example is to be processed.

- Copy the `10xDetect_and_Prep.py` script to your working directory. 

- Type `python 10xDetect_and_Prep.py config.json`.

- Copy the [`10x_results.ipynb`](https://github.com/lakigigar/scRNA-Seq-TCC-prep/blob/master/notebooks/10xResults.ipynb) from the [notebook folder](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/notebooks) to your working directory where the `config.json` file is located.

- Run the notebook locally by typing `jupyter notebook 10x_results.ipynb` or see instructions above for working with a notebook remotely.

- The result obtained should be identical to the one in [`10xResults_hgmm_small.ipynb`](https://github.com/lakigigar/scRNA-Seq-TCC-prep/blob/master/example_dataset/10xResults_hgmm_small.ipynb).

- We strongly recommend performing the first step of the workflow interactively using the [`10xGet_cell_barcodes.ipynb`](https://github.com/lakigigar/scRNA-Seq-TCC-prep/blob/master/example_dataset/10xGet_cell_barcodes_hgmm_small.ipynb). See the [README](https://github.com/lakigigar/scRNA-Seq-TCC-prep/blob/master/README.md) for how to do this with your own data.


{% include JB/setup %}


