# VisiumHD-ST-pipeline

[![Python](https://img.shields.io/badge/Python-3.10.13-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.readthedocs.io/)

This spatial transcriptomics (ST) pipeline is an end-to-end Visium HD analysis pipeline implemented using Snakemake. 
It goes from raw data (fresh from spaceranger) to running some of the most popular downstream analysis packages. 

The ST data is stored in an [annotated data object](https://anndata.readthedocs.io/en/stable/) to ensure compatibility with popular single cell/ST analysis packages in python.

<div align="center">
    <img src="meta/overview_figure.png" alt="st-pipeline-overview" width="700"/>
    <p><em>Figure 1: Overview of the Visium HD spatial transcriptomics pipeline workflow</em></p>
</div>

The pipeline has four main sections:  
1. Pre-processing of *10X Genomics* Visium HD data
    - Filtering of transcripts and bins
    - Data normalisation
2. Cell segmentation and bin aggregation (optional for non $2\mu m$ res.)
    - H&E based cell segmentation using stardist model
    - Gene expression based segmentation
    - Bin aggregation: aggregating bins within same cell to single bin representing a cell.
3. Spatially informed cell type clustering & Annotation
    - unsupervised clustering (Ficture)
    - Manual annotation or annotation using reference atlas
4. Downstream analysis
    - Infercnv: Clone discovery using copy number variation
    - Cellchat: Ligand-receptor communication network analysis
    - Monocle: Cell trajectory analysis
    - Pyscenic: Gene regulatory network analysis
    - Squidpy: Spatial statistic analysis




## Quickstart

It is highly recommended to create a new conda environment:
```bash
conda create --name st python=3.10.13
conda activate st
pip install -r requirements.txt
```

1. Configure parameters in `config.yaml` file.
2. Run the pipeline:
```bash
snakemake --cores 10
```
It is highly recommended to run this pipeline on a HPC (seee the hpc_scripts).  

## Usage

[TO DO: Detailed usage instructions]  





## Repository structure

```
.
├── input                   # Directory containing Visium HD samples
├── resources               # Main code files 
├── workflow                # Snakemake rules 
├── config.yaml             # Pipeline configuration file
└── Snakefile               # Main Snakemake file

```

