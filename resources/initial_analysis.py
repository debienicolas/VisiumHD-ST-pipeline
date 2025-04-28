import argparse
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import subprocess
import os
import json
import numpy as np
import sys
from pathlib import Path
from loguru import logger as logging 
import os
import anndata as ad
import seaborn as sns
from scipy.stats import median_abs_deviation
import scipy.sparse as sp
import json
from tqdm import tqdm


# Get absolute path to the root directory
root_dir = str(Path(__file__).resolve().parent.parent)
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
from resources.utils.utils import save_visium_data

""" 
See https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
"""
# scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False    
)


def create_overview(anndata:ad.AnnData, output_path:str):
    # create an overview output dir
    overview_output_dir = os.path.join(output_path, "overview")
    os.makedirs(overview_output_dir, exist_ok=True)
    
    
    ### Highly variable genes
    



def plot_qc_metrics(adata: ad.AnnData, output_path: str):
    # plot the 3 QC covariates
    # plot the total counts
    p1 = sns.histplot(adata.obs['total_counts'], bins=100, kde=False)
    p1.figure.savefig(os.path.join(output_path, 'total_counts.png'))
    p2 = sc.pl.violin(adata, 'pct_counts_mt', show=False)
    p2.figure.savefig(os.path.join(output_path, 'pct_counts_mt.png'))
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color='pct_counts_mt', show=False)
    p3.figure.savefig(os.path.join(output_path, 'total_counts_n_genes.png'))
    
    return

def is_outlier(adata: ad.AnnData, metric: str, nmads:int=5):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (M > np.median(M) + nmads * median_abs_deviation(M))
    return outlier

def save_adata(adata: ad.AnnData, path: str):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype(str)
    for var in adata.var.columns:
        if adata.var[var].dtype == 'object':
            adata.var[var] = adata.var[var].astype(str)
            
    if "spatial" in adata.obsm_keys():
        adata.obsm["spatial"] = np.array(adata.obsm["spatial"], dtype=np.float32)
        
    # check if the X is a sparse matrix if not convert it to a sparse matrix
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    
    adata.write_h5ad(path, compression="gzip")

def perform_initial_analysis(input_path: str, output_path: str, species:str, config:dict):
    
    MAD_folds = config["MAD_folds"]
    pct_mt_threshold = config["pct_mt_threshold"]
    min_cells = config["min_cells"]
    min_counts = config["min_counts"]
    
    os.makedirs(output_path, exist_ok=True)
    
    logging.info(f"Loading data from {input_path}")
    
    # check it the tissue_positions_list.csv exists
    tissue_positions_path = os.path.join(input_path, 'spatial', 'tissue_positions_list.csv')
    if not os.path.exists(tissue_positions_path):
        logging.info(f"Tissue positions file not found at {tissue_positions_path}")
        # turn a parquet into a csv
        tissue_positions_parq = os.path.join(input_path, 'spatial', 'tissue_positions.parquet')
        tissue_positions_raw = os.path.join(input_path, 'spatial', 'tissue_positions.csv.gz')
        tissue_positions_csv = os.path.join(input_path, 'spatial', 'tissue_positions_list.csv')
        # convert the parquet to a csv and gunzip it to 
        conversion_cmd = f"parquet-tools csv {tissue_positions_parq} | gzip -c > {tissue_positions_raw}"
        conversion_result = subprocess.run(conversion_cmd, shell=True, capture_output=True)
        subprocess.run(f"gunzip -c {tissue_positions_raw} > {tissue_positions_csv}", shell=True, capture_output=True)
        logging.info(f"Conversion command returned: {conversion_result.stdout.decode('utf-8')}")


    adata = sc.read_visium(input_path)
    # make the gene names unique to prevent issues with duplicate gene names
    adata.var_names_make_unique()
    logging.info(f"Data loaded: ")
    logging.info(adata)
    
    save_adata(adata, os.path.join(output_path, 'raw_adata.h5ad'))
    
    # check if the X is a sparse matrix if not convert it to a sparse matrix
    if not sp.issparse(adata.X):
        logging.info("Converting X to a sparse matrix")
        adata.X = sp.csr_matrix(adata.X)
    
    ### Filtering low quality cells ###
    # 3 QC covariates: 1. number of counts per barcode, 2. number of genes per barcode, 3. The fraction of counts from mitochondrial genes per barcode
    
    if species == "human":
        # mitochondrial genes ( in human: MT-*, in mouse: mt-*)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        # ribosomal genes
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
        # hemoglobin genes
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    elif species == "mouse":
        # mitochondrial genes ( in human: MT-*, in mouse: mt-*)
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
        # ribosomal genes
        adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))
        # hemoglobin genes
        adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")
    
    # calculate the QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True, percent_top=[20], log1p=True)
    logging.info("QC metrics calculated")
    logging.info(adata)
    
    # plot the 3 QC covariates
    plot_qc_metrics(adata, output_path)
    
    # automatic filtering based on MAD
    adata.obs["outlier"] = (
        is_outlier(adata, 'log1p_total_counts', MAD_folds) |
        is_outlier(adata, 'log1p_n_genes_by_counts', MAD_folds) |
        is_outlier(adata, 'pct_counts_in_top_20_genes', MAD_folds) |
        is_outlier(adata, 'pct_counts_mt', MAD_folds) |
        adata.obs["pct_counts_mt"] > 8 # percentage of mitochondrial counts greater than 8 
    )
    
    logging.info("outlier cells:")
    logging.info(adata.obs.outlier.value_counts())
    
    logging.info(f"Total number of cells: {adata.n_obs}")
    adata = adata[~adata.obs["outlier"]]
    logging.info(f"Number of cells after outlier removal: {adata.n_obs}")
    
    # plot the scatter plot of the filtered data
    p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color='pct_counts_mt', show=False)
    p1.figure.savefig(os.path.join(output_path, 'total_counts_n_genes_filtered.png'))
    
    
    ### Filtering low quality genes###
    logging.info(f"Total number of genes: {adata.n_vars}")
    # min 20 cells 
    sc.pp.filter_genes(adata, min_cells=args.config["min_cells"])
    logging.info(f"Number of genes after filtering: {adata.n_vars}, min cells: {args.config['min_cells']}")
    
    ### Filtering low quality cells ###
    logging.info(f"Total number of cells: {adata.n_obs}")
    # min 5 counts per cells
    sc.pp.filter_cells(adata, min_counts=args.config["min_counts"]) 
    logging.info(f"Number of cells after filtering: {adata.n_obs}, min counts: {args.config['min_counts']}")
    
    ### normalization ### 
    logging.info("Normalizing data")
    adata.raw = adata.copy()
    # check for nan values in the counts matrix
    # if np.isnan(adata.X.toarray()).any():
    #     logging.info("Nan values found in the counts matrix")
    # else:
    #     logging.info("No nan values found in the counts matrix")
   
    # normalize the data in batches
    batch_size = 10_000
    for i in tqdm(range(0, adata.n_obs, batch_size)):
        sc.pp.normalize_total(adata[i:i+batch_size], inplace=True, exclude_highly_expressed=True, target_sum=10_000)
        sc.pp.log1p(adata[i:i+batch_size]) # why log1p? -> handles zero counts, variance stabilizing, scale compression
    
    logging.info("Final adata:")
    logging.info(adata)
        
    # save the filtered data
    logging.info(f"Saving data to {os.path.join(output_path, 'filtered_adata.h5ad')}")
    save_adata(adata, os.path.join(output_path, 'filtered_adata.h5ad'))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform initial analysis on binned Visium data")
    parser.add_argument("--input_path", type=str, required=False, help="Path to the input binned Visium data")
    parser.add_argument("--output_path", type=str, required=False, help="Path to the output directory")
    parser.add_argument("--species", type=str, required=False, help="Species of the data")
    parser.add_argument("--config", type=str, required=False, help="Path to the config file")
    args = parser.parse_args()
    
    args.config = json.loads(args.config)
    
    logging.info(f"Config:")
    for key, value in args.config.items():
        logging.info(f"{key}: {value}")
    
    if args.config is None:
        # default config
        args.config = {
            "MAD_folds": 5,
            "pct_mt_threshold": 8,
            "min_cells": 20,
            "min_counts": 5
        }

    if args.input_path is None:
        args.input_path = "input/TG23-0227_VGEX_results_VISIUM_HD_reanalyzed_with_highquality_HE/binned_outputs/square_002um"
        args.input_path = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/binned_outputs/square_016um"
        args.output_path = "results/test_initial_analysis"
        args.species = "mouse"
        args.species = "human"
    
    perform_initial_analysis(args.input_path, args.output_path, args.species, args.config)


