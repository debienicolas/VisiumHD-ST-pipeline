"""
This script is used to run infercnvpy using the ficture output for cell types
"""
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import argparse
import pyensembl
import numpy as np
import sys
from pathlib import Path
from loguru import logger as logging
import os
import squidpy as sq
import scanpy as sc
import squidpy as sq
import scanpy as sc
import subprocess
import json


PYENSEMBL_CACHE_DIR = 'resources_data/'

# set the cache dir for pyensemble to the resources_data dir
os.environ["PYENSEMBL_CACHE_DIR"] = PYENSEMBL_CACHE_DIR

PLOT_POINT_SIZE = 10



def add_chromosome_info(visium_data:ad.AnnData, species:str) -> ad.AnnData:
    """
    Add chromosome information to the data.
    
    Args:
        visium_data: the anndata object to add chromosome information to
    
    Returns:
        visium_data: the anndata object with chromosome information added: specifically extra columns in the var dataframe: chromosome, start, end
    """
    # Create a copy to avoid view issues
    #visium_data = visium_data.copy()
    
    # if gene_names column is present, don't create it
    if "gene_names" not in visium_data.var.columns:
        logging.info("Adding gene_names column")
        visium_data.var["gene_names"] = visium_data.var.index.tolist()
        logging.info("gene_names column added from index")
        logging.info("visium_data.var.columns: ")
        logging.info(visium_data.var.columns)
    
    
    ### Add chromosome info
    
    # check if the pyensembl database is installed
    if not os.path.exists(PYENSEMBL_CACHE_DIR + "/pyensembl/GRCm38/ensembl102") and species == "mouse":
        logging.info("Pyensembl database not found, installing it")
        # cache dir is already set at the top of the script
        subprocess.run(["pyensembl", "install", "--release", "102", "--species", "mouse"])
    
    if not os.path.exists(PYENSEMBL_CACHE_DIR + "/pyensembl/GRCh38/ensembl102") and species == "human":
        logging.info("Pyensembl database not found, installing it")
        # cache dir is already set at the top of the script
        subprocess.run(["pyensembl", "install", "--release", "76", "--species", "human"])
    
    if species == "mouse":
        ensembl = pyensembl.EnsemblRelease(102, species="mus_musculus")
    else:
        ensembl = pyensembl.EnsemblRelease(76, species="homo_sapiens")
    
    gene_ids = visium_data.var.gene_ids.tolist()
    logging.info("gene_ids: ")
    logging.info(gene_ids[:5])
    logging.info("Loading the ensemble database")

    # Create the positions dataframe
    gene_positions = {}
    for gene_id in gene_ids:
        try:
            gene = ensembl.gene_by_id(gene_id)
            gene_positions[gene_id] = [gene.contig, gene.start, gene.end]
        except Exception as e:
            logging.info(f"Error adding gene {gene_id}: {e}")
            gene_positions[gene_id] = [None, None, None]
    
    # Map the positions using the gene_ids as the key
    visium_data.var["chromosome"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][0])
    visium_data.var["start"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][1])
    visium_data.var["end"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][2])
    
    logging.info("Added chromosome, start, end columns")
    logging.info("visium_data.var.columns: ")
    logging.info(visium_data.var.columns)
    logging.info("visium_data.var.head(): ")
    logging.info(visium_data.var.head())
    
    logging.info("Removing entries with nan in chromosome column")
    visium_data = visium_data[:, visium_data.var["chromosome"].notna()]
       
    # Add chr prefix required for infercnvpy
    visium_data.var["chromosome"] = visium_data.var["chromosome"].astype(str)
    visium_data.var["chromosome"] = visium_data.var["chromosome"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    logging.info(f"Chromosome values: ")
    logging.info(visium_data.var['chromosome'].value_counts())

    logging.info("Adding chr prefix (required for infercnvpy)")
    visium_data.var["chromosome"] = visium_data.var["chromosome"].astype(str)
    visium_data.var["chromosome"] = visium_data.var["chromosome"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
    
    logging.info("Chromosome information added")
    logging.info("visium_data.var.columns: ")
    logging.info(visium_data.var.columns)
    
    return visium_data

def annotate_normal(visium_data:ad.AnnData) -> ad.AnnData:
    """
    Annotate the normal cells and tumor cells in the data.

    Biological assumption: cells with Foxa1 are normal cells
    Look at the cooccurrence of other genes in cells that have Foxa1. 
    If the cooccurrence is found in other cells without Foxa1, then the cell is also considered normal.
    
    Args:
        visium_data: the anndata object to annotate
    
    Returns:
        visium_data: the anndata object with the normal/tumor annotation added to the obs dataframe. The obs column "cell_label" is added which has the value "normal" or "tumor".
    
    Exceptions:
        ValueError: if Foxa1 is not in the data
    
    """
    logging.info("Annotating normal cells")
    logging.info("visium_data.var_names: ")
    logging.info(visium_data.var_names)
    # print the var names with fox in them or Fox 
    
    # mark all the factors that have immune in them as normal
    mask = visium_data.obs["factor"].str.contains("immune", case=False)
    visium_data.obs.loc[mask, 'cell_label'] = 'normal'
    logging.info(f"Tumor/Normal ratio: ")
    logging.info(visium_data.obs['cell_label'].value_counts())
    
    # mark a cell as normal if the foxa1 is expressed
    #logging.info("Foxa1 in var_names: ")
    #logging.info([name for name in visium_data.var_names if "Foxa1" in name or "foxa1" in name or "fox" in name])
    #if "Foxa1" not in visium_data.var_names:
    #    raise ValueError("Foxa1 is not in the data")
    # mask = visium_data[:,'Foxa1'].X.toarray().flatten() > 0
    # visium_data.obs['cell_label'] = 'tumor'  # First set all to tumor
    # visium_data.obs.loc[mask, 'cell_label'] = 'normal'  # Then set normal where Foxa1 > 0
    # logging.info(f"Tumor/Normal ratio: ")
    # logging.info(visium_data.obs['cell_label'].value_counts())
    
    visium_data.obs["cell_label"] = visium_data.obs["cell_label"].astype("category")
    logging.info("Cell annotation done")
    return visium_data


def run_infercnvpy_on_individual_factors(visium_data:ad.AnnData, output_path:str, window_size:int, step:int, leiden_res:float, infercnv_clip:float):
    """
    Run infercnv on the entire anndata object, run it on each factor separately to find clones within each factor
    """


    for factor in visium_data.obs["factor"].unique():
        
        print("factor: ", factor)
        print(" visium data head: ", visium_data.obs.head())
        print(" visium data counts: ", visium_data.X.toarray().sum())
        # create the output path for the current factor
        os.makedirs(output_path + f"/{factor}", exist_ok=True)
        # select the anndata object for the current factor
        print("factor value counts: ", visium_data.obs["factor"].value_counts())
        print("Factor mask shape: ", visium_data[visium_data.obs["factor"] == factor].shape)
        # Create mask and check its properties
        print("Original visium data shape:", visium_data.shape)
        print("Factor being processed:", factor)
        
        # Create boolean mask
        mask = visium_data.obs["factor"] == factor
        print(f"Found {mask.sum()} cells for factor {factor}")
        
        factor_data = visium_data[mask]
        
        print("New factor data shape:", factor_data.shape)
        print("Var columns:", factor_data.var.columns.tolist())
        # check that the following columns are present and not empty
        assert "chromosome" in factor_data.var.columns, "chromosome column is missing"
        assert "start" in factor_data.var.columns, "start column is missing"
        assert "end" in factor_data.var.columns, "end column is missing"
        assert "cell_label" in factor_data.obs.columns, "cell_label column is missing"
        
        print("factor_data: ", factor_data)
        logging.info("Normal cells in factor_data: ")
        logging.info(factor_data.obs["cell_label"].value_counts())
        print("Normal cells in factor_data: ")
        print(factor_data.obs["cell_label"].value_counts())

        # perform infercnv
        if "normal" in factor_data.obs["cell_label"].unique():
            cnv.tl.infercnv(
                factor_data,
                reference_key="cell_label",
                window_size=window_size,
                reference_cat="normal",
                step=step
            )
        else:
            cnv.tl.infercnv(
                factor_data,
                window_size=window_size,
                step=step
            )
        print("post infercnvpy data: ", factor_data)
        # clustering by CNV profiles
        cnv.tl.pca(factor_data)
        cnv.pp.neighbors(factor_data)
        cnv.tl.leiden(factor_data, resolution=leiden_res)
    
        cnv.tl.umap(factor_data)
        cnv.tl.cnv_score(factor_data)
        logging.info(f"post infercnvpy + umap + cnv_score data: \n {factor_data}")

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
        ax4.axis("off")
        cnv.pl.umap(
            factor_data,
            color="cnv_leiden",
            legend_loc="on data",
            legend_fontoutline=2,
            ax=ax1,
            show=False,
        )
        cnv.pl.umap(factor_data, color="cnv_score", ax=ax2, show=False)
        cnv.pl.umap(factor_data, color="factor", ax=ax3, show=False)
        fig.savefig(output_path + f"/{factor}/umap.png")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 11))
        sc.pl.spatial(factor_data, color="cnv_leiden", ax=ax1, show=False)
        sc.pl.spatial(factor_data, color="factor", ax=ax2, show=False)
        # add some space between the plots
        fig.savefig(output_path + f"/{factor}/spatial.png")
        print("cnv_leiden: ", factor_data.obs["cnv_leiden"].value_counts())
        print("cnv_score: ", factor_data.obs["cnv_score"].value_counts())

        print("post infercnvpy data: ", factor_data)
        #cnv.pl.chromosome_heatmap(factor_data, groupby="factor", show=False)
        #plt.savefig(output_path + f"/{factor}/chromosome_heatmap.png")
        #plt.close()

        # convert the cell_type column to a categorical variable
        factor_data.obs['factor'] = factor_data.obs['factor'].astype(str).astype("category")

        # clip the cnv values to be between -0.25 and 0.25
        factor_data.obsm["X_cnv"].data = np.clip(factor_data.obsm["X_cnv"].data, -infercnv_clip, infercnv_clip)

        #sc.tl.dendrogram(factor_data, groupby="factor")
        fig = cnv.pl.chromosome_heatmap(factor_data, groupby="factor", dendrogram=False, show=False)
        # save the heatmap
        for key, ax in fig.items():
            ax.figure.savefig(output_path + f"/{factor}/{key}.png")

        print(f"Factor {factor} done")

        # add the clone information to the visium anndata object, this is going from a smaller anndata object to the full visium anndata object
        # set the factor subset of 
        visium_data[visium_data.obs["factor"] == factor].obs["clone_individual"] = factor_data.obs["cnv_leiden"]

    return visium_data

def main_cnv_logic(visium_data:ad.AnnData, output_path:str, window_size:int, step:int, leiden_res:float, infercnv_clip:float):
    """
    Main logic for infercnvpy. 
    
    Args:
        visium_data: the anndata object to run infercnvpy on
        output_path: the path to the output directory
    
    Returns:
        None
    """
    # check if there is a normal cell in the data
    if "normal" in visium_data.obs["cell_label"].unique():
        # perform infercnv
        cnv.tl.infercnv(
            visium_data,
            reference_key="cell_label",
            window_size=window_size,
            reference_cat="normal",
            step=step
        )
    else:
        # perform infercnv without reference
        cnv.tl.infercnv(
            visium_data,
            window_size=window_size,
            step=step
        )
        
    # clustering by CNV profiles
    cnv.tl.pca(visium_data)
    cnv.pp.neighbors(visium_data)
    cnv.tl.leiden(visium_data, resolution=leiden_res)

    cnv.tl.umap(visium_data)
    cnv.tl.cnv_score(visium_data)
    
    logging.info("Plotting the cnv umap plot")
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        visium_data,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(visium_data, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(visium_data, color="factor", ax=ax3, show=False)
    fig.savefig(output_path + "/umap.png")

    logging.info("Plotting the cnv spatial plot")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 11))
    sc.pl.spatial(visium_data, color="cnv_leiden", ax=ax1, show=False)
    sc.pl.spatial(visium_data, color="factor", ax=ax2, show=False)
    fig.savefig(output_path + "/spatial.png")
   
    logging.info("Clipping the cnv values to be between -0.25 and 0.25")
    visium_data.obsm["X_cnv"].data = np.clip(visium_data.obsm["X_cnv"].data, -infercnv_clip, infercnv_clip)

    logging.info("Plotting the cnv chromosome heatmap")
    #sc.tl.dendrogram(visium_data, groupby="factor")
    cnv.pl.chromosome_heatmap(visium_data, groupby="factor", dendrogram=False, show=False)
    plt.savefig(output_path + "/chromosome_heatmap.png")
    
    return visium_data


def run_infercnvpy(input_path, output_path, species, config):
    """
    Run infercnvpy on the ficture annotated output.
    Does both global and local(individual factors) infercnv. 

    Args:
        input_path: path to the ficture annotated anndata object
        output_path: path to the output directory
    
    Returns:
        None
        
    """
    if species not in ["mouse", "human"]:
        raise ValueError("Species must be either mouse or human")
    
    ## parse the config
    window_size = config.get("window_size", 100)
    step = config.get("step", 3)
    leiden_res = config.get("leiden_res", 0.5)
    infercnv_clip = config.get("infercnv_clip", 0.25)

    os.makedirs(output_path, exist_ok=True)

    ### Visium HD data ###

    # read the anndata object with annotated ficture clusters + set the factor from float to int for plotting
    visium_data = ad.read_h5ad(input_path)
    visium_data.obs["factor"] = visium_data.obs["factor"].astype(str).astype("category")
    logging.info("Visium data: ")
    logging.info(visium_data)
    
    # normalize the counts and log transform
    sc.pp.normalize_total(visium_data, target_sum=1e4)
    sc.pp.log1p(visium_data)
    
    ### ADD CHROMOSOME INFORMATION ###
    logging.info("Adding chromosome information")
    visium_data = add_chromosome_info(visium_data, species)
    
    ### CELL ANNOTATION ###
    logging.info("Annotating normal cells")
    visium_data = annotate_normal(visium_data)
    # if "Foxa1" in visium_data.var_names:
    #     logging.info("Foxa1 is in the data")
    #     visium_data = annotate_normal(visium_data)
    # else:
    #     visium_data.obs["cell_label"] = "tumor"
    #     logging.info("Foxa1 is not in the data, skipping cell annotation and using average")

    ### LOCAL INFERCNV ###
    run_infercnvpy_on_individual_factors(visium_data, output_path, window_size, step, leiden_res, infercnv_clip)
    
    ### GLOBAL INFERCNV ###
    visium_data = main_cnv_logic(visium_data, output_path, window_size, step, leiden_res, infercnv_clip)

    # save the anndata object
    visium_data.write_h5ad(output_path + "/infercnvpy_output.h5ad")

   
        
    
    return None

def run_infercnvpy_benchmarking(visium_data:ad.AnnData, output_path:str):
    """
    Run infercnvpy on the anndata object and save the output.
    """
    os.makedirs(output_path, exist_ok=True)

    ### Visium HD data ###

    # read the anndata object with annotated ficture clusters + set the factor from float to int for plotting
    visium_data.obs["factor"] = visium_data.obs["factor"].astype(str).astype("category")

    # normalize the counts and log transform
    sc.pp.normalize_total(visium_data, target_sum=1e4)
    sc.pp.log1p(visium_data)
    
    ### ADD CHROMOSOME INFORMATION ###
    logging.info("Adding chromosome information")
    visium_data = add_chromosome_info(visium_data)
    
    ### CELL ANNOTATION ###
    logging.info("Annotating normal cells")
    visium_data = annotate_normal(visium_data)
        

    run_infercnvpy_on_individual_factors(visium_data, output_path)
    
    ### GLOBAL INFERCNV ###
    visium_data = main_cnv_logic(visium_data, output_path)

    # save the anndata object
    visium_data.write_h5ad(output_path + "/infercnvpy_output.h5ad")

    ### LOCAL INFERCNV ### 
    for factor in visium_data.obs["factor"].unique():
        os.makedirs(output_path + f"/{factor}", exist_ok=True)
        factor_data = visium_data[visium_data.obs["factor"] == factor]
        main_cnv_logic(factor_data, output_path + f"/{factor}")
        


if __name__ == "__main__":

    # use argparse to get the input path
    parser = argparse.ArgumentParser(description="Run infercnvpy on the ficture output")
    parser.add_argument("--input", type=str, help="Path to the anndata visium object")
    parser.add_argument("--output", type=str, help="Path to the infercnvpy output directory")
    parser.add_argument("--species", type=str, help="Sample species")
    parser.add_argument("--config", type=str, help="Path to the infercnvpy config file")
    args = parser.parse_args()

    # if args.input is None or args.output is None:
    #     args.input = "results/run_4/ficture/postprocessed_ficture/ficture_anndata.h5ad"
    #     args.output = "results/run_4/infercnvpy"
    args.config = json.loads(args.config)

    if args.input is None or args.output is None:
        print("Using manual input and output")
        results_dir = "human_breast_cancer_final"
        args.input = f"results/{results_dir}/ficture/annotated_anndata.h5ad"
        args.output = f"results/{results_dir}/infercnvpy"
        args.species = "human"
    print("input: ", args.input)
    print("output: ", args.output)

    run_infercnvpy(args.input, args.output, args.species, args.config)


