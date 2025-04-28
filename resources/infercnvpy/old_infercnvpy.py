import infercnvpy as cnv
import matplotlib.pyplot as plt
import matplotlib.style as style
import anndata as ad
import pandas as pd
from infercnvpy.io import genomic_position_from_gtf
import argparse
import pyensembl
import numpy as np
import sys
from pathlib import Path
import os
import squidpy as sq
import scanpy as sc
import sys
from pathlib import Path
root_dir = str(Path(__file__).resolve().parent.parent.parent)
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
from resources.utils.profiler import comprehensive_profile
import multiprocessing as mp
# making matplotlib faster
plt.ioff()
style.use("fast")

from loguru import logger as logging
import bbknn

""" 
Resource allocation:
- extra cpus are useful as the leiden clustering can be parallelized

"""


def add_chromosome_info(visium_data:ad.AnnData):
    """
    Add chromosome information to the data.

    Args:
        data (ad.AnnData): The anndata object to add chromosome information to.

    Returns:
        ad.AnnData: The anndata object with chromosome information added.
    """
    
    logging.info("Starting chromosome info addition")
    visium_data.var["gene_names"] = visium_data.var.index.tolist()
    

    # CHROMOSOME, START and END  
    logging.info("Loading the ensembl database")
    ensembl = pyensembl.EnsemblRelease(102, species="mus_musculus")
    gene_ids = visium_data.var.gene_ids.tolist()

    # loop over the gene ids and get the chromosome, start and end position
    gene_positions = {}
    for gene_id in gene_ids:
        try:
            gene = ensembl.gene_by_id(gene_id)
            gene_positions[gene_id] = [gene.contig, gene.start, gene.end]
        except:
            gene_positions[gene_id] = [None, None, None]
    logging.info("Finished calling the ensembl database")

    # Map the positions using the gene_ids as the key
    visium_data.var["chromosome"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][0])
    visium_data.var["start"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][1])
    visium_data.var["end"] = visium_data.var["gene_ids"].map(lambda x: gene_positions[x][2])


    # remove all the entries with nan in the chromosome column
    mask = visium_data.var["chromosome"].notna()
    visium_data = visium_data[:, mask] 

    # Add chr prefix required for infercnvpy
    visium_data.var["chromosome"] = visium_data.var["chromosome"].astype(str)
    visium_data.var["chromosome"] = visium_data.var["chromosome"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    logging.info(f"Chromosome values: \n {visium_data.var['chromosome'].value_counts()}")

    return visium_data

def annotate_normal(visium_data:ad.AnnData) -> ad.AnnData:
    """
    Annotate the normal cells and tumor cells in the data.

    Biological assumption: cells with Foxa1 are normal cells
    Look at the cooccurrence of other genes in cells that have Foxa1. 
    If the cooccurrence is found in other cells without Foxa1, then the cell is also considered normal.
    """


    logging.info("Annotating normal cells")
    logging.info("visium_data.var_names: ", visium_data.var_names)
    # print the var names with fox in them or Fox 
    logging.info("Foxa1 in var_names: ", [name for name in visium_data.var_names if "Foxa1" in name or "foxa1" in name or "fox" in name])
    if "Foxa1" not in visium_data.var_names:
        raise ValueError("Foxa1 is not in the data")
    
    # mark a cell as normal if the foxa1 is expressed
    mask = visium_data[:,'Foxa1'].X.toarray().flatten() > 0
    visium_data.obs['cell_label'] = 'tumor'  # First set all to tumor
    visium_data.obs.loc[mask, 'cell_label'] = 'normal'  # Then set normal where Foxa1 > 0
    logging.info(f"Tumor/Normal ratio: {visium_data.obs['cell_label'].value_counts()}")

    logging.info("Cell annotation done")
    return visium_data


def process_factor(factor, visium_data:ad.AnnData, output_path:str):
    """
    Process a single factor
    """
    
    # create the output path for the current factor
    os.makedirs(output_path + f"/{factor}", exist_ok=True)
    
    # select the anndata object for the current factor
    factor_data = visium_data[visium_data.obs["factor"] == factor]

    # perform infercnv
    logging.info("Performing infercnv for factor: ", factor)
    cnv.tl.infercnv(factor_data, reference_key="cell_label", window_size=100,reference_cat="normal",step=3)
    logging.info("Finished infercnv")
    
    # clustering by CNV profiles
    logging.info("Performing PCA, neighbors, leiden, umap and cnv score")
    cnv.tl.pca(factor_data)
    logging.info("Finished PCA")
    #cnv.pp.neighbors(factor_data, transformer="pynndescent", n_neighbors=15)
    # create a batch column in the obs dataframe, they should be random batchesn of 10000 cells
    #bbknn.bbknn(factor_data, use_rep="X_cnv_pca")
    cnv.pp.neighbors(factor_data)
    logging.info("Finished neighbors")
    cnv.tl.leiden(factor_data, resolution=0.3, flavor="igraph", n_iterations=2,directed=False)
    logging.info("Finished leiden")
    cnv.tl.umap(factor_data)
    logging.info("Finished UMAP")
    cnv.tl.cnv_score(factor_data)
    logging.info("Finished cnv score")

    # Create and save UMAP plots for factor data
    fig = plt.figure(figsize=(11, 11))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)  
    ax3 = plt.subplot(223)
    cnv.pl.umap(
        factor_data, color="cnv_leiden", legend_loc="on data",
        legend_fontoutline=2, ax=ax1, show=False,
    )
    cnv.pl.umap(factor_data, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(factor_data, color="factor", ax=ax3, show=False)
    fig.savefig(output_path + f"/{factor}/umap.png")
    plt.close(fig)
    
    # create spatial plots
    fig = plt.figure(figsize=(11, 11))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    sc.pl.spatial(factor_data, color="cnv_leiden", size=10, ax=ax1, show=False, spot_size=1)
    sc.pl.spatial(factor_data, color="factor", size=10, ax=ax2, show=False, spot_size=1)
    fig.savefig(output_path + f"/{factor}/spatial.png")
    plt.close(fig)
    
    

    # clip the cnv values to be between -0.25 and 0.25
    factor_data.obsm["X_cnv"].data = np.clip(factor_data.obsm["X_cnv"].data, -0.25, 0.25)

    #sc.tl.dendrogram(factor_data, groupby="factor")
    fig = cnv.pl.chromosome_heatmap(factor_data, groupby="factor", dendrogram=False, show=False)
    
    # save the heatmap
    for key, ax in fig.items():
        ax.figure.savefig(output_path + f"/{factor}/{key}.png")

    logging.info(f"Factor {factor} done")
    
    return 


def run_infercnvpy_on_individual_factors(visium_data:ad.AnnData, output_path:str):
    """
    Run infercnv on the entire anndata object, run it on each factor separately to find clones within each factor
    """
    factors = visium_data.obs["factor"].unique().tolist()
    for i, factor in enumerate(factors):
        logging.info(f"Processing factor {i+1} of {len(factors)}: {factor}")
        process_factor(factor, visium_data, output_path)
    return visium_data


def sample_data(visium_data:ad.AnnData):
    """
    Take a representative subset of the data by creating a boolean mask.
    Samples cells_per_factor from each factor group.
    
    Returns a subset of the original AnnData object.
    """
    cells_per_factor = 10_000
    
    # Create an empty mask
    final_mask = np.zeros(visium_data.n_obs, dtype=bool)
    
    # For each factor, randomly select cells and add to mask
    for factor in visium_data.obs["factor"].unique():
        factor_mask = visium_data.obs["factor"] == factor
        factor_indices = np.where(factor_mask)[0]
        
        # If we have more cells than needed, randomly sample them
        if len(factor_indices) > cells_per_factor:
            selected_indices = np.random.choice(factor_indices, cells_per_factor, replace=False)
        else:
            selected_indices = factor_indices
            logging.info(f"Factor {factor} has fewer than {cells_per_factor} cells ({len(factor_indices)} cells)")
        
        final_mask[selected_indices] = True
    
    logging.info(f"Selected {np.sum(final_mask)} cells total")
    return visium_data[final_mask]
    
        
    
    
    

@comprehensive_profile
def run_infercnvpy(input_path, output_path):
    
    os.makedirs(output_path, exist_ok=True)
    
    ### Visium HD data ###
    # Read the anndata object that has the ficture output cell types already added
    visium_data = ad.read_h5ad(input_path)

    # normalize + log transform
    sc.pp.normalize_total(visium_data, target_sum=1e4)
    sc.pp.log1p(visium_data)
    
    ### ADD CHROMOSOME INFO ###
    visium_data = add_chromosome_info(visium_data)


    ### CELL ANNOTATION as Normal/Tumor ###
    # Add normal or tumor annotation in the obs dataframe
    # if Erbb2 is expressed, it is a tumor cell
    # if Foxa1 is expressed, is is a "normal cell"

    visium_data = annotate_normal(visium_data)

    
    logging.info("Running infercnv globally")
    # perform infercnv
    cnv.tl.infercnv(visium_data, reference_key="cell_label", window_size=100, reference_cat="normal", step=10,
        # n_jobs-> default uses all cores
    )
    logging.info("Finished running infercnvpy")
    
    # clustering by CNV profiles 
    # logging.info("Performing PCA")
    # cnv.tl.pca(visium_data, svd_solver="arpack", n_comps = 30)
    # logging.info("Finished PCA")
    
    # logging.info("Performing neighbors with BBKNN")
    # # visium_data.obs["batch"] = np.random.randint(0, 2, size=visium_data.n_obs)
    # # bbknn.bbknn(visium_data, use_rep="X_cnv_pca")
    
    # # # Transfer BBKNN results to the format expected by infercnvpy
    # # visium_data.uns['cnv_neighbors'] = {}
    # # visium_data.uns['cnv_neighbors']['params'] = {
    # #     'n_neighbors': 15,
    # #     'method': 'bbknn',
    # #     'metric': 'euclidean',
    # #     'use_rep': 'X_cnv_pca'
    # # }
    # # visium_data.uns['cnv_neighbors']['connectivities_key'] = 'connectivities'
    # # visium_data.uns['cnv_neighbors']['distances_key'] = 'distances'
    # cnv.pp.neighbors(visium_data)
    # logging.info("Finished neighbors")
    
    # logging.info("Performing Leiden clustering")
    # cnv.tl.leiden(visium_data, resolution=0.2, n_iterations=2, flavor="igraph", directed=False)
    # logging.info("Finished Leiden clustering")
    
    # # take a representative subset of the data (include each factor and get repr. cnv scores)
    # sampled_data = sample_data(visium_data)
    # sampled_data = visium_data

    # logging.info("Performing UMAP")
    # cnv.tl.umap(sampled_data, maxiter=200)
    # logging.info("Finished UMAP")
    # logging.info("Performing cnv score")
    # cnv.tl.cnv_score(sampled_data)
    # logging.info("Finished cnv score")
    
    # Create and save UMAP plots
    # logging.info("Creating UMAP plots")
    # fig = plt.figure(figsize=(11, 11))
    # ax1 = plt.subplot(221)
    # ax2 = plt.subplot(222)
    # ax3 = plt.subplot(223)
    # cnv.pl.umap(sampled_data, color="cnv_leiden", legend_loc="on data", legend_fontoutline=2, ax=ax1, show=False)
    # cnv.pl.umap(sampled_data, color="cnv_score", ax=ax2, show=False)
    # cnv.pl.umap(sampled_data, color="factor", ax=ax3, show=False)
    # fig.savefig(output_path + "/umap.png")
    # plt.close(fig)
    # logging.info("Finished UMAP plots")
    
    # # Create and save spatial plots
    # logging.info("Creating spatial plots")
    # fig = plt.figure(figsize=(11, 11))
    # ax1 = plt.subplot(121)
    # ax2 = plt.subplot(122)
    # sc.pl.spatial(sampled_data, color="cnv_leiden", size=10, ax=ax1, show=False, spot_size=1)
    # sc.pl.spatial(sampled_data, color="factor", size=10, ax=ax2, show=False, spot_size=1)
    # fig.savefig(output_path + "/spatial.png")
    # plt.close(fig)
    # logging.info("Finished spatial plots")
    
    logging.info("CNV value range", np.min(visium_data.obsm["X_cnv"].data), np.max(visium_data.obsm["X_cnv"].data))
    # clip the cnv values to be between -0.25 and 0.25
    visium_data.obsm["X_cnv"].data = np.clip(visium_data.obsm["X_cnv"].data, -0.25, 0.25)
    
    logging.info("Creating CNV plots")
    sc.tl.dendrogram(visium_data, groupby="factor")
    logging.info("Finished dendrogram")
    fig = cnv.pl.chromosome_heatmap(visium_data, groupby="factor", dendrogram=True, show=False)
    logging.info("Finished creating the fig")
    # save the heatmap
    for key, ax in fig.items():
        ax.figure.savefig(output_path + f"/{key}.png")
    plt.close()
    logging.info("Finished CNV plots")
    
    

    logging.info("Running infercnv on individual factors")
    # run_infercnvpy_on_individual_factors(visium_data, output_path)
    logging.info("Finished infercnv on individual factors")
    
    # save the anndata object
    logging.info("Saving the anndata object ...")
    visium_data.write_h5ad(output_path + "/infercnvpy_output.h5ad")


    return

        


if __name__ == "__main__":

    # use argparse to get the input path
    parser = argparse.ArgumentParser(description="Run infercnvpy on the ficture output")
    parser.add_argument("--input", type=str, help="Path to the anndata visium object")
    parser.add_argument("--output", type=str, help="Path to the infercnvpy output directory")
    args = parser.parse_args()

    if args.input is None or args.output is None:
        results_folder = "run_4" 
        results_folder = "8_um_nF_10"
        results_folder = "multiple_factors"
        args.input = f"results/{results_folder}/ficture/annotated_anndata.h5ad"
        args.output = f"results/{results_folder}/infercnvpy"
        args.output = "test_output_2um"
    
    
    #logging.info(f"Using {mp.cpu_count()} cores")



    run_infercnvpy(args.input, args.output)

