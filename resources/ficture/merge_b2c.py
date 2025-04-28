""" 
Merge Bin2cell output with full ficture output. Use when running ficture and b2c in parallel.

"""
import argparse

import anndata as ad
import numpy as np
from scipy.spatial import cKDTree
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

from loguru import logger as logging


def merge_b2c(b2c_output:str, ficture_output:str, output_path:str, ficture_color_table:str=None):
    """
    Merge Bin2cell output with full ficture output. Use when running ficture and b2c in parallel.
    
    Simple approach, just assign the matching position or closest bins between b2c and ficture to the same celltype.
    
    More complex approach could be using a majority vote(assign all the bins that were aggregated a vote based on their ficture label).
    This requires the other b2c output as well.
    """
    # load the adata objects
    b2c_adata = ad.read_h5ad(b2c_output)
    logging.info(f"B2C adata:")
    logging.info(b2c_adata)
    ficture_adata = ad.read_h5ad(ficture_output)
    logging.info(f"Ficture adata:")
    logging.info(ficture_adata)
    
    # get the positions of the b2c bins
    b2c_positions = b2c_adata.obsm["spatial"]
    logging.info(f"B2C positions:")
    logging.info(b2c_positions[10,:].shape)
    
    # get the positions of the ficture bins
    ficture_positions = ficture_adata.obsm["spatial"]
    logging.info(f"Ficture positions:")
    logging.info(ficture_positions.shape)
    
    # check for direct matches, rounded to 2 decimal places
    b2c_positions_rounded = np.round(b2c_positions, 2)
    ficture_positions_rounded = np.round(ficture_positions, 2)

    # create a cKDTree for the b2c positions
    b2c_tree = cKDTree(ficture_positions_rounded)
    
    # find the nearest neighbor in the ficture positions
    distances, indices = b2c_tree.query(b2c_positions_rounded, k=1)
    logging.info(f"Max distance: {np.max(distances):.2f}")
    logging.info(f"Min distance: {np.min(distances):.2f}")
    logging.info(f"Mean distance: {np.mean(distances):.2f}")
    logging.info(f"Median distance: {np.median(distances):.2f}")
    logging.info(f"Std distance: {np.std(distances):.2f}")
        
    # plot the 10 cells with the max distance on spatial plot
    max_dist_indices = np.argsort(distances)[-10:]
    logging.info(f"Max distance indices: {max_dist_indices}")
    
    # Create a new observation annotation for highlighting max distance cells
    b2c_adata.obs['max_distance_cells'] = 'Other'
    # Use the actual observation names from the AnnData object
    max_dist_obs_names = b2c_adata.obs_names[max_dist_indices]
    b2c_adata.obs.loc[max_dist_obs_names, 'max_distance_cells'] = 'Max distance'
    
    # Plot using scanpy's spatial plot function
    sc.pl.spatial(b2c_adata, 
                 color='max_distance_cells',
                 title='Cells with Maximum Distance to Nearest Ficture Point',
                 size=10,
                 groups=['Max distance'],
                 save='_max_distance_cells.png')
    
    
    
    # merge the b2c adata with the ficture adata
    # do this by taking the label from the closest ficture bin
    b2c_adata.obs["factor"] = ficture_adata.obs["factor"].iloc[indices].values
    # turn into categorical values
    b2c_adata.obs["factor"] = pd.Categorical(b2c_adata.obs["factor"])
    
    if ficture_color_table is not None:
        # read the ficture_color_table
        ficture_color_table = pd.read_csv(ficture_color_table, sep="\t")
        # create a color mapping from the ficture_color_table which has R G B columns in decimal format
        color_mapping = dict(zip(ficture_color_table['Name'], ficture_color_table[['R', 'G', 'B']].apply(lambda x: tuple(x), axis=1)))
        
    sc.pl.spatial(b2c_adata, color="factor", title="Factor", size=1, save="merged_factor.png", palette=color_mapping)
    
    # plot a scatter plot on black background of the factor values based on spatial position using the color mapping
    plt.figure(figsize=(10, 10))
    plt.scatter(b2c_adata.obs["spatial"][:,0], b2c_adata.obs["spatial"][:,1], c=b2c_adata.obs["factor"].cat.codes, cmap=color_mapping)
    plt.savefig("merged_factor_scatter.png")
    
    # save the merged adata
    b2c_adata.write(output_path)
    
    return b2c_adata
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--b2c_output", type=str, required=False)
    parser.add_argument("--ficture_output", type=str, required=False)
    parser.add_argument("--output_path", type=str, required=False)
    parser.add_argument("--ficture_color_table", type=str, required=False)
    args = parser.parse_args()
    
    if args.b2c_output is None:
        args.b2c_output = "results/human_breast_cancer_final/bin2cell/bin2cell_output.h5ad"
        args.ficture_output = "final_prior_test/preproc_2um_major/postproc/postproc_ficture/ficture_anndata.h5ad"
        args.output_path = "final_prior_test/preproc_2um_major/postproc/merged_anndata.h5ad"
        args.ficture_color_table = "final_prior_test/preproc_2um_major/analysis/nF9.d_16/figure/nF9.d_16.rgb.tsv"
    
    merge_b2c(args.b2c_output, args.ficture_output, args.output_path, args.ficture_color_table)