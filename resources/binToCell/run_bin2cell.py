import argparse, os 
import scanpy as sc
import bin2cell as b2c 
from loguru import logger as logging
import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger as logging
import json

""" 
Can have a stardist model folder to load the model from
"""


def run_bin2cell(input_adata:str, visium_dir:str, spaceranger_spatial_dir:str, highres_input_image:str, output_dir:str, config:dict=None):
    
    logging.info("Starting bin2cell")
    logging.info("Input arguments:")
    logging.info(f"input_adata: {input_adata}")
    logging.info(f"spaceranger_spatial_dir: {spaceranger_spatial_dir}")
    logging.info(f"highres_input_image: {highres_input_image}")
    logging.info(f"output_dir: {output_dir}")
    
    # check which device is being used
    device_name = tf.test.gpu_device_name()
    if device_name:
        logging.info(f"TensorFlow is using GPU: {device_name}")
    else:
        logging.info("TensorFlow is using CPU")
        
    ### Extract the params fron the config dict
    
    default_config = {
        "mpp": 0.5, 
        "buffer": 150,
        "prob_thresh": 0.02,
        "max_bin_distance": 2,
        "prob_thresh_gex": 0.05,
    }
    
    if config:
        config = default_config.update(config)
        
    
    ### Setting important parameters ###
    mpp = 0.5
    buffer = 150 # how many extra pixels to include 
    prob_thresh = 0.02 # for H&E segmentation
    max_bin_distance = 2 # for post H&E segmentation label expansion
    
    prob_thresh_gex = 0.05 # for GEX segmentation
    nms_thresh_gex = 0.5 # for GEX segmentation
        
    # create stardist output directory and regular output directory
    os.makedirs(output_dir, exist_ok=True)
    stardist_output_dir = os.path.join(output_dir, "stardist")
    os.makedirs(stardist_output_dir, exist_ok=True)
    
    # read the data
    original_adata = b2c.read_visium(visium_dir, source_image_path = highres_input_image, spaceranger_image_path = spaceranger_spatial_dir)
    adata = sc.read_h5ad(input_adata)
    # subset the original adata to the bins that are present in the input adata
    original_adata = original_adata[original_adata.obs_names.isin(adata.obs_names)]
    
    #### ---- #### using the normalized input data #### ---- ####
    # adata = original_adata
    adata.var_names_make_unique()
    
    logging.info(f"anndata object loaded:")
    logging.info(adata)
    
    # filter the data
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    logging.info(f"anndata object filtered:")
    logging.info(adata)
    
    # set the microns per pixel of the desired H&E image to create
    b2c.scaled_he_image(adata, mpp=mpp, save_path=os.path.join(stardist_output_dir, "he.tiff"), buffer=buffer)
    b2c.destripe(adata)
    
    
    ### H&E segmentation ###
    logging.info("Running H&E segmentation...")
    b2c.stardist(image_path=os.path.join(stardist_output_dir, "he.tiff"), 
                 labels_npz_path=os.path.join(stardist_output_dir, "he.npz"), 
                 stardist_model="2D_versatile_he", 
                 prob_thresh=prob_thresh)
    
    b2c.insert_labels(adata, 
                  labels_npz_path=os.path.join(stardist_output_dir, "he.npz"), 
                  basis="spatial", 
                  spatial_key="spatial_cropped_150_buffer",
                  mpp=mpp, 
                  labels_key="labels_he"
                 )
    # expand the labels based on a max_bin_distance=2 by default
    b2c.expand_labels(adata, 
                      labels_key="labels_he", 
                      expanded_labels_key="labels_he_expanded",
                      max_bin_distance=max_bin_distance)
    logging.info("H&E segmentation complete")
    
    ### GEX segmentation ###
    logging.info("Running GEX segmentation...")
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path=os.path.join(stardist_output_dir, "gex.tiff"))

    b2c.stardist(image_path=os.path.join(stardist_output_dir, "gex.tiff"), 
                 labels_npz_path=os.path.join(stardist_output_dir, "gex.npz"), 
                 stardist_model="2D_versatile_fluo", 
                 prob_thresh=prob_thresh_gex, 
             nms_thresh=nms_thresh_gex
            )
    b2c.insert_labels(adata, 
                  labels_npz_path=os.path.join(stardist_output_dir, "gex.npz"), 
                  basis="array", 
                  mpp=mpp, 
                  labels_key="labels_gex"
                 )
    logging.info("GEX segmentation complete")
    
    ### Unification of the H&E and GEX labels ###
    b2c.salvage_secondary_labels(adata, 
                             primary_label="labels_he_expanded", 
                             secondary_label="labels_gex", 
                             labels_key="labels_joint"
                            )
    
    # write the entire anndata object to a file
    adata.write_h5ad(os.path.join(output_dir, "bin2cell_unified_labels.h5ad"))
    logging.info("Finished writing the anndata object to a file")
    logging.info(adata)
    
    cdata = b2c.bin_to_cell(adata, labels_key="labels_joint", spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    # round the cdata.X to the nearest integer
    cdata.X = np.round(cdata.X)
    
    logging.info("Final anndata object:")
    logging.info(cdata)
    
    logging.info("Finished the unification of the H&E and GEX labels")
    
    ### Save the output ###
    cdata.write_h5ad(os.path.join(output_dir, "bin2cell_output.h5ad"), compression="gzip")
    logging.info("Finished saving the output")
    logging.info("The anndata object written to file: ")
    logging.info(cdata)
    
    
    ## save the rendered H&E and GEX segmentation as .npy files
    crop_he = b2c.get_crop(adata, basis="spatial", mpp = mpp, spatial_key=f"spatial_cropped_{buffer}_buffer")
    rendered_he = b2c.view_labels(image_path=os.path.join(stardist_output_dir, "he.tiff"), 
                                  labels_npz_path=os.path.join(stardist_output_dir, "he.npz"), 
                                  crop=crop_he)
    crop_gex = b2c.get_crop(adata, basis="array", mpp = mpp)        
    rendered_gex = b2c.view_labels(image_path=os.path.join(stardist_output_dir, "gex.tiff"), 
                                  labels_npz_path=os.path.join(stardist_output_dir, "gex.npz"), 
                                  stardist_normalize=True,
                                  crop=crop_gex)
    np.save(os.path.join(output_dir, "rendered_he.npy"), rendered_he)
    np.save(os.path.join(output_dir, "rendered_gex.npy"), rendered_gex)
    
    
    logging.info("Finished saving the rendered H&E and GEX segmentation as .npy files")
    
    
    logging.info("Finished running bin2cell")
    
    
    ### save the anndata object of all the bins that were discarded
    adata_discarded = adata[adata.obs["labels_joint"] <= 0]
    logging.info(f"anndata object of all the bins that were discarded:")
    logging.info(adata_discarded)
    
    adata_discarded.write_h5ad(os.path.join(output_dir, "bin2cell_discarded_bins.h5ad"), compression="gzip")
    logging.info("Finished saving the anndata object of all the bins that were discarded")
    
    return None 


def create_html_report(output_dir:str):
    """
    Create an html report of the bin2cell output

    Args:
        output_dir (str): The directory containing the bin2cell output, where the html report will be written to
    """
    # read the complete anndata object
    adata = sc.read_h5ad(os.path.join(output_dir, "bin2cell_unified_labels.h5ad"))
    logging.info(f"anndata object of all the bins:")
    logging.info(adata)
    # read the merged anndata object
    adata_merged = sc.read_h5ad(os.path.join(output_dir, "bin2cell_output.h5ad"))
    if "n_counts" not in adata_merged.obs.columns:
        adata_merged.obs["n_counts"] = np.sum(adata_merged.X, axis=1)
    if "n_genes" not in adata_merged.obs.columns:
        adata_merged.obs["n_genes"] = np.sum(adata_merged.X > 0, axis=1)
        
    logging.info(f"anndata object of merged the bins:")
    logging.info(adata_merged)
    
    # Create statistics summary table
    logging.info("Creating statistics summary table...")
    
    # Gather statistics from both anndata objects
    stats = {
        "Total bins processed": adata.n_obs,
        "Bins assigned to cells": int(np.sum(adata.obs["labels_joint"] > 0)),
        "Bins discarded": int(np.sum(adata.obs["labels_joint"] <= 0)),
        "Total cells identified": adata_merged.n_obs,
        "Average bins per cell": round(np.mean(adata_merged.obs["bin_count"]) if "bin_count" in adata_merged.obs.columns else float('nan'), 2),
        "Median bins per cell": round(np.median(adata_merged.obs["bin_count"]) if "bin_count" in adata_merged.obs.columns else float('nan'), 2),
        "Max bins per cell": int(np.max(adata_merged.obs["bin_count"]) if "bin_count" in adata_merged.obs.columns else float('nan')),
        "Average UMIs per cell": round(np.mean(adata_merged.obs["n_counts"]) if "n_counts" in adata_merged.obs.columns else float('nan'), 2),
        "Median UMIs per cell": round(np.median(adata_merged.obs["n_counts"]) if "n_counts" in adata_merged.obs.columns else float('nan'), 2),
        "Average genes per cell": round(np.mean(adata.obs["n_genes"]) if "n_genes" in adata.obs.columns else float('nan'), 2),
        "Median genes per cell": round(np.median(adata.obs["n_genes"]) if "n_genes" in adata.obs.columns else float('nan'), 2),
        #"Average UMIs per bin discarded": round(np.mean(adata[ad].obs["n_counts"]) if "n_counts" in adata_discarded.obs.columns else float('nan'), 2),
    }
    
    # Compute percent of bins used
    bins_used_percent = round((stats["Bins assigned to cells"] / stats["Total bins processed"]) * 100, 2)
    stats["Percent bins used"] = f"{bins_used_percent}%"
    
    # Create HTML table
    html_content = "<h2>bin2cell Statistics Summary</h2>\n"
    html_content += "<table border='1' style='border-collapse: collapse; width: 80%; margin: 20px auto;'>\n"
    html_content += "<tr><th style='padding: 8px; text-align: left;'>Metric</th><th style='padding: 8px; text-align: right;'>Value</th></tr>\n"
    
    for metric, value in stats.items():
        html_content += f"<tr><td style='padding: 8px;'>{metric}</td><td style='padding: 8px; text-align: right;'>{value}</td></tr>\n"
    
    html_content += "</table>\n"
    
    # Create visualizations
    logging.info("Creating summary visualizations...")
    
    # Create a pie chart of bins used vs discarded
    plt.figure(figsize=(10, 6))
    plt.pie([stats["Bins assigned to cells"], stats["Bins discarded"]], 
            labels=["Assigned to cells", "Discarded"], 
            autopct='%1.1f%%', 
            colors=['#66b3ff', '#ff9999'])
    plt.title('Distribution of Bins')
    plt.savefig(os.path.join(output_dir, "bins_distribution_pie.png"))
    
    # Add visualization to HTML report
    html_content += "<h2>Bin Distribution</h2>\n"
    html_content += f"<img src='bins_distribution_pie.png' style='display: block; margin: 0 auto; max-width: 600px;'>\n"
    
    # If n_bins is available, create a histogram of bins per cell
    if "n_bins" in adata.obs.columns:
        plt.figure(figsize=(10, 6))
        plt.hist(adata.obs["n_bins"], bins=30, color='#66b3ff', alpha=0.8)
        plt.xlabel('Number of Bins')
        plt.ylabel('Number of Cells')
        plt.title('Distribution of Bins per Cell')
        plt.grid(alpha=0.3)
        plt.savefig(os.path.join(output_dir, "bins_per_cell_hist.png"))
        
        html_content += "<h2>Bins per Cell Distribution</h2>\n"
        html_content += f"<img src='bins_per_cell_hist.png' style='display: block; margin: 0 auto; max-width: 600px;'>\n"
    
    # Save the HTML report
    with open(os.path.join(output_dir, "bin2cell_report.html"), "w") as f:
        f.write(html_content)
    
    logging.info(f"HTML report created at: {os.path.join(output_dir, 'bin2cell_report.html')}")
    logging.info("Finished creating the html report")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_adata", type=str, required=False, default=None)
    parser.add_argument("--visium_dir", type=str, required=False, default=None)
    parser.add_argument("--spaceranger_spatial_dir", type=str, required=False, default=None)
    parser.add_argument("--highres_input_image", type=str, required=False, default=None)
    parser.add_argument("--output_dir", type=str, required=False, default=None)
    parser.add_argument("--config", type=str, required=False, default=None)
    args = parser.parse_args()
    
    default_config = {
        "mpp": 0.5,
        "buffer": 150,
        "prob_thresh": 0.02,
        "max_bin_distance": 2,
        "prob_thresh_gex": 0.05,
    }
    
    config = json.loads(args.config)
    
    if args.input_adata is None:
        input_dir = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen"
        results_dir = "results/human_breast_cancer_final"
        args.input_adata = os.path.join(input_dir, "binned_outputs", "square_002um")
        args.spaceranger_spatial_dir = os.path.join(input_dir, "spatial")
        args.highres_input_image = os.path.join(input_dir, "highres_input.tif")
        args.output_dir = os.path.join(results_dir, "bin2cell")
        os.makedirs(args.output_dir, exist_ok=True)
        
    
    run_bin2cell(input_adata=args.input_adata,
                visium_dir=args.visium_dir,
                spaceranger_spatial_dir=args.spaceranger_spatial_dir,
                highres_input_image=args.highres_input_image,
                output_dir=args.output_dir,
                config=config)
    
    create_html_report(output_dir=args.output_dir)

    