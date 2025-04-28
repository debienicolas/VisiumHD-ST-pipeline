"""
    Use this script to perform after the fact bin2cell visualizations
    
"""
import os 

import anndata as ad
import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import bin2cell as b2c
import numpy as np

from loguru import logger as logging 

# set the scanpy plotting settings
sc.set_figure_params(figsize=(15, 15), dpi=1000)

if __name__ == "__main__":
    results_name = "human_breast_cancer_final"
    bin2cell_adata = f"results/{results_name}/bin2cell/bin2cell_unified_labels.h5ad"
    stardist_folder = f"results/{results_name}/bin2cell/stardist"
    output_path = f"results/{results_name}/bin2cell/bin2cell_vis"
    os.makedirs(output_path, exist_ok=True)
    
    alpha_img = 0.65
    
    # load the anndata object
    adata = ad.read_h5ad(bin2cell_adata)
    logging.info("Read the anndata object")
    logging.info(adata)
    
    # print the spatial coordinates max and min
    x_min, x_max = adata.obsm["spatial"][:, 0].min(), adata.obsm["spatial"][:, 0].max()
    y_min, y_max = adata.obsm["spatial"][:, 1].min(), adata.obsm["spatial"][:, 1].max()
    logging.info(f"Spatial coordinates range: x={x_min}-{x_max}, y={y_min}-{y_max}")
    
    # define a mask
    center_x = adata.obs['array_row'].median()
    center_y = adata.obs['array_col'].median()
    buffer = 5
    mask = ((adata.obs['array_row'] >= center_x - buffer) & 
        (adata.obs['array_row'] <= center_x + buffer) & 
        (adata.obs['array_col'] >= center_y - buffer) & 
        (adata.obs['array_col'] <= center_y + buffer)
       )
    bdata = adata[mask]
    min_x, max_x = bdata.obsm["spatial"][:, 0].min(), bdata.obsm["spatial"][:, 0].max()
    min_y, max_y = bdata.obsm["spatial"][:, 1].min(), bdata.obsm["spatial"][:, 1].max()
    logging.info(f"Spatial coordinates range: x={min_x}-{max_x}, y={min_y}-{max_y}")
    
    # spatial mask for good tissue
    # mask = ((adata.obsm["spatial"][:, 1] >= 2300) & 
    #         (adata.obsm["spatial"][:, 1] <= 2900) & 
    #         (adata.obsm["spatial"][:, 0] >= 15500) & 
    #         (adata.obsm["spatial"][:, 0] <= 16100))
    # bdata = adata[mask]
    # #bdata = adata
    # logging.info("Defined the mask")
    # logging.info(bdata)
    
    # print the spatial coordinates max and min
    int_x_min, int_x_max = int(bdata.obsm["spatial"][:, 0].min()), int(bdata.obsm["spatial"][:, 0].max())
    int_y_min, int_y_max = int(bdata.obsm["spatial"][:, 1].min()), int(bdata.obsm["spatial"][:, 1].max())
    logging.info(f"Spatial coordinates range: x={int_x_min}-{int_x_max}, y={int_y_min}-{int_y_max}")
    
    
    # plot the bins in a single color
    sc.pl.spatial(bdata, color="in_tissue", show=True, size=1)
    plt.savefig(os.path.join(output_path, "bins_in_tissue.png"))
    plt.close()
    logging.info("Plotted the bins in a single color")
    
   # plot the mask with bin grid lines
    # Extract the image from the AnnData object
    img = bdata.uns['spatial']['Visium_HD_Human_Breast_Cancer_Fresh_Frozen']['images']['hires']
    scale_factor = bdata.uns['spatial']['Visium_HD_Human_Breast_Cancer_Fresh_Frozen']['scalefactors']['tissue_hires_scalef']

    # Print debug information
    logging.info(f"Image shape: {img.shape}")
    logging.info(f"Scale factor: {scale_factor}")

    # Transform spatial coordinates to image coordinates
    img_x_min = int((int_x_min) * scale_factor)
    img_x_max = int((int_x_max) * scale_factor)
    img_y_min = int((int_y_min) * scale_factor)
    img_y_max = int((int_y_max) * scale_factor)

    logging.info(f"Image coordinates: x={img_x_min}-{img_x_max}, y={img_y_min}-{img_y_max}")

    # Ensure coordinates are within image bounds
    img_x_min = max(0, min(img_x_min, img.shape[1]))
    img_x_max = max(0, min(img_x_max, img.shape[1]))
    img_y_min = max(0, min(img_y_min, img.shape[0]))
    img_y_max = max(0, min(img_y_max, img.shape[0]))

    # Create a new figure and plot the image
    plt.figure(figsize=(15, 15))
    plt.imshow(img[img_y_min:img_y_max, img_x_min:img_x_max])

    # Get current axis
    ax = plt.gca()

    # Create grid lines using the plot's limits
    step = 2  # Adjust this value to match your bin size
    x_min, x_max = 0, img_x_max - img_x_min  # Use cropped image dimensions
    y_min, y_max = 0, img_y_max - img_y_min

    # Plot grid lines
    for x in np.arange(np.floor(x_min/step)*step, np.ceil(x_max/step)*step + step, step):
        plt.axvline(x, color="black", linewidth=0.9, alpha=1.0)
    for y in np.arange(np.floor(y_min/step)*step, np.ceil(y_max/step)*step + step, step):
        plt.axhline(y, color="black", linewidth=0.9, alpha=1.0)

    plt.axis('equal')
    plt.savefig(os.path.join(output_path, "mask_original.png"))
    plt.close()
    logging.info("Plotted the mask original")
    
    # plot the counts of the bins in the mask
    sc.pl.spatial(bdata, color="n_counts", color_map="OrRd",show=True)
    plt.savefig(os.path.join(output_path, "counts_mask.png"))
    plt.close()
    logging.info("Plotted the counts of the bins in the mask")
    
    # plot the adjusted counts of the bins in the mask
    sc.pl.spatial(bdata, color="n_counts_adjusted", color_map="OrRd",show=True)
    plt.savefig(os.path.join(output_path, "counts_mask_adjusted.png"))
    plt.close()
    logging.info("Plotted the adjusted counts of the bins in the mask")
    
    
    # plot the labels of the he segmentation
    bdata_he = bdata[bdata.obs["labels_he"] > 0]
    bdata_he.obs["labels_he"] = bdata_he.obs["labels_he"].astype(str).astype("category")
    
    # Create a color palette with one color per unique label
    n_labels = len(bdata_he.obs["labels_he"].unique())
    colors = plt.cm.tab20(np.linspace(0, 1, n_labels))  # using tab20 colormap for distinct colors
    color_dict = dict(zip(bdata_he.obs["labels_he"].unique(), colors))
    
    sc.pl.spatial(bdata_he, color="labels_he", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 palette=color_dict,
                 alpha_img=alpha_img, 
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_he.png"))
    plt.close()
    logging.info("Plotted the labels of the he segmentation")
    
    # Plot the discarded bins
    bdata_discarded = bdata[bdata.obs["labels_he"] <= 0]
    bdata_discarded.obs["labels_he"] = bdata_discarded.obs["labels_he"].astype(str).astype("category")
    
    sc.pl.spatial(bdata_discarded, color="labels_he", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_he_discarded.png"))
    plt.close()
    logging.info("Plotted the discarded bins")
    
    # plot the h&e segmentation masks
    crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped_150_buffer", mpp=0.5)
    rendered = b2c.view_labels(image_path=os.path.join(stardist_folder, "he.tiff"), 
                           labels_npz_path=os.path.join(stardist_folder, "he.npz"), 
                           crop=crop
                          )
    plt.imshow(rendered)
    plt.grid(False)
    plt.savefig(os.path.join(output_path, "labels_he_masks.png"))
    plt.close()
    logging.info("Plotted the h&e segmentation masks")
    
    
    
    ### GEX segmentation ### 
    
    
    bdata_gex = bdata[bdata.obs["labels_gex"] > 0]
    bdata_gex.obs["labels_gex"] = bdata_gex.obs["labels_gex"].astype(str).astype("category")
    
    # Create a color palette with one color per unique label
    n_labels = len(bdata_gex.obs["labels_gex"].unique())
    colors = plt.cm.tab20(np.linspace(0, 1, n_labels))  # using tab20 colormap for distinct colors
    color_dict = dict(zip(bdata_gex.obs["labels_gex"].unique(), colors))
    
    
    sc.pl.spatial(bdata_gex, color="labels_gex", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 palette=color_dict,
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_gex.png"))
    plt.close()
    logging.info("Plotted the labels of the gex segmentation")
    
    # plot the discarded bins
    bdata_discarded = bdata[bdata.obs["labels_gex"] <= 0]
    bdata_discarded.obs["labels_gex"] = bdata_discarded.obs["labels_gex"].astype(str).astype("category")
    
    sc.pl.spatial(bdata_discarded, color="labels_gex", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_gex_discarded.png"))
    plt.close()
    logging.info("Plotted the discarded bins")
    
    # plot the gex segmentation masks
    crop = b2c.get_crop(bdata, basis="array", mpp=0.5)
    rendered = b2c.view_labels(image_path=os.path.join(stardist_folder, "gex.tiff"), 
                           labels_npz_path=os.path.join(stardist_folder, "gex.npz"), 
                           crop=crop,
                           stardist_normalize=True
                          )
    plt.imshow(rendered)
    plt.grid(False)
    plt.savefig(os.path.join(output_path, "labels_gex_masks.png"))
    plt.close()
    logging.info("Plotted the gex segmentation masks")
    
    
    
    ### Combined plot - labels_joint ### 
    
    bdata_joint = bdata[bdata.obs["labels_joint"] > 0]
    bdata_joint.obs["labels_joint"] = bdata_joint.obs["labels_joint"].astype(str).astype("category")
    bdata_joint_discarded = bdata[bdata.obs["labels_joint"] <= 0]
    bdata_joint_discarded.obs["labels_joint"] = bdata_joint_discarded.obs["labels_joint"].astype(str).astype("category")
    
    sc.pl.spatial(bdata_joint, color="labels_joint_source", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_joint_source.png"))
    plt.close()
    logging.info("Plotted the labels of the joint segmentation")
    
    
    sc.pl.spatial(bdata_joint_discarded, color="labels_joint_source", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_joint_source_discarded.png"))
    plt.close()
    logging.info("Plotted the discarded bins")
    
    
    # create color palette for the labels
    n_labels = len(bdata_joint.obs["labels_joint"].unique())
    colors = plt.cm.tab20(np.linspace(0, 1, n_labels))  # using tab20 colormap for distinct colors
    color_dict = dict(zip(bdata_joint.obs["labels_joint"].unique(), colors))
    
    sc.pl.spatial(bdata_joint, color="labels_joint", 
                 img_key="0.5_mpp_150_buffer", 
                 basis="spatial_cropped_150_buffer", 
                 legend_loc=None,
                 palette=color_dict,
                 alpha_img=alpha_img,
                 show=True)
    plt.savefig(os.path.join(output_path, "labels_joint.png"))
    plt.close()
    logging.info("Plotted the labels of the joint segmentation")
    
    
    
    
    
    
    