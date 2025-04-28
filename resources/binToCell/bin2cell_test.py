import argparse, os 
import scanpy as sc
import bin2cell as b2c 
from loguru import logger as logging
import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger as logging

""" 
Can have a stardist model folder to load the model from
"""


def run_bin2cell(binned_input_dir:str, spaceranger_spatial_dir:str, highres_input_image:str, output_dir:str, 
                 prob_thresholds_gex=None, nms_thresholds_gex=None):
    
    logging.info("Starting bin2cell")
    logging.info("Input arguments:")
    logging.info(f"binned_input_dir: {binned_input_dir}")
    logging.info(f"spaceranger_spatial_dir: {spaceranger_spatial_dir}")
    logging.info(f"highres_input_image: {highres_input_image}")
    logging.info(f"output_dir: {output_dir}")
    
    # Use default thresholds if none provided
    if prob_thresholds_gex is None:
        prob_thresholds_gex = [0.05]  # Default value
    if nms_thresholds_gex is None:
        nms_thresholds_gex = [0.5]    # Default value
        
    logging.info(f"GEX probability thresholds: {prob_thresholds_gex}")
    logging.info(f"GEX NMS thresholds: {nms_thresholds_gex}")
    
    # check which device is being used
    device_name = tf.test.gpu_device_name()
    if device_name:
        logging.info(f"TensorFlow is using GPU: {device_name}")
    else:
        logging.info("TensorFlow is using CPU")
    
    
    ### Setting important parameters ###
    mpp = 0.5
    buffer = 150 # how many extra pixels to include 
    prob_thresh = 0.02 # for H&E segmentation
    max_bin_distance = 2 # for post H&E segmentation label expansion
    
    prob_thresh_gex = 0.05 # for GEX segmentation
    nms_thresh_gex = 0.5 # for GEX segmentation
    
    # set the stardist model directory
    stardist_model_dir = "resources/binToCell/stardist_model"
    
    # create stardist output directory and regular output directory
    os.makedirs(output_dir, exist_ok=True)
    stardist_output_dir = os.path.join(output_dir, "stardist")
    os.makedirs(stardist_output_dir, exist_ok=True)
    
    # read the data
    adata = b2c.read_visium(binned_input_dir, source_image_path = highres_input_image, spaceranger_image_path = spaceranger_spatial_dir)
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
    
    ### GEX segmentation with multiple thresholds ###
    logging.info("Running GEX segmentation with multiple thresholds...")
    
    # Create the base GEX image once
    gex_image_path = os.path.join(stardist_output_dir, "gex.tiff")
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path=gex_image_path)
    
    # Create a directory for threshold experiments
    thresholds_dir = os.path.join(output_dir, "threshold_experiments")
    os.makedirs(thresholds_dir, exist_ok=True)
    
    # Get crop for rendering
    crop_he = b2c.get_crop(adata, basis="spatial", mpp=mpp, spatial_key=f"spatial_cropped_{buffer}_buffer")
    rendered_he = b2c.view_labels(image_path=os.path.join(stardist_output_dir, "he.tiff"), 
                                        labels_npz_path=os.path.join(stardist_output_dir, "he.npz"), 
                                        crop=crop_he
                                    )
    # save rendered he to output directory
    np.save(os.path.join(thresholds_dir, "rendered_he.npy"), rendered_he)
    
    # Store the best segmentation results based on default thresholds
    best_labels_key = None
    best_prob_thresh = None
    best_nms_thresh = None
    
    # Run segmentation for each combination of thresholds
    for prob_thresh_gex in prob_thresholds_gex:
        for nms_thresh_gex in nms_thresholds_gex:
            # Create a unique name for this threshold combination
            thresh_name = f"prob{prob_thresh_gex}_nms{nms_thresh_gex}"
            labels_key = f"labels_gex_{thresh_name}"
            
            # Run stardist with these thresholds
            labels_npz_path = os.path.join(thresholds_dir, f"gex_{thresh_name}.npz")
            
            logging.info(f"Running GEX segmentation with prob_thresh={prob_thresh_gex}, nms_thresh={nms_thresh_gex}")
            b2c.stardist(image_path=gex_image_path, 
                         labels_npz_path=labels_npz_path, 
                         stardist_model="2D_versatile_fluo", 
                         prob_thresh=prob_thresh_gex, 
                         nms_thresh=nms_thresh_gex)
            
            b2c.insert_labels(adata, 
                          labels_npz_path=labels_npz_path, 
                          basis="array", 
                          mpp=mpp, 
                          labels_key=labels_key)
            
            crop_gex = b2c.get_crop(adata, basis="array", mpp=mpp)
            rendered_gex = b2c.view_labels(image_path=gex_image_path, 
                                    labels_npz_path=labels_npz_path, 
                                    crop=crop_gex
                                   )
            
            
            # Save the rendered image
            np.save(os.path.join(thresholds_dir, f"rendered_gex_{thresh_name}.npy"), rendered_gex)
            
            
            # plot the same mask on the 
            mask = ((adata.obs['array_row'] >= 1400) & 
                    (adata.obs['array_row'] <= 1700) & 
                    (adata.obs['array_col'] >= 400) & 
                    (adata.obs['array_col'] <= 700)
                )

            bdata = adata[mask]
            fig, axs = plt.subplots(1, 3, figsize=(10, 5))
            sc.pl.spatial(bdata, color=None, img_key=f"{mpp}_mpp_{buffer}_buffer", basis=f"spatial_cropped_{buffer}_buffer", show=False, ax=axs[0])
            axs[1].imshow(rendered_he)
            axs[2].imshow(rendered_gex)
            plt.savefig(os.path.join(thresholds_dir, f"plot_{thresh_name}.png"))
            plt.close()
            
            # If this is the default threshold combination, mark it as the best
            if prob_thresh_gex == 0.05 and nms_thresh_gex == 0.5:
                best_labels_key = labels_key
                best_prob_thresh = prob_thresh_gex
                best_nms_thresh = nms_thresh_gex
    
    logging.info("GEX segmentation with multiple thresholds complete")
    
    # Use the best segmentation for the rest of the pipeline
    logging.info(f"Using GEX segmentation with prob_thresh={best_prob_thresh}, nms_thresh={best_nms_thresh} for unification")
    
    
   
    
    logging.info("Finished running bin2cell")
    
    return None 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--um_binned_input_dir", type=str, required=False, default=None)
    parser.add_argument("--spaceranger_spatial_dir", type=str, required=False, default=None)
    parser.add_argument("--highres_input_image", type=str, required=False, default=None)
    parser.add_argument("--output_dir", type=str, required=False, default=None)
    parser.add_argument("--prob_thresholds_gex", type=str, required=False, default=None, 
                        help="Comma-separated list of probability thresholds for GEX segmentation")
    parser.add_argument("--nms_thresholds_gex", type=str, required=False, default=None,
                        help="Comma-separated list of NMS thresholds for GEX segmentation")
    args = parser.parse_args()

    if args.um_binned_input_dir is None:
        args.um_binned_input_dir = "input/TG23-0227_VGEX_results_VISIUM_HD_reanalyzed_with_highquality_HE/binned_outputs/square_002um"
        args.spaceranger_spatial_dir = "input/TG23-0227_VGEX_results_VISIUM_HD_reanalyzed_with_highquality_HE/spatial"
        args.highres_input_image = "input/TG23-0227_VGEX_results_VISIUM_HD_reanalyzed_with_highquality_HE/highres_input.tif"
        args.output_dir = "results/highres/bin2cell_threshold_experiments"
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Default threshold values to try if not specified
        args.prob_thresholds_gex = "0.03,0.05,0.07"
        args.nms_thresholds_gex = "0.3,0.5,0.7"
    
    # Parse threshold lists
    prob_thresholds_gex = [float(x) for x in args.prob_thresholds_gex.split(",")] if args.prob_thresholds_gex else None
    nms_thresholds_gex = [float(x) for x in args.nms_thresholds_gex.split(",")] if args.nms_thresholds_gex else None
    
    run_bin2cell(binned_input_dir=args.um_binned_input_dir,
                 spaceranger_spatial_dir=args.spaceranger_spatial_dir,
                 highres_input_image=args.highres_input_image,
                 output_dir=args.output_dir,
                 prob_thresholds_gex=prob_thresholds_gex,
                 nms_thresholds_gex=nms_thresholds_gex)
    
    

    