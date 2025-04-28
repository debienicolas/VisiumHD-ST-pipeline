"""
Postprocessing of the Ficture output.

Choose the best amount of clusters and train width by using the coherence score.  

Merge the original visium HD data with the ficture output by adding an obs column with the cell type.

"""

import argparse
import gzip
import os
import sys
import json
from pathlib import Path
import subprocess
import anndata as ad

from loguru import logger
#from resources.utils.utils import load_initial_data, save_visium_data

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import squidpy as sq

def load_initial_data(input_path: str):
    """
    Load Visium data from 10x Genomics output, using the filtered feature barcode matrix
    """

    # check if there is a tissue_positions_list.csv file, if not, use the spatial folder parquet file to transform the coordinates
    if not os.path.exists(input_path + "/spatial/tissue_positions_list.csv"):
        command = f"parquet-tools csv {os.path.join(input_path, 'spatial', 'tissue_positions.parquet')} > {os.path.join(input_path, 'spatial', 'tissue_positions_list.csv')}"
        subprocess.run(command, shell=True)
        
    adata = sq.read.visium(path = input_path,
                            counts_file = "filtered_feature_bc_matrix.h5")
    
    # make the var names unique
    adata.var_names_make_unique()

    # some spatial coordinates are strings, convert them to floats
    adata.obsm["spatial"] = adata.obsm["spatial"].astype(float) # type: ignore
    
    return adata



# def postprocess_ficture(
#     input_path_ficture: str, input_path_visium: str, output_path: str
# ) -> None:

#     DISTANCE_THRESHOLD = 2

#     # load the visium data
#     visium_data = load_initial_data(input_path_visium)

#     # load the ficture metadata from the ficture output file
#     config = {}
#     with gzip.open(input_path_ficture, "rt") as f:
#         # First line: K and TOPK
#         line = f.readline().strip("#\n").split(";")
#         for param in line:
#             key, value = param.split("=")
#             config[key] = int(value)
#         # Second line: BLOCK parameters
#         line = f.readline().strip("#\n").split(";")
#         for param in line:
#             key, value = param.split("=")
#             config[key] = int(value) if value.isdigit() else value
#         # Third line: Image parameters
#         line = f.readline().strip("#\n").split(";")
#         for param in line:
#             key, value = param.split("=")
#             config[key] = float(value)

#     pixel_sorted = pd.read_csv(
#         input_path_ficture, sep="\t", compression="gzip", skiprows=4, header=None
#     )
#     pixel_sorted.columns = ["Block", "X", "Y", "K1", "K2", "K3", "P1", "P2", "P3"]

#     # scale the ficture output to the visium data scale
#     pixel_sorted["X"] = pixel_sorted["X"] / config["SCALE"] + config["OFFSET_X"]
#     pixel_sorted["Y"] = pixel_sorted["Y"] / config["SCALE"] + config["OFFSET_Y"]

#     # to convert to the tissue position in the original image, we need to use the scalefactors
#     microns_per_pixel = visium_data.uns["spatial"]["TG23-0227_VGEX_result_VISIUM"][
#         "scalefactors"
#     ]["microns_per_pixel"]
#     pixel_sorted["X"] = pixel_sorted["X"] / microns_per_pixel
#     pixel_sorted["Y"] = pixel_sorted["Y"] / microns_per_pixel

#     # take the most confident K for cell_type assignment
#     pixel_sorted["cell_type"] = pixel_sorted["K1"]

#     # drop the irrelevant columns
#     pixel_sorted = pixel_sorted.drop(
#         columns=["Block", "K1", "K2", "K3", "P1", "P2", "P3"]
#     )

#     # for some reason the X and Y column are interchanged in the ficture output
#     pixel_sorted = pixel_sorted.rename(columns={"X": "Y", "Y": "X"})

#     ### Assign cell types to the visium data
#     visium_positions_rounded = np.round(visium_data.obsm["spatial"])  # type: ignore
#     ficture_positions_rounded = np.round(pixel_sorted[["X", "Y"]].to_numpy())

#     tree = cKDTree(ficture_positions_rounded)
#     distances, indices = tree.query(visium_positions_rounded, k=1)

#     visium_data.obs["factor"] = pixel_sorted["cell_type"].iloc[indices].values

#     # distance threshold to remove outliers
#     visium_data.obs.loc[distances > DISTANCE_THRESHOLD, "cell_type"] = np.nan

#     save_visium_data(visium_data, output_path)

def postprocess_ficture_bin2cell(
    input_path_ficture: str,
    input_path_visium: str,
    output_path: str
) -> None:
    """
    Postprocess the Ficture output for the bin2cell data.
    """
    pass

def spatula_postprocess(
    ficture_output_analysis_folder: str,
    preprocessed_ficture_file: str,
    input_path_visium: str,
    output_folder: str,
    cluster_coherence_file: str,
    chosen_id: str,
) -> None:
    """
    Postprocess the Ficture output using Spatula's join-pixel-tsv method.
    Spatula join-pixel-tsv: (https://seqscope.github.io/spatula/tools/join_pixel_tsv/)

    The preprocessed ficture input is Y-sorted, so we also need to sort the ficture output by Y. (major axis)

    input_path_ficture: the path to the ficture output file: .pixel.sorted.tsv.gz
    preprocessed_ficture_file: the path to the preprocessed ficture file: .sorted.tsv.gz
    input_path_visium: the path to the visium data directory
    output_folder: the path the the ficture postprocessing output folder
    """
    logger.info("Starting spatula postprocessing")
    logger.info(f"Ficture output analysis folder: {ficture_output_analysis_folder}")
    logger.info(f"Preprocessed ficture file: {preprocessed_ficture_file}")
    logger.info(f"Input path visium: {input_path_visium}")
    logger.info(f"Output folder: {output_folder}")
    logger.info(f"Cluster coherence file: {cluster_coherence_file}")
    logger.info(f"Chosen ID file: {chosen_id}")

    #### Must first decide on which n_factors and train_width to use
    # if there is only a single dir in the ficture_output_analysis_folder, then use that
    if len(os.listdir(ficture_output_analysis_folder)) == 1:
        id = os.listdir(ficture_output_analysis_folder)[0]
    else:
        # read the cluster coherence file
        logger.info("Reading cluster coherence file")
        cluster_coh = pd.read_csv(cluster_coherence_file, sep="\t")
        # currently only have the perplexity score so just sort by that
        cluster_coh = cluster_coh.sort_values(by="perplexity", ascending=True)
        # get the first row id
        id = cluster_coh.iloc[0]["id"]

    # now build the output_file_ficture path
    output_file_ficture = Path(ficture_output_analysis_folder) / id
    # find the .pixel.sorted.tsv.gz file in the ficture_output_analysis_folder
    logger.info(f"Looking for pixel sorted file in {output_file_ficture}")
    pixel_sorted_file = list(Path(output_file_ficture).glob("*.pixel.sorted.tsv.gz"))
    assert (
        len(pixel_sorted_file) == 1
    ), "Found multiple .pixel.sorted.tsv.gz files in the ficture_output_analysis_folder"
    output_file_ficture = pixel_sorted_file[0]
    output_file_ficture = str(output_file_ficture)
    logger.info(f"Found pixel sorted file: {output_file_ficture}")

    logger.info(f"Best n_factors and train_width combination: {id}")

    ### Now the postprocessing can be done using the best n_factors and train_width combination

    # sort the ficture output file by the major axis (Y)
    # the new path is the same as the input path, but remove the suffix.sorted.tsv.gz and add the suffix sorted_by_major_axis.tsv.gz
    sorted_ficture_path = output_file_ficture.replace(
        ".sorted.tsv.gz", ".sorted_by_major_axis.tsv.gz"
    )
    logger.info(f"Sorting ficture output by major axis (Y)")
    cmd = f"(gzip -cd {output_file_ficture} | head | grep ^#; gzip -cd {output_file_ficture} | grep -v ^# | sort -S 1G -gk3) | gzip -c > {sorted_ficture_path}"
    logger.debug(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True)

    logger.info(f"Sorted ficture path: {sorted_ficture_path}")
    # first gunzip the sorted_ficture_path file
    logger.info("Removing spaces from the third line of the sorted ficture file")
    with gzip.open(sorted_ficture_path, "rt") as file:
        lines = file.readlines()
        lines[2] = lines[2].replace(" ", "")
    with gzip.open(sorted_ficture_path, "wt") as file:
        file.writelines(lines)

    logger.info("Removed spaces from the third line")

    scale_factors = input_path_visium + "/spatial/scalefactors_json.json"
    logger.info(f"Reading scale factors from {scale_factors}")
    pixels_per_um = json.load(open(scale_factors))["microns_per_pixel"]
    logger.info(f"Pixels per um: {pixels_per_um}")
    
    os.makedirs(output_folder, exist_ok=True)
    logger.info(f"Created output folder: {output_folder}")
    
    output_file_path = output_folder + "/transcripts_ficture_joined"

    ## join the files together
    cmd = f"resources/ficture/spatula/bin/spatula join-pixel-tsv \
            --mol-tsv {preprocessed_ficture_file} \
            --pix-prefix-tsv nF15__,{sorted_ficture_path} \
            --out-prefix {output_file_path}\
            --sort-axis Y"
    
    if not os.path.exists(output_file_path + ".tsv.gz"):
        logger.info("Joining files with spatula")
        logger.debug(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True)
        logger.info("Files joined successfully")
    else:
        logger.info(f"Output file already exists: {output_file_path}.tsv.gz")

    # this creates two files in the output folder:
    # - dist.hist.tsv
    # - tsv.gz -> has a factor number column and prob./confidence column

    # load the initial visium data
    logger.info(f"Loading initial visium data from {input_path_visium}")
    visium_data = load_initial_data(input_path_visium)
    logger.info(f"Visium data loaded with shape: {visium_data.shape}")
    logger.info(visium_data)

    # load the joined ficture data
    logger.info(f"Loading joined ficture data from {output_file_path}.tsv.gz")
    joined_ficture_data = pd.read_csv(
        output_file_path + ".tsv.gz", sep="\t", compression="gzip"
    )
    logger.info(f"Joined ficture data loaded with shape: {joined_ficture_data.shape}")
    logger.info(joined_ficture_data.head())
    # set the column name for the factor number
    factor_col = joined_ficture_data.columns[
        joined_ficture_data.columns.str.endswith("K1")
    ][0]
    logger.info(f"Factor column: {factor_col}")
    
    logger.info("Converting coordinates to pixel space")
    joined_ficture_data["X"] = joined_ficture_data["X"] / pixels_per_um
    joined_ficture_data["Y"] = joined_ficture_data["Y"] / pixels_per_um
    
    logger.info("Grouping by X and Y coordinates")
    joined_ficture_data = (
        joined_ficture_data.groupby(["X", "Y"])
        .agg({"gene": "first", factor_col: "first", "nF15__P1": "first"})
        .reset_index()
    )
    logger.info(f"After grouping, shape: {joined_ficture_data.shape}")

    logger.info("Preparing spatial data for KD-tree")
    visium_pos = visium_data.obsm["spatial"].copy()  # type: ignore
    # sort by the y column from smallest to largest and then by x column from smallest to largest
    sorted_indices = np.lexsort((visium_pos[:, 1], visium_pos[:, 0]))
    visium_pos = visium_pos[sorted_indices]

    # sort the ficture data by the y column from smallest to largest and then by x column from smallest to largest
    logger.info("Sorting ficture data by Y and X")
    joined_ficture_data = joined_ficture_data.sort_values(by=["Y", "X"])
    ficture_pos = joined_ficture_data[["Y", "X"]].to_numpy()
    
    logger.info(f"Visium pos max and min x and y: {np.max(visium_pos[:, 0])} {np.min(visium_pos[:, 0])} {np.max(visium_pos[:, 1])} {np.min(visium_pos[:, 1])}")
    logger.info(f"Visium array col max and min: {np.max(visium_data.obs['array_col'])} {np.min(visium_data.obs['array_col'])}")
    logger.info(f"Visium array row max and min: {np.max(visium_data.obs['array_row'])} {np.min(visium_data.obs['array_row'])}")
    logger.info(f"Ficture pos max and min x and y: {np.max(joined_ficture_data['X'])} {np.min(joined_ficture_data['X'])} {np.max(joined_ficture_data['Y'])} {np.min(joined_ficture_data['Y'])}")
    

    logger.info("Building KD-tree for nearest neighbor search")
    tree = cKDTree(ficture_pos)
    
    logger.info("Finding nearest neighbors for visium spots")
    distances, indices = tree.query(visium_pos, k=1)
    
    logger.info("Assigning factors to visium data")
    factors_array = np.full((len(visium_pos)), np.nan)
    valid_indices = distances < 1
    logger.info(f"Found {np.sum(valid_indices)} valid matches out of {len(visium_pos)} spots")
    
    factors_array[valid_indices] = joined_ficture_data[factor_col].values[
        indices[valid_indices]
    ]

    visium_data.obs["factor"] = np.full((len(visium_data.obs)), np.nan)
    visium_data.obs["factor"] = factors_array[
        np.argsort(sorted_indices)
    ]  # Restore original order

    # remove the rows with NA in the factor column
    logger.info(f"Before filtering NA factors, shape: {visium_data.shape}")
    visium_data = visium_data[visium_data.obs["factor"].notna()]  # type: ignore
    logger.info(f"After filtering NA factors, shape: {visium_data.shape}")

    # save_visium_data(visium_data, output_folder + "/ficture_anndata.pkl")

    logger.info("Preparing data for saving")
    visium_data.obs["in_tissue"] = visium_data.obs["in_tissue"].astype(bool)
    visium_data.obs["array_row"] = visium_data.obs["array_row"].astype(int)
    visium_data.obs["array_col"] = visium_data.obs["array_col"].astype(int)
    
    logger.info(f"Saving anndata to {output_folder}/ficture_anndata.h5ad")
    visium_data.write(output_folder + "/ficture_anndata.h5ad")  # type: ignore

    # write the chosen id to the file
    logger.info(f"Writing chosen ID '{id}' to {chosen_id}")
    with open(chosen_id, "w") as file:
        file.write(id)

    logger.info("Spatula postprocessing completed successfully")
    return

def bin2cell_postprocess(
    ficture_output_analysis_folder: str,
    preprocessed_ficture_file: str,
    input_path_visium: str,
    output_folder: str,
    cluster_coherence_file: str,
    chosen_id: str,
    scale_factors: str,
) -> None:
    """
    Postprocess the Ficture output using Spatula's join-pixel-tsv method.
    Spatula join-pixel-tsv: (https://seqscope.github.io/spatula/tools/join_pixel_tsv/)

    The preprocessed ficture input is Y-sorted, so we also need to sort the ficture output by Y. (major axis)

    input_path_ficture: the path to the ficture output file: .pixel.sorted.tsv.gz
    preprocessed_ficture_file: the path to the preprocessed ficture file: .sorted.tsv.gz
    input_path_visium: the path to the visium data directory
    output_folder: the path the the ficture postprocessing output folder
    """
    logger.info("Starting bin2cell postprocessing")
    logger.info(f"Ficture output analysis folder: {ficture_output_analysis_folder}")
    logger.info(f"Preprocessed ficture file: {preprocessed_ficture_file}")
    logger.info(f"Input path visium: {input_path_visium}")
    logger.info(f"Output folder: {output_folder}")
    logger.info(f"Cluster coherence file: {cluster_coherence_file}")
    logger.info(f"Chosen ID file: {chosen_id}")
    
    # visium_data = ad.read_h5ad(input_path_visium)
    # logger.info(visium_data['obsm'])
    # logger.info(visium_data['uns'])

    #### Must first decide on which n_factors and train_width to use
    # if there is only a single dir in the ficture_output_analysis_folder, then use that
    if len(os.listdir(ficture_output_analysis_folder)) == 1:
        id = os.listdir(ficture_output_analysis_folder)[0]
    else:
        #raise ValueError("There is more than one directory in the ficture_output_analysis_folder")

        # read the cluster coherence file
        logger.info("Reading cluster coherence file")
        cluster_coh = pd.read_csv(cluster_coherence_file, sep="\t")
        # currently only have the perplexity score so just sort by that
        cluster_coh = cluster_coh.sort_values(by="perplexity", ascending=True)
        # get the id from the chosen_id file
        with open(chosen_id, "r") as file:
            id = file.read()

    # now build the output_file_ficture path
    output_file_ficture = Path(ficture_output_analysis_folder) / id
    # find the .pixel.sorted.tsv.gz file in the ficture_output_analysis_folder
    logger.info(f"Looking for pixel sorted file in {output_file_ficture}")
    pixel_sorted_file = list(Path(output_file_ficture).glob("*.pixel.sorted.tsv.gz"))
    assert (
        len(pixel_sorted_file) == 1
    ), "Found multiple .pixel.sorted.tsv.gz files in the ficture_output_analysis_folder"
    output_file_ficture = pixel_sorted_file[0]
    output_file_ficture = str(output_file_ficture)
    logger.info(f"Found pixel sorted file: {output_file_ficture}")

    logger.info(f"Best n_factors and train_width combination: {id}")

    ### Now the postprocessing can be done using the best n_factors and train_width combination

    # sort the ficture output file by the major axis (Y)
    # the new path is the same as the input path, but remove the suffix.sorted.tsv.gz and add the suffix sorted_by_major_axis.tsv.gz
    sorted_ficture_path = output_file_ficture.replace(
        ".sorted.tsv.gz", ".sorted_by_major_axis.tsv.gz"
    )
    logger.info(f"Sorting ficture output by major axis (Y)")
    cmd = f"(gzip -cd {output_file_ficture} | head | grep ^#; gzip -cd {output_file_ficture} | grep -v ^# | sort -S 1G -gk3) | gzip -c > {sorted_ficture_path}"
    logger.debug(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True)

    logger.info(f"Sorted ficture path: {sorted_ficture_path}")
    # first gunzip the sorted_ficture_path file
    logger.info("Removing spaces from the third line of the sorted ficture file")
    with gzip.open(sorted_ficture_path, "rt") as file:
        lines = file.readlines()
        lines[2] = lines[2].replace(" ", "")
    with gzip.open(sorted_ficture_path, "wt") as file:
        file.writelines(lines)

    logger.info("Removed spaces from the third line")

    scale_factors = json.load(open(scale_factors))["microns_per_pixel"]
    logger.info(f"Reading scale factors from {scale_factors}")
    pixels_per_um = scale_factors
    logger.info(f"Pixels per um: {pixels_per_um}")
    
    os.makedirs(output_folder, exist_ok=True)
    logger.info(f"Created output folder: {output_folder}")
    
    output_file_path = output_folder + "/transcripts_ficture_joined"

    ## join the files together
    cmd = f"resources/ficture/spatula/bin/spatula join-pixel-tsv \
            --mol-tsv {preprocessed_ficture_file} \
            --pix-prefix-tsv nF15__,{sorted_ficture_path} \
            --out-prefix {output_file_path}\
            --sort-axis Y"
    
    if not os.path.exists(output_file_path + ".tsv.gz"):
        logger.info("Joining files with spatula")
        logger.debug(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True)
        logger.info("Files joined successfully")
    else:
        logger.info(f"Output file already exists: {output_file_path}.tsv.gz")

    # this creates two files in the output folder:
    # - dist.hist.tsv
    # - tsv.gz -> has a factor number column and prob./confidence column

    # load the initial visium data
    logger.info(f"Loading initial visium data from {input_path_visium}")
    visium_data = ad.read_h5ad(input_path_visium)
    logger.info(f"Visium data loaded with shape: {visium_data.shape}")
    print(visium_data)

    # load the joined ficture data
    logger.info(f"Loading joined ficture data from {output_file_path}.tsv.gz")
    joined_ficture_data = pd.read_csv(
        output_file_path + ".tsv.gz", sep="\t", compression="gzip"
    )
    logger.info(f"Joined ficture data loaded with shape: {joined_ficture_data.shape}")
    
    # set the column name for the factor number
    factor_col = joined_ficture_data.columns[
        joined_ficture_data.columns.str.endswith("K1")
    ][0]
    logger.info(f"Factor column: {factor_col}")
    
    logger.info("Converting coordinates to pixel space")
    joined_ficture_data["X"] = joined_ficture_data["X"] / pixels_per_um
    joined_ficture_data["Y"] = joined_ficture_data["Y"] / pixels_per_um
    
    logger.info("Grouping by X and Y coordinates")
    joined_ficture_data = (
        joined_ficture_data.groupby(["X", "Y"])
        .agg({"gene": "first", factor_col: "first", "nF15__P1": "first"})
        .reset_index()
    )
    logger.info(f"After grouping, shape: {joined_ficture_data.shape}")

    logger.info("Preparing spatial data for KD-tree")
    # check if need to use obsm spatial or obs array_row and array_col by checking the max and min of both with joined_ficture_data
    # check if the max and min of the joined_ficture_data is roughly the same as the max and min of the visium_data 
    if round(np.max(joined_ficture_data['X'])) == round(np.max(visium_data.obs['array_col'])) and round(np.min(joined_ficture_data['X'])) == round(np.min(visium_data.obs['array_col'])) and round(np.max(joined_ficture_data['Y'])) == round(np.max(visium_data.obs['array_row'])) and round(np.min(joined_ficture_data['Y'])) == round(np.min(visium_data.obs['array_row'])):
        logger.info("Using obsm spatial")
        visium_pos = visium_data.obsm["spatial"].copy()  # type: ignore
    else:
        logger.info("Using obs array_row and array_col")
        visium_pos = np.column_stack((visium_data.obs['array_col'], visium_data.obs['array_row']))
    
    logger.info(f"Visium pos max and min x and y: {np.max(visium_pos[:, 0])} {np.min(visium_pos[:, 0])} {np.max(visium_pos[:, 1])} {np.min(visium_pos[:, 1])}")
    logger.info(f"Ficture pos max and min x and y: {np.max(joined_ficture_data['X'])} {np.min(joined_ficture_data['X'])} {np.max(joined_ficture_data['Y'])} {np.min(joined_ficture_data['Y'])}")
    
    # sort by the y column from smallest to largest and then by x column from smallest to largest
    sorted_indices = np.lexsort((visium_pos[:, 1], visium_pos[:, 0]))
    visium_pos = visium_pos[sorted_indices]

    # sort the ficture data by the y column from smallest to largest and then by x column from smallest to largest
    logger.info("Sorting ficture data by Y and X")
    joined_ficture_data = joined_ficture_data.sort_values(by=["Y", "X"])
    ficture_pos = joined_ficture_data[["Y", "X"]].to_numpy()

    logger.info("Building KD-tree for nearest neighbor search")
    tree = cKDTree(ficture_pos)
    
    logger.info("Finding nearest neighbors for visium spots")
    distances, indices = tree.query(visium_pos, k=1)
    
    logger.info("Assigning factors to visium data")
    factors_array = np.full((len(visium_pos)), np.nan)
    valid_indices = distances < 1
    logger.info(f"Found {np.sum(valid_indices)} valid matches out of {len(visium_pos)} spots")
    
    factors_array[valid_indices] = joined_ficture_data[factor_col].values[
        indices[valid_indices]
    ]

    visium_data.obs["factor"] = np.full((len(visium_data.obs)), np.nan)
    visium_data.obs["factor"] = factors_array[
        np.argsort(sorted_indices)
    ]  # Restore original order

    # remove the rows with NA in the factor column
    logger.info(f"Before filtering NA factors, shape: {visium_data.shape}")
    visium_data = visium_data[visium_data.obs["factor"].notna()]  # type: ignore
    logger.info(f"After filtering NA factors, shape: {visium_data.shape}")

    # save_visium_data(visium_data, output_folder + "/ficture_anndata.pkl")

    logger.info("Preparing data for saving")
    # check if the visium_data has these columns
    if "in_tissue" not in visium_data.obs.columns:
        visium_data.obs["in_tissue"] = True
    if "array_row" not in visium_data.obs.columns:
        visium_data.obs["array_row"] = visium_data.obs["array_row"].astype(int)
    if "array_col" not in visium_data.obs.columns:
        visium_data.obs["array_col"] = visium_data.obs["array_col"].astype(int)
    
    logger.info(f"Saving anndata to {output_folder}/ficture_anndata.h5ad")
    visium_data.write(output_folder + "/ficture_anndata.h5ad")  # type: ignore

    # write the chosen id to the file
    logger.info(f"Writing chosen ID '{id}' to {chosen_id}")
    with open(chosen_id, "w") as file:
        file.write(id)

    logger.info("Spatula postprocessing completed successfully")
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Postprocess Ficture output")
    parser.add_argument("--ficture_output_analysis_folder",type=str,required=False,help="Path to Ficture output analysis folder")
    parser.add_argument("--preprocessed_ficture_file",type=str,required=False,help="Path to preprocessed Ficture file")
    parser.add_argument("--input_path_visium",type=str,required=False,help="Path to visium data directory")
    parser.add_argument("--output_folder",type=str,required=False,help="Path to path where the anndata object is/will be saved")
    parser.add_argument("--cluster_coherence_file",type=str,required=False,help="Path to the cluster coherence file")
    parser.add_argument("--chosen_id", type=str, required=False, help="Path to the chosen id file")
    parser.add_argument("--n_factors",type=str,required=False,help="The factors you want to analyze down the line, multiple are comma separated")
    parser.add_argument("--scale_factors",type=str,required=False,help="Path to the scale factors file")
    args = parser.parse_args()

    output_name = "8_um_nF_10"
    if args.output_folder is None:
        results_dir = "final_prior_test/preproc_2um_major"
        output_dir = os.path.join(results_dir, "postproc")
        args.ficture_output_analysis_folder = os.path.join(results_dir, "analysis")
        args.preprocessed_ficture_file = os.path.join(results_dir, "transcripts.sorted.tsv.gz")
        args.input_path_visium = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/binned_outputs/square_002um"
        args.output_folder = os.path.join(output_dir, "postproc_ficture")
        args.cluster_coherence_file = os.path.join(results_dir, "coherence.tsv")
        args.chosen_id = os.path.join(results_dir, "chosen_id.txt")
        args.scale_factors = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/binned_outputs/square_002um/spatial/scalefactors_json.json"
    if args.input_path_visium.endswith(".h5ad"):
        bin2cell_postprocess(
            args.ficture_output_analysis_folder,
            args.preprocessed_ficture_file,
            args.input_path_visium,
            args.output_folder,
            args.cluster_coherence_file,
            args.chosen_id,
            args.scale_factors
        )
    else:
        spatula_postprocess(
            args.ficture_output_analysis_folder,
            args.preprocessed_ficture_file,
            args.input_path_visium,
            args.output_folder,
            args.cluster_coherence_file,
            args.chosen_id
        )
        