"""
1. find microns_per_pixel
2. combine info from multiple files and create a single input file, sorted by a single axis
3. coordinates of input file is in pixels, but keep track of min and mac in microns
4. set plot_um_per_pixel to 2 when visualising output (ficture plot_pixel_full)
"""

import sys, os, re, copy, gzip, time, logging, pickle, argparse
import numpy as np
import pandas as pd
import subprocess
import json
import argparse
import scanpy as sc
from loguru import logger
from scipy.io import mmwrite


"""
Format raw Visium output data into correct input format for Ficture.
"""


def format_visium(input_path=None, output_path=None):

    if input_path is None or output_path is None:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--input",
            type=str,
            help="Input directory of certain Visium resolution binned output. ",
        )
        parser.add_argument("--output", type=str, help="Output directory")
        args = parser.parse_args()

        input_path = args.input
        output_path = args.output

    # if output directory does not exist, create it
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # set path variables
    tissue_positions_parq = os.path.join(
        input_path, "spatial", "tissue_positions.parquet"
    )
    tissue_positions_raw = os.path.join(output_path, "tissue_positions.csv.gz")
    scale_factors_path = os.path.join(input_path, "spatial", "scalefactors_json.json")
    microns_per_pixel = json.load(open(scale_factors_path))["microns_per_pixel"]

    # turn a parquet into a csv
    conversion_cmd = (
        f"parquet-tools csv {tissue_positions_parq} | gzip -c > {tissue_positions_raw}"
    )
    subprocess.run(conversion_cmd, shell=True)

    logger.info(f"Converted parquet to csv: {tissue_positions_raw}")

    # run the spatula convert-sge command: https://seqscope.github.io/spatula/tools/convert_sge/#converting-10x-genomics-visium-hd-feature-barcode-matrix
    convert_sge_cmd = f"resources/ficture/spatula/bin/spatula convert-sge \
        --in-sge {os.path.join(input_path, 'filtered_feature_bc_matrix')} \
        --pos {tissue_positions_raw}\
        --units-per-um {1/microns_per_pixel} \
        --colnames-count Count \
        --out-tsv {output_path} \
        --icols-mtx 1"
    subprocess.run(convert_sge_cmd, shell=True)

    logger.info(f"Converted sge to tsv: {output_path}")

    # this creates a file called transcripts.unsorted.tsv.gz

    unsorted_file = os.path.join(output_path, "transcripts.unsorted.tsv.gz")
    sorted_file = os.path.join(output_path, "transcripts.sorted.tsv.gz")

    logger.info(f"Sorting unsorted transcripts file")
    # Update the sort command to use the sampled file
    sort_cmd = f"""(gzip -cd {unsorted_file} \
        | head -1; gzip -cd {unsorted_file} \
        | tail -n +2 | sort -S 1G -gk2) \
        | gzip -c > {sorted_file}"""
    subprocess.run(sort_cmd, shell=True)

    logger.info(f"Formatted visium data for FICTURE")

    return

def format_anndata(input_adata, input_path, output_path):
    """
    Format the output from bin2cell to be used as input for FICTURE.

    input requirements for spatula:
    - filtered_feature_bc_matrix:
        - barcodes.tsv.gz: single column with barcode names
        - features.tsv.gz: gene id, gene name , "Gene expression" column
        - matrix.mtx.gz: count matrix
    - tissue_positions.csv.gz:

    """
    os.makedirs(output_path, exist_ok=True)

    # read the bin2cell output
    adata = sc.read_h5ad(input_adata)
    logger.info(f"Read bin2cell output: {input_adata}")
    logger.info(adata)
    # create the filtered_feature_bc_matrix directory
    filtered_feature_bc_matrix_dir = os.path.join(
        output_path, "filtered_feature_bc_matrix"
    )
    os.makedirs(filtered_feature_bc_matrix_dir, exist_ok=True)

    # create the barcodes.tsv.gz file
    barcodes = adata.obs_names.to_frame()
    barcodes.to_csv(
        os.path.join(filtered_feature_bc_matrix_dir, "barcodes.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
        header=False,
    )

    # create the features.tsv.gz file
    features = adata.var
    logger.info(f"Features: {features.head()}")
    # drop the genome column and n_cells column
    features = features.drop(columns=["genome", "n_cells"])
    # create a gene_name column that is the index
    features["gene_name"] = features.index
    logger.info(f"Features: {features.head()}")
    # reorder the columns
    features = features[["gene_ids", "gene_name", "feature_types"]]
    features.to_csv(
        os.path.join(filtered_feature_bc_matrix_dir, "features.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
        header=False,
    )

    # round the matrix to the nearest integer
    adata.X = adata.X.astype(int)  # type: ignore
    # create the matrix.mtx.gz file
    matrix = adata.X.transpose()  # type:ignore
    # add metadata to the matrix
    metadata = "mtx file from bin2cell output"
    mmwrite(
        os.path.join(filtered_feature_bc_matrix_dir, "matrix.mtx"),
        matrix,
        comment=metadata,
    )
    subprocess.run(
        f"gzip {os.path.join(filtered_feature_bc_matrix_dir, 'matrix.mtx')}", shell=True
    )

    # create the tissue_positions.csv.gz file
    positions = adata.obs[["in_tissue", "array_row", "array_col"]]
    # create a column for in_tissue with all 1s
    positions["object_id"] = adata.obs.index
    # reorder the columns
    positions = positions[["object_id", "in_tissue", "array_row", "array_col"]]
    # rename the columns
    positions.columns = ["barcode", "in_tissue", "array_row", "array_col"]
    # save the positions
    tissue_positions_raw = os.path.join(output_path, "tissue_positions.csv.gz")
    positions.to_csv(tissue_positions_raw, sep=",", index=False, compression="gzip")

    scale_factors_path = os.path.join(
        input_path, "spatial", "scalefactors_json.json"
    )
    microns_per_pixel = json.load(open(scale_factors_path))["microns_per_pixel"]

    # run the spatula convert-sge command: https://seqscope.github.io/spatula/tools/convert_sge/#converting-10x-genomics-visium-hd-feature-barcode-matrix
    convert_sge_cmd = f"resources/ficture/spatula/bin/spatula convert-sge \
        --in-sge {filtered_feature_bc_matrix_dir} \
        --pos {tissue_positions_raw}\
        --units-per-um {1/microns_per_pixel} \
        --colnames-count Count \
        --out-tsv {output_path} \
        --icols-mtx 1 \
        --pos-colname-x array_row \
        --pos-colname-y array_col"
    subprocess.run(convert_sge_cmd, shell=True)

    unsorted_file = os.path.join(output_path, "transcripts.unsorted.tsv.gz")
    sorted_file = os.path.join(output_path, "transcripts.sorted.tsv.gz")

    logger.info(f"Sorting unsorted transcripts file")
    # Update the sort command to use the sampled file
    sort_cmd = f"""(gzip -cd {unsorted_file} \
        | head -1; gzip -cd {unsorted_file} \
        | tail -n +2 | sort -S 1G -gk2) \
        | gzip -c > {sorted_file}"""
    subprocess.run(sort_cmd, shell=True)

    logger.info(f"Formatted bin2cell output for FICTURE")

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",type=str,help="Input directory of certain Visium resolution binned output.", default=None)
    parser.add_argument("--input_adata",type=str,help="Input anndata object. ", default=None)
    parser.add_argument("--output", type=str, help="Output directory")
    args = parser.parse_args()
    logger.info(f"Formatting visium data for FICTURE")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # if input is an anndata object
    if args.input_adata is not None:
        format_anndata(args.input_adata,args.input, args.output)
    else:
        format_visium(args.input, args.output)
