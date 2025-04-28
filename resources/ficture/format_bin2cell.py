import os
import json
import argparse
import subprocess

import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger
from scipy.io import mmwrite


def format_bin2cell_output(input_path: str, visium_binned_dir: str, output_path: str):
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
    adata = sc.read_h5ad(input_path)
    logger.info(f"Read bin2cell output: {input_path}")

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
    positions = adata.obs[["object_id", "array_row", "array_col"]]
    # create a column for in_tissue with all 1s
    positions["in_tissue"] = 1
    # reorder the columns
    positions = positions[["object_id", "in_tissue", "array_row", "array_col"]]
    # rename the columns
    positions.columns = ["barcode", "in_tissue", "array_row", "array_col"]
    # save the positions
    tissue_positions_raw = os.path.join(output_path, "tissue_positions.csv.gz")
    positions.to_csv(tissue_positions_raw, sep=",", index=False, compression="gzip")

    scale_factors_path = os.path.join(
        visium_binned_dir, "spatial", "scalefactors_json.json"
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
    parser.add_argument("--input", type=str, help="path to bin2cell output anndata object. ")
    parser.add_argument("--visium_binned_dir", type=str, help="path to Visium binned output. ")
    parser.add_argument("--output", type=str, help="Output directory")
    args = parser.parse_args()
    
    if args.input is None:
        args.input = "results/human_breast_cancer_final/bin2cell/bin2cell_output.h5ad"
        args.visium_binned_dir = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/binned_outputs/square_002um"
        args.output = "final_prior_test/2_um_post_b2c/b2c"
        os.makedirs(args.output, exist_ok=True)

    logger.info("Formatting bin2cell output for FICTURE")
    logger.info(f"Input: {args.input}")
    logger.info(f"Visium binned dir: {args.visium_binned_dir}")
    logger.info(f"Output: {args.output}")
    
    format_bin2cell_output(args.input, args.visium_binned_dir, args.output)
