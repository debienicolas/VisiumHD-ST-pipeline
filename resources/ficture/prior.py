"""
    This script contains the methods to run ficture using prior knowledge.
    Relevant documentation: https://seqscope.github.io/ficture/localrun/#optional-initializing-lda-model-from-pseudo-bulk-data
"""
import os
import re
import anndata as ad
import pandas as pd
from ficture.scripts.run_together import run_together
import argparse
import subprocess
from ficture.utils.minimake import minimake
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger as logging



from matplotlib import cm

def create_pseudo_bulk(input_adata:str, reference_adata:str, output_dir:str, cell_type_column:str="celltype_subset"):
    """
    Create a pseudo-bulk dataset from input data and reference data.
    This creates a tsv file with genes as rows and cell types as columns.
    Args:
        input_adata: path to the input adata object
        reference_adata: path to the reference adata object
        output_dir: path to the output directory
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    # load the input adata object
    input_adata = ad.read_h5ad(input_adata)
    logging.info(f"Input adata object loaded")
    logging.info(input_adata)
    # load the reference adata object
    reference_adata = ad.read_h5ad(reference_adata)
    logging.info(f"Reference adata object loaded")
    logging.info(reference_adata)
    
    # get the genes that are present in both adata objects
    common_genes = list(set(input_adata.var_names) & set(reference_adata.var_names))
    logging.info(f"Amount of common genes: {len(common_genes)}")
    
    reference_adata = reference_adata[:, common_genes].copy()
    
    cell_types = reference_adata.obs[cell_type_column].unique().tolist()
    logging.info(f"Amount of cell types: {len(cell_types)}")
    # create a dataframe with the first column called gene and the rest of the columns called the cell types
    pseudo_bulk_df = pd.DataFrame(0,index=common_genes, columns=cell_types)
    pseudo_bulk_df.index.name = "gene"
    
    # the reference adata object is already normalized to 10_000
    # average the counts across the cell types
    for cell_type in reference_adata.obs[cell_type_column].unique():
        logging.info(f"Processing cell type: {cell_type}")
        cell_type_adata = reference_adata[reference_adata.obs[cell_type_column] == cell_type, :].copy()
        cell_type_adata = cell_type_adata[:, common_genes].copy()
        counts_mean = cell_type_adata.X.toarray().mean(axis=0)
        logging.info(f"Counts mean: {counts_mean.shape}")
        # add the counts to the pseudo_bulk_df
        pseudo_bulk_df[cell_type] = counts_mean
        
    
    # ficture requires the names to be of format name_1, name_2, name_3, etc.
    pseudo_bulk_df.columns = [f"{col}_{i}" for i,col in enumerate(pseudo_bulk_df.columns)]
    
    # save the pseudo_bulk_df
    pseudo_bulk_df.to_csv(os.path.join(output_dir, f"pseudo_bulk_{cell_type_column}.tsv"), sep="\t")
    logging.info(f"Pseudo-bulk dataset saved to {os.path.join(output_dir, f'pseudo_bulk_{cell_type_column}.tsv')}")
    logging.info(f"Amount of columns/celltypes {len(pseudo_bulk_df.columns)}")
    
    
    return os.path.join(output_dir, f"pseudo_bulk_{cell_type_column}.tsv")


def create_pseudo_bulk_unfiltered(reference_adata:str, output_dir:str, cell_type_column:str="celltype_major"):
    """
    Create a pseudo-bulk dataset from reference data.
    This creates a tsv file with genes as rows and cell types as columns.
    Args:
        reference_adata: path to the reference adata object
        output_dir: path to the output directory
    """
    os.makedirs(output_dir, exist_ok=True)

    # load the reference adata object
    reference_adata = ad.read_h5ad(reference_adata)
    logging.info(f"Reference adata object loaded")
    logging.info(reference_adata)
    
    # unnormalize the counts: reverse log1p
    #reference_adata.X = np.expm1(reference_adata.X)
    
    # get the genes that are present in both adata objects
    
    
    cell_types = reference_adata.obs[cell_type_column].unique().tolist()
    logging.info(f"Amount of cell types: {len(cell_types)}")
    # create a dataframe with the first column called gene and the rest of the columns called the cell types
    pseudo_bulk_df = pd.DataFrame(0,index=reference_adata.var_names, columns=cell_types)
    pseudo_bulk_df.index.name = "gene"
    
    # the reference adata object is already normalized to 10_000
    # average the counts across the cell types
    for cell_type in reference_adata.obs[cell_type_column].unique():
        logging.info(f"Processing cell type: {cell_type}")
        cell_type_adata = reference_adata[reference_adata.obs[cell_type_column] == cell_type, :].copy()
        counts_mean = cell_type_adata.X.toarray().mean(axis=0)
        logging.info(f"Counts mean: {counts_mean.shape}")
        # add the counts to the pseudo_bulk_df
        pseudo_bulk_df[cell_type] = counts_mean
        
    
    # ficture requires the names to be of format name_1, name_2, name_3, etc. no spaces allowed and no - 
    pseudo_bulk_df.columns = [col.replace(" ", "").replace("-", "") for col in pseudo_bulk_df.columns]
    # remove any numbers and + signs
    pseudo_bulk_df.columns = [re.sub(r'\d+', '', col) for col in pseudo_bulk_df.columns]
    pseudo_bulk_df.columns = [re.sub(r'\+', '', col) for col in pseudo_bulk_df.columns]
    pseudo_bulk_df.columns = [f"{col}_{i}" for i,col in enumerate(pseudo_bulk_df.columns)]
    
    # save the pseudo_bulk_df
    pseudo_bulk_df.to_csv(os.path.join(output_dir, f"pseudo_bulk_{cell_type_column}.tsv.gz"), sep="\t", compression="gzip")
    logging.info(f"Pseudo-bulk dataset saved to {os.path.join(output_dir, f'pseudo_bulk_{cell_type_column}.tsv.gz')}")
    logging.info(f"Amount of columns/celltypes {len(pseudo_bulk_df.columns)}")
    logging.info(f"Columns: {pseudo_bulk_df.columns}")
    # print the order of the columns
    logging.info(f"Order of columns: {pseudo_bulk_df.columns.tolist()}")
    
    return os.path.join(output_dir, f"pseudo_bulk_{cell_type_column}.tsv.gz"), pseudo_bulk_df.columns.tolist()



def run_ficture_prior(in_tsv:str,psuedo_bulk_tsv:str, out_dir:str, in_minmax:str,
                      bgzip_path:str, tabix_path:str, n_factors:str, anchor_res:float,
                      train_width:str, mu_scale:int=1, n_jobs:int=10):
    """
    Run ficture using prior knowledge.
    - plot each factor
    - decode sub um per pixel
    - decode plot um per pixel
    - minibatch buffer (can be skipped)
    """
    
    # parse input parameters
    train_widths = [int(x) for x in train_width.split(",")]
    n_factors = [int(x) for x in n_factors.split(",")]
    
    batch_tsv = f"{out_dir}/batched.matrix.tsv"
    batch_out = f"{out_dir}/batched.matrix.tsv.gz"
    minmax_out = in_minmax if in_minmax is not None else f"{out_dir}/coordinate_minmax.tsv"
    
    # create output directory
    os.makedirs(out_dir, exist_ok=True)
    
    min_ct_per_unit = 0
    min_ct_per_feature = 0
    
    mm=minimake()
    
    # preprocessing
    logging.info("Preprocessing...")
    logging.info(f"Creating minibatch from {in_tsv}")
    
    cmds = []
    cmds.append(rf"$(info --------------------------------------------------------------)")
    cmds.append(rf"$(info Creating minibatch from {in_tsv}...)")
    cmds.append(rf"$(info --------------------------------------------------------------)")
    ## create minibatch
    cmds.append(f"ficture make_spatial_minibatch --input {in_tsv} --output {batch_tsv} --mu_scale {mu_scale} --batch_size 500 --batch_buff 100 --major_axis Y")
    cmds.append(f"sort -k 2,2n -k 1,1g {batch_tsv} | {bgzip_path} -c > {batch_out}")
    cmds.append(f"rm {batch_tsv}")

    mm.add_target(batch_out, [in_tsv], cmds)
    
    ### segment
    for train_width in train_widths:
        dge_out = f"{out_dir}/hexagon.d_{train_width}.tsv"
        cmds = []
        cmds.append(rf"$(info --------------------------------------------------------------)")
        cmds.append(rf"$(info Creating DGE for {train_width}um...)")
        cmds.append(rf"$(info --------------------------------------------------------------)")
        cmds.append(f"ficture make_dge --key Count --input {in_tsv} --output {dge_out} --hex_width {train_width} --mu_scale {mu_scale} --min_ct_per_unit {min_ct_per_unit} --precision 2 --major-axis Y")
        cmds.append(f"sort -k 1,1n {dge_out} | {bgzip_path} -c > {dge_out}.gz")
        cmds.append(f"rm {dge_out}")
        mm.add_target(f"{dge_out}.gz", [in_tsv], cmds)
    
    ### lda
    for train_width in train_widths:
        for n_factor in n_factors:
            model_id=f"nF{n_factor}.d_{train_width}"
            model_path=f"{out_dir}/analysis/{model_id}"
            figure_path=f"{model_path}/figure"
            hexagon = f"{out_dir}/hexagon.d_{train_width}.tsv.gz"
            model_prefix=f"{model_path}/{model_id}"
            model=f"{model_prefix}.model.p"

            cmds = []
            cmds.append(rf"$(info --------------------------------------------------------------------------)")
            cmds.append(rf"$(info Creating LDA from pseudobulk for {train_width}um and {n_factor} factors...)")
            cmds.append(rf"$(info --------------------------------------------------------------------------)")
            cmds.append(f"mkdir -p {model_path}/figure")
            # --model argument needs to be used (check ficture github issue)
            ### import parameters for filtering that is done by ficture
            # min_ct_per_feature default is 50 -> too high for 2um resolution
            # min_ct_per_unit default is 50 -> too high for 2um resolution
            cmds.append(f"ficture init_model_from_pseudobulk --input {hexagon} --output {model_prefix} --model {pseudo_bulk_tsv} --key Count --min_ct_per_feature {min_ct_per_feature} --min_ct_per_unit {min_ct_per_unit}")
            cmds.append(f"cp {model_prefix}.posterior.count.tsv.gz {model_prefix}.model_matrix.tsv.gz")

            fit_tsv=f"{model_path}/{model_id}.fit_result.tsv.gz"
            fig_prefix=f"{figure_path}/{model_id}"
            cmap=f"{figure_path}/{model_id}.rgb.tsv"
            cmds.append(f"ficture choose_color --input {fit_tsv} --output {fig_prefix} --cmap_name turbo")

            fillr = (train_width / 2 + 1)
            cmds.append(f"ficture plot_base --input {fit_tsv} --output {fig_prefix}.coarse --fill_range {fillr} --color_table {cmap} --plot_discretized")
            cmds.append(f"touch {model_prefix}.done")

            mm.add_target(f"{model_prefix}.done", [in_tsv, hexagon], cmds)
    
    
    ### decode
    script_path = f"{out_dir}/sort_decode.sh"
    with open(script_path, "w") as f:
        f.write(r"""#!/bin/bash
input=$1
output=$2
coor=$3
model_id=$4
bsize=$5
scale=$6
topk=$7
gzip=$8
#bgzip=$8
#tabix=$9

K=$( echo $model_id | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' )
while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
bsize=2000
scale=100
header="##K=${K};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" | sort -S 1G -k1,1g -k3,3g ) | ${gzip} -c > ${output}

#${tabix} -f -s1 -b3 -e3 ${output}
rm ${input}
""")
        for train_width in train_widths:
            for n_factor in n_factors:
                batch_in = f"{out_dir}/batched.matrix.tsv.gz"
                model_id=f"nF{n_factor}.d_{train_width}"
                model_path=f"{out_dir}/analysis/{model_id}"
                figure_path=f"{model_path}/figure"
                model_prefix=f"{model_path}/{model_id}"
                cmap=f"{figure_path}/{model_id}.rgb.tsv"
                model=f"{out_dir}/analysis/{model_id}/{model_id}.model_matrix.tsv.gz"

                fit_widths = [int(train_width)]
                
                for fit_width in fit_widths:
                    cmds = []

                    fit_nmove = int(fit_width / int(anchor_res))
                    anchor_info=f"prj_{fit_width}.r_{anchor_res}"
                    radius = int(anchor_res) + 1

                    prj_prefix = f"{model_path}/{model_id}.{anchor_info}"
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(rf"$(info Creating projection for {train_width}um and {n_factor} factors, at {fit_width}um)")
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(f"ficture transform --input {in_tsv} --output_pref {prj_prefix} --model {model} --key Count --hex_width {fit_width} --n_move {fit_nmove} --mu_scale {mu_scale} --thread 5")

                    batch_input=f"{out_dir}/batched.matrix.tsv.gz"
                    anchor=f"{prj_prefix}.fit_result.tsv.gz"
                    decode_basename=f"{model_id}.decode.{anchor_info}_{radius}"
                    decode_prefix=f"{model_path}/{decode_basename}"

                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(rf"$(info Performing pixel-level decoding..)")
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(f"ficture slda_decode --input {batch_in} --output {decode_prefix} --model {model} --anchor {anchor} --anchor_in_um --neighbor_radius {radius} --mu_scale {mu_scale} --key Count --precision 0.1 --lite_topk_output_pixel 3 --lite_topk_output_anchor 3 --thread 3")

                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(rf"$(info Sorting and reformatting the pixel-level output..)")
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(f"bash {script_path} {decode_prefix}.pixel.tsv.gz {decode_prefix}.pixel.sorted.tsv.gz {minmax_out} {model_id} 100 100 3 {bgzip_path}")

                    de_input=f"{decode_prefix}.posterior.count.tsv.gz"
                    de_output=f"{decode_prefix}.bulk_chisq.tsv"

                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(rf"$(info Performing pseudo-bulk differential expression analysis..)")
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(f"ficture de_bulk --input {de_input} --output {de_output}")

                    cmap=f"{figure_path}/{model_id}.rgb.tsv"
                    cmds.append(f"ficture factor_report --path {model_path} --pref {decode_basename} --color_table {cmap}")

                    decode_tsv=f"{decode_prefix}.pixel.sorted.tsv.gz"
                    decode_png=f"{model_path}/figure/{decode_basename}.pixel.png"

                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(rf"$(info Drawing pixel-level output image...)")
                    cmds.append(rf"$(info --------------------------------------------------------------)")
                    cmds.append(f"ficture plot_pixel_full --input {decode_tsv} --color_table {cmap} --output {decode_png}  --full")

                    # plot each factor
                    if True:
                        sub_prefix=f"{model_path}/figure/sub/{decode_basename}.pixel"
                        cmds.append(f"mkdir -p {model_path}/figure/sub")
                        cmds.append(f"ficture plot_pixel_single --input {decode_tsv} --output {sub_prefix}  --full --all")

                    cmds.append(f"touch {decode_prefix}.done")
                    mm.add_target(f"{decode_prefix}.done", [batch_in, hexagon,f"{model_prefix}.done"], cmds)

    mm.write_makefile(f"{out_dir}/Makefile")
    os.system(f"make -f {out_dir}/Makefile -j {n_jobs}")


def create_cmap(factor_names_order:list, output_dir:str):
    """
    Create a cmap of the cell types
    
    Args:
        factor_names_order: list of factor names in order
        output_dir: path to the output directory
        
    Returns:
        path to the created color map file
    """
    
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Number of factors
    n_factors = len(factor_names_order)
    logging.info(f"Creating color map for {n_factors} factors")
    
    # Create a dataframe for the color map
    cmap_df = pd.DataFrame({
        'Name': range(n_factors),
        'Annotation': factor_names_order,
        'Color_index': range(n_factors),  # Can be adjusted if specific color indices are preferred
        'R': np.zeros(n_factors),
        'G': np.zeros(n_factors),
        'B': np.zeros(n_factors)
    })
    
    # Generate colors using a colormap (turbo is often used in ficture)
    # If matplotlib is available, use it for better colors
    colormap = cm.get_cmap('turbo', n_factors)
    
    for i in range(n_factors):
        color = colormap(i / max(1, n_factors - 1))
        cmap_df.loc[i, 'R'] = color[0]
        cmap_df.loc[i, 'G'] = color[1]
        cmap_df.loc[i, 'B'] = color[2]
        
       
    # Output path
    cmap_path = os.path.join(output_dir, "cmap.rgb.tsv")
    
    # Save the color map
    cmap_df.to_csv(cmap_path, sep='\t', index=False)
    logging.info(f"Color map saved to {cmap_path}")
    
    return cmap_path


def preprocess_ficture_prior(reference_adata_path:str, output_dir:str, cell_type_column:str="celltype_major"):
    """
    Preprocess the reference adata object for ficture prior.
    - create a pseudo-bulk dataset
    - create a cmap
    """
    
    # create a pseudo-bulk dataset
    pseudo_bulk_tsv_path, factor_names_order = create_pseudo_bulk_unfiltered(reference_adata_path, output_dir, cell_type_column)
    logging.info(f"Created pseudo-bulk dataset: {pseudo_bulk_tsv_path}")
    logging.info(f"Factor names order: {factor_names_order}")
    # create a cmap of the cell types
    cmap_path = create_cmap(factor_names_order, output_dir)
    logging.info(f"Created color map: {cmap_path}")
    
    logging.info(f"Preprocessing complete. Pseudo-bulk data: {pseudo_bulk_tsv_path}, Color map: {cmap_path}")
    
    return pseudo_bulk_tsv_path, cmap_path
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference_adata_path", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default=None)
    parser.add_argument("--cell_type_column", type=str, default="celltype_major")
    args = parser.parse_args()
    
    if args.reference_adata_path is None:
        args.reference_adata_path = "resources_data/BrCa_Atlas_Count_out/annotated_brca_atlas.h5ad"
        args.output_dir = "results/human_breast_cancer_final"
        args.cell_type_column = "celltype_major"
    
    # Preprocess the reference data to create pseudo-bulk and color map
    pseudo_bulk_tsv_path, cmap_path = preprocess_ficture_prior(args.reference_adata_path, args.output_dir, args.cell_type_column)
    


    
    
    



