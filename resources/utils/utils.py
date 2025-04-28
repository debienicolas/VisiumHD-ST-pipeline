import squidpy as sq
import os
import subprocess
import anndata as ad
import pickle
# from gprofiler import GProfiler
import anndata as ad

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


def save_visium_data(adata: ad.AnnData, output_path: str):
    """
    Save Visium data to anndata object
    """

    
    # use pickle to save the anndata object
    with open(output_path, 'wb') as f:
        pickle.dump(adata, f)



def load_visium_data(input_path: str):
    """
    Load Visium data from pickle file
    """

    with open(input_path, 'rb') as f:
        adata = pickle.load(f)
        
    return adata



# def mouse_to_human_gene_names(gene_names: list[str]) -> dict[str, str]:
#     """
#     Convert mouse gene names to human gene names
#     """

#     gp = GProfiler(return_dataframe=True)
#     ortholog_mapping = gp.orth(organism='mmusculus', query=gene_names, target='hsapiens')
#     mouse_to_human_dict = dict(zip(ortholog_mapping['incoming'], ortholog_mapping['name']))
#     mouse_to_human_dict = {k: v if v != 'N/A' else k for k, v in mouse_to_human_dict.items()}
    
#     return mouse_to_human_dict