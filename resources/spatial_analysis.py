import argparse
import anndata as ad
import squidpy as sq
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from loguru import logger as logging
import scanpy as sc
import scipy.sparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CPUS = int(os.getenv("SLURM_CPUS_PER_TASK", 3))
logging.info(f"using {CPUS} CPUs")

DPI = 300
FIGSIZE = (10, 10)

def main(input:ad.AnnData, output_dir:str, cluster_key:str):
    
    logging.info(f"input: {input}")
    logging.info(f"output_dir: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"cluster_key: {cluster_key}")
    
    # load the input data
    adata = ad.read_h5ad(input)
    adata.var_names_make_unique()
    logging.info(f"adata: {adata}")
        
    # spatial enrichment analysis -> based on neighborhood connectivity
    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, dpi=DPI, figsize=FIGSIZE)
    plt.savefig(os.path.join(output_dir, "nhood_enrichment.png"), dpi=DPI)
    logging.info(f"nhood_enrichment done")
    
    # spatial enrichment analysis based on coordinates with increasing radius
    sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    for cluster in adata.obs[cluster_key].unique():
        sq.pl.co_occurrence(
            adata,
            cluster_key=cluster_key,
            clusters=cluster,
            figsize=FIGSIZE,
            dpi=DPI
        )
        plt.savefig(os.path.join(output_dir, f"co_occurrence_{cluster}.png"), dpi=DPI)
    logging.info(f"co_occurrence done")
    
    # ligand-receptor analysis -> implementation of cellphonedb 
    
    # spatially variable genes with moran's i -> results saved in adata.uns['moranI']
    #genes = adata[:, adata.var.highly_variable].var_names.values[:1000]
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
    #    genes=genes,
        n_perms=100,
        n_jobs=CPUS,
        figsize=FIGSIZE,
        dpi=DPI
    )
    print("Moran's I:")
    print(adata.uns["moranI"].head(10))
    moran_df = adata.uns["moranI"]
    moran_df.to_csv(os.path.join(output_dir, "moranI.csv"), index=False)
    # select the indices of the top 5 genes
    top_genes = adata.uns["moranI"].sort_values(by="I", ascending=True).head(5).index.tolist()
    top_scores = adata.uns["moranI"].sort_values(by="I", ascending=True).head(5)["I"].tolist()
    print("Moran's I top 5 genes:")
    print(top_genes)
    # plot the top 5 genes and the clusters
    cols = top_genes + [cluster_key]
    sq.pl.spatial_scatter(
        adata,
        color=cols,
        figsize=FIGSIZE,
        dpi=DPI,
        title=f"Top 5 Moran's I genes:{top_scores}",
        cmap="hot"
    )
    plt.savefig(os.path.join(output_dir, "spatial_scatter_moran.png"), dpi=DPI)
    
    
    lowest_genes = adata.uns["moranI"].sort_values(by="I", ascending=False).head(5).index.tolist()
    lowest_scores = adata.uns["moranI"].sort_values(by="I", ascending=False).head(5)["I"].tolist()
    print("Moran's I lowest 5 genes:")
    print(lowest_genes)
    # plot the lowest 5 genes and the clusters
    cols = lowest_genes + [cluster_key]
    sq.pl.spatial_scatter(
        adata,
        color=cols,
        figsize=FIGSIZE,
        dpi=DPI,
        title=f"Lowest 5 Moran's I genes:{lowest_scores}",
        cmap="hot"
    )
    plt.savefig(os.path.join(output_dir, "spatial_scatter_moran_lowest.png"), dpi=DPI)
    
    return None
    

def from_source_spatial(binned_output_dir:str, output_dir:str):
    os.makedirs(output_dir, exist_ok=True)
    # load the binned output
    binned_adata = sc.read_visium(binned_output_dir)
    binned_adata.var_names_make_unique()
    logging.info(f"binned data: {binned_output_dir}")
    logging.info(f"binned_adata:")
    logging.info(binned_adata)
    
    
    
    # check if Foxa1 and Il33 are in the genes
    if "FOXA1" not in binned_adata.var_names or "IL33" not in binned_adata.var_names:
        # find the closest gene to foxa1 and il33
        logging.info("Foxa1 and Il33 are not in the genes")
        # find similar genes using string matching
        foxa1_similar = [gene for gene in binned_adata.var_names if "FOXA1" in gene.upper() or "Foxa1" in gene]
        il33_similar = [gene for gene in binned_adata.var_names if "IL33" in gene.upper() or "Il33" in gene]
        
        logging.info(f"Similar genes to Foxa1: {foxa1_similar if foxa1_similar else 'None found'}")
        logging.info(f"Similar genes to Il33: {il33_similar if il33_similar else 'None found'}")
    
        
    binned_adata = binned_adata[:, ["FOXA1", "IL33"]]
    print("binned_adata: ", binned_adata)
    
    # Get the expression values and convert to dense arrays
    foxa1_expr = binned_adata.X[:, 0].toarray().flatten()
    il33_expr = binned_adata.X[:, 1].toarray().flatten()
    
    # Create mask using the dense arrays
    mask = (foxa1_expr > 0) & (il33_expr > 0)
    binned_adata = binned_adata[mask, :]
    
    # Get indices for FOXA1 and IL33 in var_names
    # foxa1_idx = np.where(binned_adata.var_names == "FOXA1")[0][0]
    # il33_idx = np.where(binned_adata.var_names == "IL33")[0][0]
    # print("foxa1_idx: ", foxa1_idx)
    # print("il33_idx: ", il33_idx)
    # # Get expression values from the count matrix
    # foxa1_expr = binned_adata.X[:, foxa1_idx]
    # il33_expr = binned_adata.X[:, il33_idx]
    
    # # If the matrix is sparse, convert to dense array
    # if scipy.sparse.issparse(binned_adata.X):
    #     foxa1_expr = foxa1_expr.toarray().flatten()
    #     il33_expr = il33_expr.toarray().flatten()
    
    # # Keep only rows where both genes have non-zero expression
    # mask = (foxa1_expr > 0) & (il33_expr > 0)
    # binned_adata = binned_adata[mask, :]
    
    print("binned_adata: ", binned_adata)
    
    # loop over the binned_adata and print the counts
    for i in range(binned_adata.shape[0]):
        print(f"cell {i}, foxa1: {binned_adata.X[i, 0]}, il33: {binned_adata.X[i, 1]}")
    
    # plot the foxa1 and il33 counts
    sq.pl.spatial_scatter(
        binned_adata,
        color=["FOXA1", "IL33"],
        figsize=FIGSIZE,
        dpi=DPI,
        img_alpha=0.7,
        size=10
    )
    plt.savefig(os.path.join(output_dir, "anndata_subset.png"), dpi=DPI)
    plt.close()
    logging.info(f"spatial scatter done")
    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=False, default=None)
    parser.add_argument("--cluster_key", type=str, required=False, default="factor")
    parser.add_argument("--output_dir", type=str, required=False, default=None)
    args = parser.parse_args()
    if args.input is None:
        args.input = "results/human_breast_cancer_final/ficture/annotated_anndata.h5ad"
        args.cluster_key = "factor"
        args.output_dir = "results/human_breast_cancer_final/spatial_analysis_008"
    
    main(args.input, args.output_dir, args.cluster_key)
    #from_source_spatial("input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/binned_outputs/square_008um", args.output_dir)