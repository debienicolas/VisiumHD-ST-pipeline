""" 
Comparison of annotations.

"""
import os
import argparse
from collections import defaultdict
import sys

import anndata as ad
import numpy as np
from scipy.spatial import cKDTree
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import scipy
import json

from loguru import logger as logging

# append the parent directory to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from annotation.post_annotation import find_marker_genes


def compare_cells(anndata_1:str, anndata_2:str, annotation_column:str, output_path:str):
    """
    Compare the cells of two anndata objects based on a given annotation column.
    Keep track of:
        - differences in total cells
        - differences in cell counts per annotation category
    
    Plot the differences spatially.
    """
    os.makedirs(output_path, exist_ok=True)
    
    # load the anndata objects
    adata_1 = ad.read_h5ad(anndata_1)
    adata_2 = ad.read_h5ad(anndata_2)
    
    if annotation_column not in adata_1.obs.columns or annotation_column not in adata_2.obs.columns:
        raise ValueError(f"Annotation column {annotation_column} not found in anndata objects")
    
    logging.info(f"Comparing {len(adata_1.obs)} cells from {anndata_1} with {len(adata_2.obs)} cells from {anndata_2}")
    
    # build a tree based on the spatial positions of the largest anndata object
    if len(adata_1.obs) > len(adata_2.obs):
        tree = cKDTree(adata_1.obsm["spatial"])
    else:
        logging.info("Swapping the anndata objects")
        adata_1, adata_2 = adata_2, adata_1
        tree = cKDTree(adata_1.obsm["spatial"])
    
    distances, indices = tree.query(adata_2.obsm["spatial"], k=1)
    
    buffer = 10     
    perfect_matches = len(np.where(distances == 0)[0])
    logging.info(f"Perfect spatial matches: {perfect_matches}")
    within_buffer = len(np.where(distances < buffer)[0]) - perfect_matches
    logging.info(f"Within buffer: {within_buffer}")
    
    # turn the annotation column into an int and then string
    logging.info(f"Factor types: {adata_1.obs[annotation_column].unique()}, {adata_1.obs[annotation_column].dtype}")
    logging.info(f"Factor types: {adata_2.obs[annotation_column].unique()}, {adata_2.obs[annotation_column].dtype}")
    adata_1.obs[annotation_column] = adata_1.obs[annotation_column].astype(int).astype(str)
    adata_2.obs[annotation_column] = adata_2.obs[annotation_column].astype(int).astype(str)
    
    annotations_1 = adata_1.obs[annotation_column].unique().tolist()
    annotations_2 = adata_2.obs[annotation_column].unique().tolist()
    
    # sort the annotations
    annotations_1.sort()
    annotations_2.sort()
    print(annotations_1)
    print(annotations_2)
    
    
    
    # create a decision matrix
    decision_matrix = defaultdict(int)
    # print the factor types
    adata_2.obs["diff_annot"] = "True"
    for i, cell in enumerate(adata_2.obs.index):
        index_1 = indices[i]
        decision_matrix[str(adata_1.obs[annotation_column].iloc[index_1]) + "_" + str(adata_2.obs[annotation_column].iloc[i])] += 1
        if adata_1.obs[annotation_column].iloc[index_1] == adata_2.obs[annotation_column].iloc[i]:
            adata_2.obs.loc[cell, "diff_annot"] = "False"
            
    print(decision_matrix)
    
    # plot the differences spatially -  only plot the cells that have a difference 
    sc.pl.spatial(adata_2, color="diff_annot", title="Difference in annotation", color_map="bwr")
    plt.savefig(os.path.join(output_path, "diff_annot.png"))

def within_cluster_variance(X:np.ndarray, labels:pd.Series, output_dir:str):
    """
    Calculate the within-cluster variance of a dataset
    
    This might need dimensionality reduction to calculate a representative within-cluster variance and silhouette score (curse of dimensionality).
    
    Args:
        X: numpy array of shape (n_samples, n_features)
        labels: pandas series of shape (n_samples,)
        output_dir: path to the output directory
    """
    # create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    results = {}
    expression_profiles = {}
    
    # calculate the within-cluster variance
    variances = []
    for label in labels.unique():
        cluster_data = X[labels == label]
        mean = np.mean(cluster_data, axis=0)
        variance = np.mean(np.sum((cluster_data - mean) ** 2, axis=1))
        variances.append(variance)
        results[f"{label}_variance"] = variance
        
        expression_profiles[label] = mean
    
    results["total_variance"] = sum(variances)
    
    # calculate the silhouette score 
    sil_score = silhouette_score(X, labels)
    results["silhouette_score"] = sil_score
    
    # save the results
    with open(os.path.join(output_dir, "within_cluster_variance.json"), "w") as f:
        json.dump(results, f)
    
    # plot the expression profiles in a smoothed lineplot
    
    for label, profile in expression_profiles.items():
        plt.plot(profile, label=label)    
    plt.legend()
    plt.savefig(os.path.join(output_dir, "expression_profiles.png"))
    
    return



def compare_annotations(supervised_paths:dict, unsupervised_paths:dict, label_column:str,output_dir:str):
    """
    Compare the annotations/clustering of supervised and unsupervised methods
    
    1. Determine gene markers for each cluster
    2. Comparison of expression profiles through cluster centroids
    
    Args:
        supervised_paths: dictionary with the method name as key and the path to the anndata object as value
        unsupervised_paths: dictionary with the method name as key and the path to the anndata object as value
        output_dir: path to the output directory
    """
    # create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
        
    # Determine the gene markers for each method and cluster
    for method, path in supervised_paths.items():
        
        logging.info(f"Finding gene markers for {method} in {path}")
        
        # create the output and figures directories
        method_output_dir, figures_dir = os.path.join(output_dir, method), os.path.join(output_dir, method, "figures")
        os.makedirs(method_output_dir, exist_ok=True)
        os.makedirs(figures_dir, exist_ok=True)
        
        # load the anndata object
        adata = ad.read_h5ad(path)
        logging.info(f"Loaded anndata object with {adata.n_obs} cells and {adata.n_vars} genes")
        
        # check that the label column is present in the anndata object
        if label_column not in adata.obs.columns:
            raise ValueError(f"Label column {label_column} not found in {method} anndata object")
        
        # find the marker genes
        find_marker_genes(adata, label_column, figures_dir, method_output_dir)
        
        # calculate the within-cluster variance and silhouette score
        counts = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
        within_cluster_variance(counts, adata.obs[label_column], method_output_dir)
        
        
    pass



if __name__ == "__main__":
    # anndata_1 is the reference
    #{"celltypist": "celltypist_results/celltypist_annotated_visium.h5ad",
    compare_annotations(supervised_paths={"ficture_sup": "final_prior_test/2_um_post_b2c/postproc/postproc_ficture/ficture_anndata.h5ad"}, 
                       unsupervised_paths={}, 
                       label_column="factor",
                       output_dir="results/human_breast_cancer_final/annotation_comparison")