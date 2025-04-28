"""
Annotate the cells using the celltypist model.
If the model is not provided, train the model on the reference atlas and column.
"""

import sys
import os
import argparse

import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns
import celltypist
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score, ConfusionMatrixDisplay
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from loguru import logger as logging

# add the project root to the python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


### BRCA ATLAS - 49 celltypes adata object
BRCA_ATLAS_ADATA = "resources_data/BrCa_Atlas_Count_out/annotated_brca_atlas.h5ad"
CPUS = os.getenv("SLURM_CPUS_PER_TASK", 5)
CPUS = int(CPUS)
logging.info(f"Using {CPUS} CPUs")

# set scanpy plot settings
sc.set_figure_params(figsize=(15, 15), dpi=100)


def train_celltypist_model(atlas_adata_path:str,cell_label_name:str="celltype_subset",output_dir:str="resources_data/celltypist"):
    """
    Train a celltypist model on an ATLAS adata object.
    The atlas adata object should have the following columns in the obs:
     - celltype_major -> 9 major celltypes (B-cells, T-cells, Endothelial, ...)
     - celltype_subset -> 49 subcelltypes
    
    The model returns the test set predictions in anndata object. 
     
    Args:
        atlas_adata: adata object with the atlas data
        cell_label_name: name of the column in the obs that contains the cell labels to train the model on
        output_dir: directory to save the model
        replace: if True, replace the existing model, redo the training
    """
    
    model_output_path = os.path.join(output_dir, f"model_{cell_label_name}.pkl")
    test_output_path = os.path.join(output_dir, f"test_preds_{cell_label_name}.h5ad")
    
    atlas_adata = ad.read_h5ad(atlas_adata_path)
    
    # check that the anndata has been normalized to 10_000 and log1p transformed for the celltypist training
    if not np.allclose(atlas_adata.X.toarray().sum(axis=1), 10_000):
        logging.info("The atlas adata object has not been normalized to 10_000")
        sc.pp.normalize_total(atlas_adata, inplace=True, target_sum=10_000)
        sc.pp.log1p(atlas_adata)
    else:
        logging.info("The atlas adata object has been normalized to 10_000 and log1p transformed")

    # split the data into training and test set - stratified by celltype_subset
    train_data_df = atlas_adata.obs.copy()
    train_df, test_df = train_test_split(train_data_df, test_size=0.55, stratify=train_data_df[cell_label_name])
    train_adata, test_adata = atlas_adata[train_df.index, :], atlas_adata[test_df.index, :]
    
    new_model = celltypist.train(train_adata, labels=cell_label_name, n_jobs=CPUS)
    new_model.write(model_output_path)
    
    # predict the cell types from the test set
    predictions = celltypist.annotate(test_adata, model=model_output_path, majority_voting=True)
    results = predictions.to_adata()
    
    # save the adata object with the test set predictions
    results.write_h5ad(test_output_path)
    
    return model_output_path

def evaluate_celltypist_model(results_adata:ad.AnnData, output_dir:str="resources_data/celltypist/model_eval", cell_label_name:str="celltype_subset"):
    """
    Evaluate the celltypist model on the test set.
    
    Args:
        results_adata: adata object with the test set predictions, has the original cell labels and the predicted cell labels
        output_dir: directory to save the evaluation plots
        cell_label_name: name of the column in the obs that contains the cell labels to train the model on
    """
    
    true_labels = results_adata.obs[cell_label_name]
    predicted_labels_majority = results_adata.obs['majority_voting']
    predicted_labels = results_adata.obs['predicted_labels']
    
    ### Create the umap plot with celltype_subset, predicted_labels and majority_voting
    sc.tl.umap(results_adata)
    sc.pl.umap(results_adata, color=['celltype_subset', 'predicted_labels', 'majority_voting'])
    plt.savefig(os.path.join(output_dir, "celltypist_umap.png"))
    plt.close()
    
    # print some statistics about the predicted labels and accuracy
    logging.info(f"Accuracy predicted_labels: {accuracy_score(true_labels, predicted_labels)}")
    logging.info(f"Accuracy majority_voting: {accuracy_score(true_labels, predicted_labels_majority)}")
    
    
    # generate a confusion matrix
    labels = true_labels.unique()
    cm = confusion_matrix(true_labels, predicted_labels, labels=labels)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    disp.plot(include_values=True, cmap='Blues', ax=None, xticks_rotation='horizontal')
    plt.show()
    plt.savefig(os.path.join(output_dir, "celltypist_confusion_matrix.png"))
    plt.close()
    

def celltypist_annotate(adata_to_annotate:str, model_path:str, reference_atlas_path:str, cell_type_column:str, output_dir:str):
    """
    Annotate the input adata object with the celltypist model.
    If the model is not provided or present, train the model on the reference atlas and cell type column.
    """
    
    # check if the model is provided and exists
    if model_path is None or not os.path.exists(model_path):
        logging.info(f"Model not provided or does not exist, training the model on the reference atlas")
        # train the model
        model_path = train_celltypist_model(reference_atlas_path, cell_type_column)

    # annotate the input adata object
    adata_to_annotate = ad.read_h5ad(adata_to_annotate) 
    # check if the adata object is normalized to 10_000 and log1p transformed
    if not np.allclose(np.expm1(adata_to_annotate.X.toarray()).sum(axis=1), 10_000):
        logging.info("The adata object has not been normalized to 10_000")
        sc.pp.normalize_total(adata_to_annotate, inplace=True, target_sum=10_000)
        sc.pp.log1p(adata_to_annotate)
    else:
        logging.info("The adata object has been normalized to 10_000 and log1p transformed")
        
    logging.info(f"Annotating the input adata object with the celltypist model")
    predictions = celltypist.annotate(adata_to_annotate, model=model_path, majority_voting=True)
    results = predictions.to_adata()
    logging.info(f"Saving the results to {output_dir}")
    # save the results
    results.write_h5ad(os.path.join(output_dir, "celltypist_annotated.h5ad"))
    logging.info(f"Results saved to {output_dir}")
    
    return None



def celltypist_annotate_visium(input_path:str, output_path:str, model_path:str="resources_data/celltypist/model_celltype_subset.pkl", ficture_color_table:str=None):
    print(input_path)
    os.makedirs(output_path, exist_ok=True)
    # check if the input path ends with .h5ad
    if not input_path.endswith(".h5ad"):
        visium = sc.read_visium(input_path)
    else:
        visium = ad.read_h5ad(input_path)
        
    if not os.path.exists(os.path.join(output_path, "celltypist_annotated_visium.h5ad")):
        print("genes in visium")
        print(visium.var_names)
        print("Read the visium data")
        print(visium)
        # basic filtering
        sc.pp.filter_genes(visium, min_cells=10)
        sc.pp.filter_cells(visium, min_counts=10)
        print("Filtered the visium data")
        print(visium)
    
        # if the data is not normalized to 10_000, normalize it
        if not np.allclose(np.expm1(visium.X.toarray()).sum(axis=1), 10_000):
            sc.pp.normalize_total(visium, inplace=True, target_sum=10_000)
            sc.pp.log1p(visium)
            
        print("Normalized the visium data")
        print(visium)
        
        
        # annotate the visium object - skipping the majority voting as it requires too much memory for 2um data
        predictions = celltypist.annotate(visium, model=model_path, majority_voting=True)
        decision_matrix = predictions.decision_matrix
        decision_matrix.to_csv(os.path.join(output_path, "celltypist_decision_matrix.csv"))
        prob_matrix = predictions.probability_matrix
        prob_matrix.to_csv(os.path.join(output_path, "celltypist_probability_matrix.csv"))
        
        results = predictions.to_adata()
    else:
        results = ad.read_h5ad(os.path.join(output_path, "celltypist_annotated_visium.h5ad"))
    
    
    
    print(results)
    # print the valuecounts of the predicted labels
    print(results.obs['predicted_labels'].value_counts())
    
    # turn the predicted labels into a categorical variable
    results.obs['predicted_labels'] = pd.Categorical(results.obs['predicted_labels'])
    
    # if the ficture_color_table is provided, use it to color the predicted labels
    if ficture_color_table is not None:
        # read the ficture_color_table
        ficture_color_table = pd.read_csv(ficture_color_table, sep="\t")
        # create a color mapping from the ficture_color_table which has R G B columns in decimal format
        color_mapping = dict(zip(ficture_color_table['Annotation'], ficture_color_table[['R', 'G', 'B']].apply(lambda x: tuple(x), axis=1)))
        
    
    # spatial scatter plot of the predicted labels, disable tight layout
    plt.figure(figsize=(20, 15))  # Increase figure size
    sc.pl.spatial(results, color=['predicted_labels'], show=False, alpha_img=0.2, size=2, palette=color_mapping)
    plt.tight_layout()  # Add tight_layout to adjust spacing
    plt.savefig(os.path.join(output_path, "celltypist_spatial_scatter.png"), bbox_inches='tight', dpi=300)  # Add bbox_inches='tight' to prevent legend cutoff
    plt.close()
    
    plt.figure(figsize=(20, 15))  # Increase figure size
    sc.pl.spatial(results, color=['majority_voting'], show=False, alpha_img=0.2, size=2, palette=color_mapping)
    plt.tight_layout()  # Add tight_layout to adjust spacing
    plt.savefig(os.path.join(output_path, "celltypist_spatial_scatter_majority_voting.png"), bbox_inches='tight', dpi=300)
    plt.close()
    
    # plot the umap of the results
    sc.tl.umap(results)
    sc.pl.umap(results, color=['predicted_labels', 'majority_voting'], size=2,palette=color_mapping)
    plt.savefig(os.path.join(output_path, "celltypist_umap.png"))
    plt.close()
    
    
    # save the results to file
    if not os.path.exists(os.path.join(output_path, "celltypist_annotated_visium.h5ad")):
        results.write_h5ad(os.path.join(output_path, "celltypist_annotated_visium.h5ad"))
    else:
        return None
    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_anndata_path", type=str, default=None)
    parser.add_argument("--model_path", type=str, default=None)
    parser.add_argument("--reference_atlas_path", type=str, default=None)
    parser.add_argument("--cell_type_column", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default=None)
    args = parser.parse_args()
    
    
    celltypist_annotate(
        adata_to_annotate=args.input_anndata_path,
        model_path=args.model_path,
        reference_atlas_path=args.reference_atlas_path,
        cell_type_column=args.cell_type_column,
        output_dir=args.output_dir
    )
    
    
    
    
    