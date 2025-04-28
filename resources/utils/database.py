import sys
import os

import anndata as ad
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score
import numpy as np
import scipy.sparse as sp


from loguru import logger as logging

# add the project root to the python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def create_annotated_adata(adata:str, metadata_csv:str, output_path:str="resources_data/annotated_brca_atlas.h5ad"):
    
    logging.info(f"Creating AnnData object from {adata}")
    
    adata = sc.read_10x_mtx(adata)
    
    # add the metadata to the adata object
    metadata_df = pd.read_csv(metadata_csv, low_memory=False)
    
    # Remove the first row which appears to be a secondary header
    metadata_df = metadata_df.iloc[1:].reset_index(drop=True)
    
    # Create a dictionary mapping cell barcodes to metadata rows
    # Assuming the NAME column contains the cell barcodes that match adata.obs_names
    metadata_df.set_index('NAME', inplace=True)
    
    # Check if all cell barcodes from anndata exist in metadata
    missing_cells = set(adata.obs_names) - set(metadata_df.index)
    if missing_cells:
        print(f"Warning: {len(missing_cells)} cells in AnnData are not in metadata")
    
    # Add metadata to AnnData object for cells that exist in metadata
    for col in metadata_df.columns:
        adata.obs[col] = pd.Series([metadata_df.loc[cell, col] if cell in metadata_df.index else None 
                                 for cell in adata.obs_names], index=adata.obs_names)
    
    print(adata)
    
    # Save the annotated AnnData object if needed
    adata.write(output_path)


def check_if_normalized(adata:ad.AnnData):
    # get a sample of the matrix (first 1000 cells)
    if sp.issparse(adata.X):
        sample_matrix = adata.X[:1000, :].toarray() if adata.X.shape[0] > 1000 else adata.X.toarray()
    else:
        sample_matrix = adata.X[:1000, :]
    
    # check if the matrix is normalized
    print("Min value: ", np.min(sample_matrix))
    print("Max value: ", np.max(sample_matrix))
    print("Mean value: ", np.mean(sample_matrix))
    print("Median value: ", np.median(sample_matrix))
    print("Standard deviation: ", np.std(sample_matrix))
    
    
    # Check if the data sums to a specific value per cell
    if sp.issparse(adata.X):
        cell_sums = adata.X.sum(axis=1).A1
    else:
        cell_sums = adata.X.sum(axis=1)
    
    print(f"Mean sum per cell: {np.mean(cell_sums)}")
    print(f"Standard deviation of sums: {np.std(cell_sums)}")
    print(f"Are cell sums similar? {np.std(cell_sums) / np.mean(cell_sums) < 0.1}")
    
    

def train_celltypist_model(input_adata:ad.AnnData, output_path:str="resources_data/celltypist/brca_atlas_model.pkl", replace:bool=False):
    test_output_path = output_path.replace(".pkl", "test__predictions.h5ad")
    
    # check if the model already exists
    if os.path.exists(output_path) and os.path.exists(test_output_path) and not replace: 
        logging.info(f"Model already exists at {output_path}")
        test_adata = sc.read_h5ad(test_output_path)
        return test_adata
    
    logging.info(f"Training model at {output_path}")
    # train the model
    # the input counts matrix should be normalized to 10_000 and log1p transformed
    sc.pp.normalize_total(input_adata, inplace=True, target_sum=10000)
    sc.pp.log1p(input_adata)
    
    # build a test and train set, with the test set being 20% of the data stratified by celltype_subset
    training_data_df = input_adata.obs.copy()
    train_data_df, test_data_df = train_test_split(training_data_df, test_size=0.2, stratify=training_data_df['celltype_subset'])
    
    # convert the dataframes to annodata
    train_adata, test_adata = input_adata[train_data_df.index, :], input_adata[test_data_df.index, :]
    
    new_model = celltypist.train(train_adata, labels="celltype_subset", n_jobs=2, feature_selection=True, use_SGD=True)
    # save the model
    new_model.write(output_path)
    
    # predict the cell types: majority voting True -> informed using subcluster
    predictions = celltypist.annotate(test_adata, model=output_path, majority_voting=True)
    results = predictions.to_adata()
    
    results.write_h5ad(test_output_path)
    return results


if __name__ == "__main__":
    if not os.path.exists("resources_data/annotated_brca_atlas.h5ad"):
        create_annotated_adata(adata="resources_data/BrCa_Atlas_Count_out", metadata_csv="resources_data/Whole_miniatlas_meta.csv", output_path="resources_data/annotated_brca_atlas.h5ad")
    else:
        logging.info("Annotated AnnData object already exists")
        
    # print the unique values celltype_major obs column
    adata = sc.read_h5ad("resources_data/annotated_brca_atlas.h5ad")
    
    # 49 celltypes found in celltype_subset column
    
    results = train_celltypist_model(adata, "resources_data/celltypist/brca_atlas_model.pkl")
    
    
    sc.tl.umap(results)
    sc.pl.umap(results, color=['celltype_subset', 'predicted_labels', 'majority_voting'], save="_celltypist_umap.png")
    
    ### generate a confusion matrix
    from sklearn.metrics import confusion_matrix
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # Get the true and predicted labels
    true_labels = results.obs['celltype_subset']
    predicted_labels_majority = results.obs['majority_voting']
    predicted_labels = results.obs['predicted_labels']
    
    logging.info(f"True labels unique: {len(true_labels.unique())}")
    logging.info(f"Predicted labels unique: {len(predicted_labels.unique())}")
    logging.info(f"Predicted labels majority unique: {len(predicted_labels_majority.unique())}")
    
    logging.info(f"True labels overlap with predicted labels majority: {len(set(true_labels) & set(predicted_labels_majority))}")
    logging.info(f"True labels overlap with predicted labels: {len(set(true_labels) & set(predicted_labels))}")
    
    # print some statistics about the predicted labels and accuracy, not about the majority voting
    logging.info(f"Accuracy: {accuracy_score(true_labels, predicted_labels)}")
    
    
    # Create a confusion matrix
    cm = confusion_matrix(true_labels, predicted_labels)
    
    # Plot the confusion matrix
    plt.figure(figsize=(40, 40))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=predicted_labels.unique(), yticklabels=true_labels.unique())
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.show()
    plt.savefig("resources_data/celltypist/confusion_matrix.png")
    plt.close()
    
    
    
    
    
    
    





