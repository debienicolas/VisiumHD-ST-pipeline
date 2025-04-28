import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import anndata
import os
import warnings
from matplotlib.colors import LinearSegmentedColormap
import scipy.sparse
from loguru import logger as logging
import json

def analyze_spatial_transcriptomics(adata_path, label_column, atlas_path, atlas_label_column, output_dir='./results'):
    """
    Perform post-annotation analysis of spatial transcriptomics data.
    
    Parameters:
    -----------
    adata_path : str
        Path to the annotated AnnData object
    label_column : str
        Column name in adata.obs containing cell type labels
    output_dir : str
        Directory to save analysis results
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create figures directory inside output directory
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    # Load the annotated data
    print("Loading annotated data...")
    adata = sc.read_h5ad(adata_path)  # Using read_h5ad instead of read to avoid deprecation warning
    print(adata)
    
    # Verify label column exists
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found in data. Available columns: {list(adata.obs.columns)}")
    
    # Log-transform data if not already done (to avoid warning in rank_genes_groups)
    # if not adata.raw:
    #     print("Log-transforming data for differential expression analysis...")
    #     adata.raw = adata  # Store raw counts
    #     sc.pp.log1p(adata)
    
    # Basic QC and statistics of annotations
    #annotation_stats(adata, label_column, figures_dir, output_dir)
    
    # Identify marker genes for each cell type
    #find_marker_genes(adata, label_column, figures_dir, output_dir)
    
    # Spatial distribution of cell types
    #analyze_spatial_distribution(adata, label_column, figures_dir, output_dir)
    
    # Expression patterns of top markers
    #visualize_marker_expression(adata, label_column, figures_dir, output_dir)
    
    # Evaluate annotations
    eval_annotations(adata, atlas_path, atlas_label_column, label_column, output_dir)
    
    print(f"Analysis complete. Results saved to {output_dir}")

def compare_to_reference(adata:anndata.AnnData, atlas_path:str, atlas_label_column:str, label_column:str, output_dir:str):
    """
    Compare the annotations of the atlas and the annotated data.
    Do this by calculating cosine similarities between the cell type vectors of the atlas and the annotated data.
    save the results in a csv file.
    """
    # load the atlas
    atlas = sc.read_h5ad(atlas_path)
    print("loaded atlas")
    print(atlas)
    
    pass

def eval_annotations(adata:anndata.AnnData, atlas_path:str, atlas_label_column:str, label_column:str, output_dir:str):
    """Evaluate the annotations of the atlas and the annotated data."""
    print("Evaluating annotations...")
    
    # Load the atlas
    atlas = sc.read_h5ad(atlas_path)
    print("loaded atlas")
    print(atlas)
    # normalize the atlas to 10000
    sc.pp.normalize_total(atlas, inplace=True, target_sum=10_000)
    sc.pp.log1p(atlas)
    print("normalized atlas")
    print(atlas)
    
    # First, find common genes between the two datasets
    common_genes = list(set(adata.var_names).intersection(set(atlas.var_names)))
    print(f"Number of common genes: {len(common_genes)}")
    
    # Subset both datasets to common genes
    adata_subset = adata[:, common_genes]
    atlas_subset = atlas[:, common_genes]
    
    # Ensure matrices are in dense format for concatenation
    X_adata = adata_subset.X.toarray() if scipy.sparse.issparse(adata_subset.X) else adata_subset.X
    # print min and max of X_adata
    print(f"min of X_adata: {np.min(X_adata)}")
    print(f"max of X_adata: {np.max(X_adata)}")
    # do umap on the data
    sc.pp.neighbors(adata_subset)
    sc.tl.umap(adata_subset)
    sc.pl.umap(adata_subset, color=label_column, show=False)
    plt.savefig(f"{output_dir}/umap_adata.png")
    plt.close()
    print("umap done")
    
    
    X_atlas = atlas_subset.X.toarray() if scipy.sparse.issparse(atlas_subset.X) else atlas_subset.X
    # print min and max of X_atlas
    print(f"min of X_atlas: {np.min(X_atlas)}")
    print(f"max of X_atlas: {np.max(X_atlas)}")
    
    # do umap on the atlas
    sc.pp.neighbors(atlas_subset)
    sc.tl.umap(atlas_subset)
    sc.pl.umap(atlas_subset, color=atlas_label_column, show=False)
    plt.savefig(f"{output_dir}/umap_atlas.png")
    plt.close()
    print("umap done")

    # Create proper DataFrame for obs
    obs_df = pd.DataFrame({
        label_column: pd.concat([
            adata_subset.obs[label_column],
            atlas_subset.obs[atlas_label_column].rename(label_column)
        ], axis=0),
        'origin': np.concatenate([np.zeros(len(adata_subset.obs)), np.ones(len(atlas_subset.obs))])
    })
    
    # Create merged AnnData object
    adata_eval = anndata.AnnData(
        X=np.concatenate([X_adata, X_atlas], axis=0),
        var=adata_subset.var,
        obs=obs_df
    )
    
    print("Merged atlas and annotated data")
    print(adata_eval)
    
    # run umap on the data
    sc.pp.neighbors(adata_eval)
    sc.tl.umap(adata_eval)
    print("umap done")
    print(adata_eval)
    
    # plot the umap: colored by the cell type major and Factor column and the origin should be indicated by the shape of the points
    plt.figure(figsize=(10, 10))
    
    # Add UMAP coordinates to the observation DataFrame
    plot_df = adata_eval.obs.copy()
    plot_df['UMAP1'] = adata_eval.obsm['X_umap'][:, 0]
    plot_df['UMAP2'] = adata_eval.obsm['X_umap'][:, 1]
    
    # save the plot_df to a csv file
    plot_df.to_csv(f"{output_dir}/plot_df.csv", index=False)
    
    # Create the plot using the DataFrame
    sns.scatterplot(data=plot_df, x='UMAP1', y='UMAP2', hue=label_column, style='origin', size=100)
    plt.title('UMAP visualization of integrated data')
    plt.savefig(f"{output_dir}/umap_eval.png")
    plt.close()
    
def annotation_stats(adata, label_column, figures_dir, output_dir):
    """Calculate and visualize statistics about cell type annotations."""
    print("Calculating annotation statistics...")
    
    # Count cells per cell type
    cell_counts = adata.obs[label_column].value_counts()
    
    # Calculate percentage
    cell_percentages = 100 * cell_counts / len(adata)
    
    # Combine counts and percentages into a DataFrame
    stats_df = pd.DataFrame({
        'Cell Count': cell_counts,
        'Percentage': cell_percentages
    })
    
    # Save to CSV
    stats_df.to_csv(f"{output_dir}/cell_type_stats.csv")
    
    # Create bar plot with fixed ticks to avoid warning
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(x=stats_df.index, y='Cell Count', data=stats_df)
    plt.xticks(range(len(stats_df.index)), labels=stats_df.index, rotation=45, ha='right')
    plt.title(f'Number of Cells per {label_column}')
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/cell_type_counts.png", dpi=300)
    plt.close()
    
    # Calculate basic statistics for each cell type
    stats_dict = {
        'Total cells': len(adata),
        'Number of cell types': len(cell_counts),
        'Most abundant cell type': cell_counts.idxmax(),
        'Cells in most abundant type': cell_counts.max(),
        'Least abundant cell type': cell_counts.idxmin(),
        'Cells in least abundant type': cell_counts.min(),
    }
    
    # Save summary stats
    with open(f"{output_dir}/annotation_summary.txt", 'w') as f:
        for key, value in stats_dict.items():
            f.write(f"{key}: {value}\n")
    
    return stats_df

def find_marker_genes(adata:str|anndata.AnnData, label_column:str, figures_dir:str, output_dir:str):
    """Identify marker genes for each cell type and save results."""
    logging.info("Finding marker genes for each cell type...")
    
    if isinstance(adata, str):
        logging.info(f"Loading adata from file_path: {adata}")
        adata = sc.read_h5ad(adata)
    
    # Convert the label column to categorical type if it's not already
    logging.info(f"Label column data type: {adata.obs[label_column].dtype.name}")
    if adata.obs[label_column].dtype.name != 'str':
        logging.info(f"Converting {label_column} to string data type")
        adata.obs[label_column] = adata.obs[label_column].astype('str')
    
    # Get all cell types
    cell_types = adata.obs[label_column].unique().tolist()
    logging.info(f"Cell types found: {cell_types}")
    logging.info(f"Cell types data type: {type(cell_types[0])}")
    logging.info(f"cell types value counts: {adata.obs[label_column].value_counts()}")
    
    
    groups = adata.obs[label_column].unique().tolist()
    value_counts = adata.obs[label_column].value_counts()
    for group, count in value_counts.items():
        if count <= 1:
            groups.remove(group)
            logging.info(f"Removing cell type: {group} because it has only {count} cells")
    
    logging.info(f"Groups used for ranking: {groups}")
    
    # Use log-transformed data for rank_genes_groups
    sc.tl.rank_genes_groups(adata, groupby=label_column, groups=groups, method='wilcoxon', use_raw=False, pts=True)
    logging.info("Rank genes groups done")
    
    # Dictionary to store top markers for each cell type
    all_markers = {}
    
    # Extract results and make marker dataframes for each cell type
    for cell_type in groups:
        logging.info(f"Processing cell type: {cell_type}")
        logging.info("cell type data type: ", type(cell_type))
        # Get results for this cell type
        genes = adata.uns['rank_genes_groups']['names'][cell_type]
        scores = adata.uns['rank_genes_groups']['scores'][cell_type]
        pvals = adata.uns['rank_genes_groups']['pvals'][cell_type]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cell_type]
        logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cell_type]
        pts = adata.uns['rank_genes_groups']['pts'][cell_type]
        
        # Create DataFrame for this cell type
        markers_df = pd.DataFrame({
            'gene': genes,
            'score': scores,
            'logfoldchange': logfoldchanges,
            'pvalue': pvals,
            'pvalue_adj': pvals_adj,
            'pts': pts
        })
        
        # Sort by score
        markers_df = markers_df.sort_values('score', ascending=False)
        
        # Save markers for this cell type
        markers_df.to_csv(f"{output_dir}/markers_{cell_type.replace('/', '_')}.csv", index=False)
        logging.info(f"Markers for {cell_type} saved")
        # Store top 10 markers for visualization
        all_markers[cell_type] = markers_df.head(10)['gene'].tolist()
    
    # save all_markers to a json file
    with open(f"{output_dir}/all_markers.json", "w") as f:
        json.dump(all_markers, f)
    
    # Create a dotplot of top markers across cell types
    top_genes = []
    for cell_type, genes in all_markers.items():
        top_genes.extend(genes[:5])  # Take top 5 from each cell type
    
    # Remove duplicates but preserve order
    top_genes_unique = []
    for gene in top_genes:
        if gene not in top_genes_unique:
            top_genes_unique.append(gene)
    
    # Limit to top 50 genes to avoid overcrowded plots
    top_genes_unique = top_genes_unique[:50]
    
    # Generate dotplot with fixed filename
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.dotplot(adata, var_names=top_genes_unique, groupby=label_column, 
                    save=f"top_markers.png", show=False)
    # Move the file from scanpy's default location to our figures directory
    if os.path.exists('figures/dotplot_top_markers.png'):
        os.rename('figures/dotplot_top_markers.png', 
                 f"{figures_dir}/dotplot_top_markers.png")
    plt.close()
    logging.info("Dotplot top markers done")
    
    # Generate heatmap with fixed filename
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.heatmap(adata, var_names=top_genes_unique, groupby=label_column, 
                     save=f"top_markers.png", show=False)
    # Move the file from scanpy's default location to our figures directory
    if os.path.exists('figures/heatmap_top_markers.png'):
        os.rename('figures/heatmap_top_markers.png', 
                 f"{figures_dir}/heatmap_top_markers.png")
    plt.close()
    logging.info("Heatmap top markers done")
    return all_markers

def analyze_spatial_distribution(adata, label_column, figures_dir, output_dir):
    """Analyze the spatial distribution of cell types."""
    print("Analyzing spatial distribution of cell types...")
    
    # Check if spatial coordinates are available
    if 'spatial' not in adata.obsm.keys():
        print("Spatial coordinates not found. Skipping spatial analysis.")
        return
    
    # Create a spatial plot colored by cell type with fixed filename
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.spatial(adata, color=label_column, 
                     save=f"cell_types.png", show=False)
    # Move the file from scanpy's default location to our figures directory
    if os.path.exists(f'figures/spatial_cell_types.png'):
        os.rename(f'figures/spatial_cell_types.png', 
                 f"{figures_dir}/spatial_cell_types.png")
    
    # Calculate spatial statistics
    cell_types = adata.obs[label_column].unique()
    
    # Create a dataframe for spatial coordinates
    spatial_df = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'])
    spatial_df[label_column] = adata.obs[label_column].values
    
    # Calculate spatial statistics for each cell type
    spatial_stats = []
    for cell_type in cell_types:
        subset = spatial_df[spatial_df[label_column] == cell_type]
        
        # Calculate centroid
        centroid_x = subset['x'].mean()
        centroid_y = subset['y'].mean()
        
        # Calculate dispersion
        dispersion_x = subset['x'].std()
        dispersion_y = subset['y'].std()
        
        # Calculate nearest neighbor distances within cell type
        nn_distances = []
        if len(subset) > 1:
            from sklearn.neighbors import NearestNeighbors
            coords = subset[['x', 'y']].values
            nbrs = NearestNeighbors(n_neighbors=2).fit(coords)
            distances, _ = nbrs.kneighbors(coords)
            nn_distances = distances[:, 1]  # Distances to 1st neighbor (not self)
        
        spatial_stats.append({
            'cell_type': cell_type,
            'count': len(subset),
            'centroid_x': centroid_x,
            'centroid_y': centroid_y,
            'dispersion_x': dispersion_x,
            'dispersion_y': dispersion_y,
            'avg_nearest_neighbor': np.mean(nn_distances) if len(nn_distances) > 0 else np.nan
        })
    
    # Create spatial stats dataframe
    spatial_stats_df = pd.DataFrame(spatial_stats)
    spatial_stats_df.to_csv(f"{output_dir}/spatial_statistics.csv", index=False)
    
    # Create a plot showing cell type centroids with dispersion
    plt.figure(figsize=(10, 8))
    for _, row in spatial_stats_df.iterrows():
        plt.scatter(row['centroid_x'], row['centroid_y'], s=100, label=row['cell_type'])
        # Draw dispersion ellipse
        from matplotlib.patches import Ellipse
        ellipse = Ellipse((row['centroid_x'], row['centroid_y']), 
                          width=2*row['dispersion_x'], height=2*row['dispersion_y'],
                          alpha=0.2)
        plt.gca().add_patch(ellipse)
    
    plt.title(f'{label_column} Centroids and Dispersion')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/cell_type_centroids.png", dpi=300)
    plt.close()
    
    return spatial_stats_df

def visualize_marker_expression(adata, label_column, figures_dir, output_dir):
    """Visualize the expression patterns of top marker genes."""
    print("Visualizing expression patterns of marker genes...")
    
    # Create red-yellow heatmap colormap for spatial plots
    red_yellow_cmap = LinearSegmentedColormap.from_list('red_yellow', ['#FFFF00', '#FF0000'])
    
    # Use log-transformed data for rank_genes_groups
    sc.tl.rank_genes_groups(adata, groupby=label_column, method='wilcoxon', use_raw=False)
    
    cell_types = adata.obs[label_column].unique()
    
    # For each cell type, plot expression of top markers
    for cell_type in cell_types:
        # Get top 3 markers
        markers = adata.uns['rank_genes_groups']['names'][cell_type][:3]
        
        # If spatial data available, create spatial plots
        if 'spatial' in adata.obsm.keys():
            for gene in markers:
                if gene in adata.var_names:
                    # Use red-yellow colormap for spatial plots
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        sc.pl.spatial(adata, img=None, color=gene, cmap=red_yellow_cmap,
                                    save=f"{cell_type.replace('/', '_')}_{gene}_spatial.png", show=False)
                    # Move the file from scanpy's default location to our figures directory
                    if os.path.exists(f'figures/spatial_{cell_type.replace("/", "_")}_{gene}_spatial.png'):
                        os.rename(f'figures/spatial_{cell_type.replace("/", "_")}_{gene}_spatial.png', 
                                f"{figures_dir}/spatial_{cell_type.replace('/', '_')}_{gene}.png")
        
        # Create violin plots with updated parameters to avoid deprecation warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig, ax = plt.subplots(figsize=(10, 5))
            # Use updated violin plot params to avoid deprecation warnings
            sc.pl.violin(adata, markers, groupby=label_column, ax=ax,
                       save=f"{cell_type.replace('/', '_')}_markers.png", show=False)
        # Move the file from scanpy's default location to our figures directory
        if os.path.exists(f'figures/violin_{cell_type.replace("/", "_")}_markers.png'):
            os.rename(f'figures/violin_{cell_type.replace("/", "_")}_markers.png', 
                    f"{figures_dir}/violin_{cell_type.replace('/', '_')}_markers.png")
    
    # Create a combined feature plot of top markers for each cell type
    top_genes = {}
    for cell_type in cell_types:
        top_genes[cell_type] = adata.uns['rank_genes_groups']['names'][cell_type][0]
    
    # UMAP visualization if not already computed
    if 'X_umap' not in adata.obsm.keys():
        try:
            print("Computing PCA and UMAP for visualization...")
            # Perform PCA with specific number of components to avoid warning
            sc.pp.pca(adata, n_comps=50)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        except Exception as e:
            print(f"Could not compute UMAP: {e}")
    
    # If UMAP is available, plot it
    if 'X_umap' in adata.obsm.keys():
        # Plot UMAP colored by cell type
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pl.umap(adata, color=label_column, 
                      save=f"celltypes.png", show=False)
        # Move the file from scanpy's default location to our figures directory
        if os.path.exists('figures/umap_celltypes.png'):
            os.rename('figures/umap_celltypes.png', 
                    f"{figures_dir}/umap_celltypes.png")
        
        # Plot UMAP with top marker for each cell type
        for cell_type, gene in top_genes.items():
            if gene in adata.var_names:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sc.pl.umap(adata, color=gene, 
                             save=f"{cell_type.replace('/', '_')}_top_marker.png", show=False)
                # Move the file from scanpy's default location to our figures directory
                if os.path.exists(f'figures/umap_{cell_type.replace("/", "_")}_top_marker.png'):
                    os.rename(f'figures/umap_{cell_type.replace("/", "_")}_top_marker.png', 
                            f"{figures_dir}/umap_{cell_type.replace('/', '_')}_top_marker.png")

def fix_deprecated_plots():
    """Update seaborn violin plots to avoid deprecation warnings."""
    # This function patches scanpy's plotting functions to avoid deprecation warnings
    # It can be called before running the analysis
    
    original_violin = sc.pl.violin
    
    def updated_violin(*args, **kwargs):
        if 'scale' in kwargs:
            kwargs['density_norm'] = 'width'
            del kwargs['scale']
        return original_violin(*args, **kwargs)
    
    sc.pl.violin = updated_violin

if __name__ == "__main__":
    import argparse
    
    # Fix seaborn deprecation warnings
    fix_deprecated_plots()
    
    parser = argparse.ArgumentParser(description='Analyze annotated spatial transcriptomics data')
    parser.add_argument('--input', type=str, required=False, default=None, 
                        help='Path to annotated AnnData object (h5ad file)')
    parser.add_argument('--label_column', type=str, required=False, default="factor",
                        help='Column name in adata.obs containing cell type labels')
    parser.add_argument('--atlas_path', type=str, required=False, default="resources_data/BrCa_Atlas_Count_out/annotated_brca_atlas.h5ad",
                        help='Path to atlas file (h5ad file)')
    parser.add_argument('--atlas_label_column', type=str, required=False, default="celltype_major",
                        help='Column name in atlas.obs containing cell type labels')
    parser.add_argument('--output', type=str, default='./results',
                        help='Directory to save analysis results')
    
    args = parser.parse_args()
    if args.input is None:
        args.input = "celltypist_results/celltypist_annotated_visium_0.h5ad"
        args.label_column = "predicted_labels"
        args.atlas_label_column = "celltype_major"
        args.output = "results/human_breast_cancer_final/post_annotation"
    
    analyze_spatial_transcriptomics(args.input, args.label_column, args.atlas_path, args.atlas_label_column, args.output)