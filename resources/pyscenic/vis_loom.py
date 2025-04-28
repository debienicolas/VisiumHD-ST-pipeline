import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Path to your loom file
loom_file = "results/run_4/pyscenic/pbmc10k_scenic_integrated-output.loom"

# Connect to the loom file
try:
    ds = lp.connect(loom_file, mode='r')
    print("Successfully connected to loom file")
    
    # Basic info about the file
    print(f"Number of cells: {ds.shape[1]}")
    print(f"Number of genes: {ds.shape[0]}")
    
    # Extract cell metadata and embeddings
    cell_ids = ds.ca.CellID
    
    # Print all column attributes for debugging
    print("\nAll column attributes:")
    for attr in ds.ca.keys():
        attr_data = ds.ca[attr]
        print(f"  - {attr}: shape={attr_data.shape}, type={type(attr_data)}, dtype={attr_data.dtype}")
        # Print a sample value if possible
        if len(attr_data) > 0:
            try:
                print(f"    Sample value: {attr_data[0]}")
            except:
                print("    Cannot display sample value")
    
    # Extract UMAP coordinates
    if 'UMAP_1' in ds.ca and 'UMAP_2' in ds.ca:
        umap_x = ds.ca.UMAP_1
        umap_y = ds.ca.UMAP_2
    else:
        print("\nUMAP coordinates not found with expected names, trying alternatives...")
        # Try to find any attributes that might contain UMAP coordinates
        coord_attrs = [attr for attr in ds.ca.keys() if any(dim in attr for dim in ['UMAP', 'umap', '_1', '_2'])]
        print(f"Potential coordinate attributes: {coord_attrs}")
        
        if len(coord_attrs) >= 2:
            # Try to pair them logically
            x_coords = [attr for attr in coord_attrs if any(x in attr for x in ['_1', 'X', 'x'])]
            y_coords = [attr for attr in coord_attrs if any(y in attr for y in ['_2', 'Y', 'y'])]
            
            if x_coords and y_coords:
                umap_x = ds.ca[x_coords[0]]
                umap_y = ds.ca[y_coords[0]]
                print(f"Using {x_coords[0]} and {y_coords[0]} for plotting")
            else:
                # Just use the first two coordinates found
                umap_x = ds.ca[coord_attrs[0]]
                umap_y = ds.ca[coord_attrs[1]]
                print(f"Using {coord_attrs[0]} and {coord_attrs[1]} for plotting")
        else:
            # If no coordinates found, create dummy coordinates
            print("No suitable coordinates found, creating dummy coordinates")
            umap_x = np.random.rand(ds.shape[1])
            umap_y = np.random.rand(ds.shape[1])
    
    # Extract cluster information
    if 'Leiden_clusters_Scanpy' in ds.ca:
        clusters = ds.ca.Leiden_clusters_Scanpy
    elif 'ClusterID' in ds.ca:
        clusters = ds.ca.ClusterID
    else:
        print("No clustering information found, using zeros")
        clusters = np.zeros(len(cell_ids))
    
    # Create a DataFrame for plotting
    plot_df = pd.DataFrame({
        'UMAP_1': umap_x,
        'UMAP_2': umap_y,
        'Cluster': clusters
    })
    
    # Plot UMAP colored by clusters
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=plot_df, x='UMAP_1', y='UMAP_2', hue='Cluster', 
                   palette='tab20', s=10, alpha=0.7)
    plt.title('UMAP Visualization of Cells')
    plt.savefig('umap_clusters.png', dpi=300)
    plt.close()
    print("Saved cluster visualization to umap_clusters.png")
    
    # Find regulon data - try multiple approaches
    print("\nLooking for regulon data...")
    
    # Approach 1: Look for AUC_ prefix
    auc_cols = [col for col in ds.ca.keys() if col.startswith('AUC_')]
    
    if auc_cols:
        print(f"Found {len(auc_cols)} potential regulon columns: {auc_cols[:5]}...")
        
        # Create a DataFrame with AUC values
        auc_data = {}
        for col in auc_cols:
            try:
                # Get the data and replace NaN values with 0
                data = ds.ca[col]
                data = np.nan_to_num(data, nan=0.0)  # Replace NaN with 0
                auc_data[col] = data
            except Exception as e:
                print(f"Error accessing column {col}: {e}")
        
        if auc_data:
            auc_mtx = pd.DataFrame(auc_data, index=cell_ids)
            
            # Check if we have any data with variance
            col_vars = auc_mtx.var()
            if (col_vars > 0.0001).any():
                # Get top variable regulons
                regulon_var = col_vars.sort_values(ascending=False)
                top_regulons = regulon_var.index[:6]  # Top 6 most variable regulons
                
                print(f"\nTop regulons by variance:")
                for i, reg in enumerate(top_regulons):
                    print(f"{i+1}. {reg}: variance = {regulon_var[reg]:.6f}")
                
                # Plot top regulons
                fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                axes = axes.flatten()
                
                for i, regulon in enumerate(top_regulons):
                    plot_df['Regulon_Activity'] = auc_mtx[regulon]
                    
                    # Print some stats about this regulon's values
                    values = plot_df['Regulon_Activity']
                    print(f"\nRegulon {regulon} stats:")
                    print(f"  Min: {values.min()}")
                    print(f"  Max: {values.max()}")
                    print(f"  Mean: {values.mean()}")
                    print(f"  Std: {values.std()}")
                    print(f"  Unique values: {len(values.unique())}")
                    print(f"  NaN values: {np.isnan(values).sum()}")
                    
                    # Check if we have enough variation to plot
                    if values.std() > 0.0001:
                        sns.scatterplot(data=plot_df, x='UMAP_1', y='UMAP_2', 
                                      hue='Regulon_Activity', palette='viridis', 
                                      s=10, alpha=0.7, ax=axes[i])
                        axes[i].set_title(f'Regulon: {regulon}')
                    else:
                        axes[i].text(0.5, 0.5, f"Insufficient variation in {regulon}", 
                                    ha='center', va='center', fontsize=12)
                        axes[i].set_title(f'Regulon: {regulon} (no variation)')
                    
                plt.tight_layout()
                plt.savefig('top_regulons.png', dpi=300)
                plt.close()
                print("Saved regulon visualization to top_regulons.png")
            else:
                print("No columns with sufficient variance found for plotting")
        else:
            print("Could not create AUC matrix from the columns")
    else:
        print("No regulon AUC scores found in the loom file")
    
    # Try to directly access the original AUC matrix from the pyscenic output
    print("\nTrying to directly visualize the original pyscenic output...")
    try:
        pyscenic_output = "results/run_4/pyscenic/pyscenic_output.loom"
        with lp.connect(pyscenic_output, mode='r') as lf:
            if 'RegulonsAUC' in lf.ca:
                auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
                print(f"Successfully loaded original AUC matrix with shape {auc_mtx.shape}")
                
                # Check for NaN values
                nan_count = auc_mtx.isna().sum().sum()
                print(f"Number of NaN values in original AUC matrix: {nan_count}")
                
                # Replace NaN with 0
                auc_mtx = auc_mtx.fillna(0)
                
                # Get top variable regulons
                regulon_var = auc_mtx.var().sort_values(ascending=False)
                top_regulons = regulon_var.index[:6]  # Top 6 most variable regulons
                
                # Create a new plot DataFrame
                # We need to match cell IDs between the two loom files
                common_cells = list(set(cell_ids).intersection(set(auc_mtx.index)))
                if len(common_cells) > 0:
                    print(f"Found {len(common_cells)} common cells between the two loom files")
                    
                    # Create a new DataFrame with UMAP coordinates and cluster info for common cells
                    common_idx = [list(cell_ids).index(cell) for cell in common_cells]
                    new_plot_df = pd.DataFrame({
                        'UMAP_1': umap_x[common_idx],
                        'UMAP_2': umap_y[common_idx],
                        'Cluster': clusters[common_idx]
                    }, index=common_cells)
                    
                    # Plot top regulons
                    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                    axes = axes.flatten()
                    
                    for i, regulon in enumerate(top_regulons):
                        new_plot_df['Regulon_Activity'] = auc_mtx.loc[common_cells, regulon]
                        
                        sns.scatterplot(data=new_plot_df, x='UMAP_1', y='UMAP_2', 
                                      hue='Regulon_Activity', palette='viridis', 
                                      s=10, alpha=0.7, ax=axes[i])
                        axes[i].set_title(f'Regulon: {regulon}')
                    
                    plt.tight_layout()
                    plt.savefig('original_top_regulons.png', dpi=300)
                    plt.close()
                    print("Saved original regulon visualization to original_top_regulons.png")
                else:
                    print("No common cells found between the two loom files")
            else:
                print("No RegulonsAUC found in original pyscenic output")
    except Exception as e:
        print(f"Error accessing original pyscenic output: {e}")
    
    # Close the connection
    ds.close()
    print("\nAnalysis complete")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()