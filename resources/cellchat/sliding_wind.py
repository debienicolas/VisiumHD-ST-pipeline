import scanpy as sc
import numpy as np
from typing import Tuple, List
import pandas as pd

def calculate_window_parameters(
    coords: np.ndarray,
    n_windows_per_dim: int
) -> Tuple[float, float]:
    """
    Calculate window size and step size based on desired number of windows.
    
    Parameters:
    -----------
    coords : np.ndarray
        Array of spatial coordinates
    n_windows_per_dim : int
        Desired number of windows per dimension (x and y)
    
    Returns:
    --------
    tuple
        (window_size, step_size)
    """
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    
    # Calculate total span in each dimension
    x_span = x_max - x_min
    y_span = y_max - y_min
    
    # Use the larger span to determine window size
    max_span = max(x_span, y_span)
    
    # Calculate window size and step size
    window_size = max_span / (n_windows_per_dim - 1)  # Subtract 1 to ensure overlap ??? not really sure that this is necessary
    step_size = window_size / 2  # 50% overlap between windows
    
    return window_size, step_size

def create_spatial_windows(
    adata: sc.AnnData,
    n_windows_per_dim: int = 10,
    spatial_key: Tuple[str, int, int] = ('spatial', 0, 1),
) -> List[np.ndarray]:
    """
    Create sliding windows based on spatial coordinates in AnnData object.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with spatial coordinates
    n_windows_per_dim : int
        Number of windows per dimension (default: 10)
    spatial_key : tuple
        Keys for accessing spatial coordinates (default: ('spatial', 0, 1))
    min_cells : int
        Minimum number of cells required in a window
    
    Returns:
    --------
    list
        List of boolean masks for each window
    """
    print(f"\nCreating spatial windows with:")
    print(f"- Number of windows per dimension: {n_windows_per_dim}")
    print(f"- Spatial key: {spatial_key}")
    
    # Get spatial coordinates
    print("\nExtracting spatial coordinates...")
    coords = adata.obsm[spatial_key[0]]
    coords = coords[:, [spatial_key[1], spatial_key[2]]]
    print(f"Coordinate shape: {coords.shape}")
    
    # Calculate window parameters
    window_size, step_size = calculate_window_parameters(coords, n_windows_per_dim)
    print(f"\nCalculated parameters:")
    print(f"- Window size: {window_size:.2f}")
    print(f"- Step size: {step_size:.2f}")
    
    # Get coordinate ranges
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    print(f"Coordinate ranges:")
    print(f"X: {x_min:.2f} to {x_max:.2f}")
    print(f"Y: {y_min:.2f} to {y_max:.2f}")
    
    # Initialize list to store window masks
    windows = []
    
    # Slide window across x and y coordinates
    window_count = 0
    x_pos = x_min
    while x_pos <= x_max:
        y_pos = y_min
        while y_pos <= y_max:
            # Create window boundaries
            x_lower, x_upper = x_pos, x_pos + window_size
            y_lower, y_upper = y_pos, y_pos + window_size
            
            # Find cells in current window
            mask = (
                (coords[:, 0] >= x_lower) &
                (coords[:, 0] < x_upper) &
                (coords[:, 1] >= y_lower) &
                (coords[:, 1] < y_upper)
            )
            
            cells_in_window = np.sum(mask)
            windows.append(mask)
            window_count += 1
            if window_count % 10 == 0:
                print(f"Created window {window_count} at position ({x_pos:.2f}, {y_pos:.2f}) with {cells_in_window} cells")
        
            y_pos += step_size
        x_pos += step_size
    
    print(f"\nTotal windows created: {len(windows)}")
    return windows

def analyze_spatial_windows(
    adata: sc.AnnData,
    n_windows_per_dim: int = 10,
    spatial_key: str = 'spatial',
    min_cells: int = 10
) -> pd.DataFrame:
    """
    Perform sliding window analysis on spatial data.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with spatial coordinates
    n_windows_per_dim : int
        Number of windows per dimension (default: 10)
    spatial_key : str
        Key for accessing spatial coordinates in adata.obsm
    min_cells : int
        Minimum number of cells required in a window
    
    Returns:
    --------
    pd.DataFrame
        DataFrame containing statistics for each window
    """
    print("\nStarting spatial window analysis...")
    
    # Create windows
    print("Creating windows...")
    windows = create_spatial_windows(
        adata,
        n_windows_per_dim=n_windows_per_dim,
        spatial_key=(spatial_key, 0, 1),
    )
    
    # Initialize list to store window statistics
    window_stats = []
    print("\nAnalyzing windows...")
    
    # Analyze each window
    for idx, window_mask in enumerate(windows):
        window_cells = adata[window_mask]
        
        # Calculate basic statistics for the window
        stats = {
            'window_id': idx,
            'n_cells': len(window_cells),
            'x_center': window_cells.obsm[spatial_key][:, 0].mean(),
            'y_center': window_cells.obsm[spatial_key][:, 1].mean()
        }
        
        window_stats.append(stats)
        
        if idx % 10 == 0:
            print(f"Analyzed window {idx} containing {stats['n_cells']} cells")
    
    print(f"\nCompleted analysis of {len(window_stats)} windows")
    return pd.DataFrame(window_stats)

if __name__ == "__main__":
    print("Loading AnnData file...")
    adata = sc.read("results/run_4/ficture/annotated_anndata.h5ad")
    print(f"Loaded AnnData with shape: {adata.shape}")
    
    # Perform sliding window analysis
    window_stats = analyze_spatial_windows(
        adata,
        n_windows_per_dim=4,  # This will create approximately 100 windows (10x10 grid)
        spatial_key='spatial',
    )
    
    print("\nWindow Statistics Summary:")
    print(window_stats.describe())
    print("\nFirst few rows of window statistics:")
    print(window_stats.head())
