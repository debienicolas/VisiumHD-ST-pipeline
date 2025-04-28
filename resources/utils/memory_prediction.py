import numpy as np
import scipy.sparse as sp

def predict_anndata_memory(n_cells, n_genes, density=None):
    """
    Predict memory usage of AnnData object based on number of cells and genes.
    
    Parameters
    ----------
    n_cells : int
        Number of cells in the dataset
    n_genes : int
        Number of genes in the dataset
    density : float, optional
        Expected density of the count matrix (fraction of non-zero values).
        If None, both dense and sparse estimates are provided.
        If provided, only the appropriate format estimate is returned.
    
    Returns
    -------
    dict
        Memory estimates in bytes and human-readable format
    """
    # Basic memory usage for metadata structures
    base_memory = 1024 * 1024  # ~1MB for basic AnnData structure
    
    # Dense X matrix memory (float32 typically used)
    dense_x_memory = n_cells * n_genes * 4  # 4 bytes per float32 value
    
    # Sparse CSR memory calculation
    if density is None:
        # Default to assuming 10% density if not specified
        density = 0.1
        
    nnz = int(n_cells * n_genes * density)  # number of non-zero elements
    # CSR format: data array (nnz * 4), indices array (nnz * 4), indptr array ((n_cells+1) * 4)
    sparse_x_memory = (nnz * 4) + (nnz * 4) + ((n_cells + 1) * 4)
    
    # Total memory estimates
    dense_total = base_memory + dense_x_memory
    sparse_total = base_memory + sparse_x_memory
    
    # Convert to human-readable format
    def bytes_to_human_readable(bytes_val):
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if bytes_val < 1024 or unit == 'TB':
                return f"{bytes_val:.2f} {unit}"
            bytes_val /= 1024
    
    result = {
        "dense": {
            "bytes": dense_total,
            "human_readable": bytes_to_human_readable(dense_total)
        },
        "sparse": {
            "bytes": sparse_total,
            "human_readable": bytes_to_human_readable(sparse_total),
            "density": density
        }
    }
    
    # If density was explicitly provided, return only relevant format
    if density <= 0.5:
        print(f"With {density:.1%} density, sparse representation is recommended")
    else:
        print(f"With {density:.1%} density, dense representation is recommended")
        
    return result


def compare_memory_usage(n_cells, n_genes, density=0.1, measure_actual=False):
    """
    Compare memory usage between dense and sparse formats for different dataset sizes.
    
    Parameters
    ----------
    n_cells : int or list
        Number of cells (or list of cell counts to compare)
    n_genes : int or list
        Number of genes (or list of gene counts to compare)
    density : float
        Expected density of the count matrix
    measure_actual : bool
        Whether to also measure actual memory usage
        
    Returns
    -------
    dict
        Comparison of memory usage for different configurations
    """
    results = {}
    
    # Convert to lists if single values provided
    if isinstance(n_cells, int):
        n_cells = [n_cells]
    if isinstance(n_genes, int):
        n_genes = [n_genes]
    
    for cells in n_cells:
        for genes in n_genes:
            key = f"{cells} cells × {genes} genes"
            result = {
                "predicted": predict_anndata_memory(cells, genes, density)
            }
            if measure_actual:
                try:
                    result["actual"] = measure_actual_memory(cells, genes, density)
                except MemoryError:
                    result["actual"] = "Memory Error - dataset too large"
                except Exception as e:
                    result["actual"] = f"Error: {str(e)}"
            results[key] = result
    
    return results


def measure_actual_memory(n_cells, n_genes, density=0.1):
    """
    Create actual AnnData object and measure its memory usage.
    
    Parameters
    ----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    density : float
        Density of the matrix
        
    Returns
    -------
    dict
        Actual memory usage in bytes and human-readable format
    """
    import anndata as ad
    import psutil
    import os
    
    # Get initial memory usage
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss
    
    # Create sparse matrix with given density
    data = sp.random(n_cells, n_genes, density=density, format='csr', dtype=np.float32)
    adata_sparse = ad.AnnData(data)
    
    # Measure sparse memory
    sparse_memory = process.memory_info().rss - initial_memory
    
    # Create dense matrix
    data_dense = data.toarray()
    adata_dense = ad.AnnData(data_dense)
    
    # Measure dense memory
    dense_memory = process.memory_info().rss - initial_memory
    
    # Convert to human readable
    def bytes_to_human_readable(bytes_val):
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if bytes_val < 1024 or unit == 'TB':
                return f"{bytes_val:.2f} {unit}"
            bytes_val /= 1024
    
    return {
        "dense": {
            "bytes": dense_memory,
            "human_readable": bytes_to_human_readable(dense_memory)
        },
        "sparse": {
            "bytes": sparse_memory,
            "human_readable": bytes_to_human_readable(sparse_memory)
        }
    }


if __name__ == "__main__":
    # Example usage with verification
    print("Memory prediction vs actual usage for a small dataset:")
    n_cells, n_genes = 1000, 20000
    predicted = predict_anndata_memory(n_cells, n_genes)
    actual = measure_actual_memory(n_cells, n_genes)
    
    print("\nPredicted memory usage:")
    print(f"Dense matrix: {predicted['dense']['human_readable']}")
    print(f"Sparse matrix (density={predicted['sparse']['density']:.1%}): {predicted['sparse']['human_readable']}")
    
    print("\nActual memory usage:")
    print(f"Dense matrix: {actual['dense']['human_readable']}")
    print(f"Sparse matrix: {actual['sparse']['human_readable']}")
    
    print("\nComparison of different dataset sizes (with actual measurements):")
    comparison = compare_memory_usage(
        n_cells=[1000, 10_000, 100_000],
        n_genes=[10_000, 20_000],
        density=0.05,
        measure_actual=True
    )
    for config, result in comparison.items():
        print(f"\n{config}:")
        print("  Predicted:")
        print(f"    Dense: {result['predicted']['dense']['human_readable']}")
        print(f"    Sparse: {result['predicted']['sparse']['human_readable']}")
        if isinstance(result['actual'], dict):
            print("  Actual:")
            print(f"    Dense: {result['actual']['dense']['human_readable']}")
            print(f"    Sparse: {result['actual']['sparse']['human_readable']}")
        else:
            print(f"  Actual: {result['actual']}") 