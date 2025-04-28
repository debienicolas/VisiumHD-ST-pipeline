import pandas as pd
import dask.dataframe as dd
import numpy as np

def aggregate_spot_counts(input_file, output_file, chunksize=1000000):
    """
    Memory-efficient aggregation of gene counts for each spatial spot.
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
        chunksize (int): Number of rows to process at once
    """
    # Use dtype optimization for memory efficiency
    dtypes = {
        'random_index': 'category',
        'X': np.float32,
        'Y': np.float32,
        'Count': np.int32
    }
    
    # Initialize empty aggregation dictionary
    spot_data = {}
    
    # Process the file in chunks
    for chunk in pd.read_csv(input_file, chunksize=chunksize, dtype=dtypes):
        # Group by random_index for this chunk
        chunk_agg = chunk.groupby('random_index').agg({
            'X': 'first',
            'Y': 'first',
            'Count': 'sum'
        })
        
        # Update our running totals
        for idx, row in chunk_agg.iterrows():
            if idx in spot_data:
                spot_data[idx]['Count'] += row['Count']
            else:
                spot_data[idx] = {
                    'X': row['X'],
                    'Y': row['Y'],
                    'Count': row['Count']
                }
    
    # Convert aggregated data to dataframe
    result = pd.DataFrame.from_dict(spot_data, orient='index')
    result.index.name = 'random_index'
    result.reset_index(inplace=True)
    
    # Save the result with efficient compression
    result.to_csv(output_file, index=False, compression='gzip')
    print(f"Aggregation complete. Output saved to {output_file}")
    print(f"Total spots processed: {len(result)}")
    print(f"Count range: {result['Count'].min()} - {result['Count'].max()}")

if __name__ == "__main__":
    input_file = "path_to_your_input.csv"
    output_file = "aggregated_counts.csv.gz"
    aggregate_spot_counts(input_file, output_file) 