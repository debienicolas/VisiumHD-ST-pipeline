#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip
import io
import os
import sys
import csv
# add the project root to the python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def plot_unique_transcript_positions(file_path, output_path='transcript_positions.png'):
    """
    Read transcript data from a gzipped TSV file and plot unique X,Y positions.
    
    Args:
        file_path: Path to the gzipped TSV file
        output_path: Path to save the output plot
    """
    print(f"Reading data from {file_path}...")
    
    # Read the gzipped TSV file
    with gzip.open(file_path, 'rt') as f:
        # Read header to determine column names
        header = f.readline().strip().split('\t')
        print(f"Found columns: {header}")
        
        # Find the correct column names for X and Y coordinates
        x_col = None
        y_col = None
        for col in header:
            if col.lower() == 'x':
                x_col = col
            elif col.lower() == 'y':
                y_col = col
        
        if not x_col or not y_col:
            raise ValueError(f"Could not find X and Y columns in the header: {header}")
            
        print(f"Using columns: {x_col} and {y_col} for coordinates")
        
        # Read data in chunks to avoid memory issues
        chunk_size = 100000
        unique_coords = set()
        total_lines = 0
        
        while True:
            chunk = []
            for _ in range(chunk_size):
                line = f.readline()
                if not line:
                    break
                chunk.append(line)
            
            if not chunk:
                break
                
            # Process chunk
            total_lines += len(chunk)
            if total_lines % 1000000 == 0:
                print(f"Processed {total_lines} lines...")
                
            # Parse the chunk manually to avoid pandas issues
            for line in chunk:
                parts = line.strip().split('\t')
                if len(parts) >= 2:  # Ensure we have at least X and Y columns
                    try:
                        x = float(parts[header.index(x_col)])
                        y = float(parts[header.index(y_col)])
                        unique_coords.add((x, y))
                    except (ValueError, IndexError):
                        continue  # Skip lines with parsing errors
    
    # Convert set to numpy array for plotting
    unique_coords_array = np.array(list(unique_coords))
    
    print(f"Found {len(unique_coords)} unique positions out of {total_lines} total transcripts")
    
    # Plot
    plt.figure(figsize=(10, 10))
    plt.scatter(unique_coords_array[:, 0], unique_coords_array[:, 1], s=1, alpha=0.5)
    plt.title('Unique Transcript Positions')
    plt.xlabel(f'{x_col} Position')
    plt.ylabel(f'{y_col} Position')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f'Plot saved as {output_path}')
    plt.show()

def plot_csv_transcript_positions(file_path, output_path='transcript_positions_csv.png'):
    """
    Read transcript data from a gzipped CSV file and plot unique X,Y positions.
    
    Args:
        file_path: Path to the gzipped CSV file
        output_path: Path to save the output plot
    """
    print(f"Reading CSV data from {file_path}...")
    
    # Read the gzipped CSV file
    with gzip.open(file_path, 'rt') as f:
        # Read header to determine column names
        csv_reader = csv.reader(f)
        header = next(csv_reader)
        print(f"Found columns: {header}")
        
        # Find the correct column names for X and Y coordinates
        x_col_idx = None
        y_col_idx = None
        
        for idx, col in enumerate(header):
            if col.lower() == 'x_location':
                x_col_idx = idx
            elif col.lower() == 'y_location':
                y_col_idx = idx
        
        if x_col_idx is None or y_col_idx is None:
            raise ValueError(f"Could not find x_location and y_location columns in the header: {header}")
            
        print(f"Using columns: {header[x_col_idx]} and {header[y_col_idx]} for coordinates")
        
        # Process data
        unique_coords = set()
        total_lines = 0
        
        for row in csv_reader:
            total_lines += 1
            if total_lines % 1000000 == 0:
                print(f"Processed {total_lines} lines...")
                
            try:
                x = float(row[x_col_idx])
                y = float(row[y_col_idx])
                unique_coords.add((x, y))
            except (ValueError, IndexError):
                continue  # Skip lines with parsing errors
    
    # Convert set to numpy array for plotting
    unique_coords_array = np.array(list(unique_coords))
    
    print(f"Found {len(unique_coords)} unique positions out of {total_lines} total transcripts")
    
    # Plot
    plt.figure(figsize=(10, 10))
    plt.scatter(unique_coords_array[:, 0], unique_coords_array[:, 1], s=1, alpha=0.5)
    plt.title('Unique Transcript Positions (CSV)')
    plt.xlabel(f'{header[x_col_idx]} Position')
    plt.ylabel(f'{header[y_col_idx]} Position')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f'Plot saved as {output_path}')
    plt.show()

def plot_parquet_transcript_positions(file_path, output_path='transcript_positions_parquet.png'):
    """
    Read transcript data from a Parquet file and plot unique X,Y positions.
    
    Args:
        file_path: Path to the Parquet file
        output_path: Path to save the output plot
    """
    print(f"Reading Parquet data from {file_path}...")
    
    try:
        # Check if pyarrow is installed
        import pyarrow.parquet as pq
        
        # Read the parquet file metadata to get column names
        parquet_file = pq.ParquetFile(file_path)
        schema = parquet_file.schema
        column_names = schema.names
        print(f"Found columns: {column_names}")
        
        # Find the correct column names for X and Y coordinates
        x_col = None
        y_col = None
        
        for col in column_names:
            col_lower = col.lower()
            if col_lower == 'x' or col_lower == 'x_location':
                x_col = col
            elif col_lower == 'y' or col_lower == 'y_location':
                y_col = col
        
        if not x_col or not y_col:
            raise ValueError(f"Could not find X and Y columns in the schema: {column_names}")
            
        print(f"Using columns: {x_col} and {y_col} for coordinates")
        
        # Process data in batches to avoid memory issues
        unique_coords = set()
        total_rows = 0
        
        # Read and process the file in chunks
        for batch in parquet_file.iter_batches(batch_size=100000, columns=[x_col, y_col]):
            df_batch = batch.to_pandas()
            total_rows += len(df_batch)
            
            if total_rows % 1000000 == 0:
                print(f"Processed {total_rows} rows...")
            
            # Add coordinates to the set
            for x, y in zip(df_batch[x_col], df_batch[y_col]):
                if pd.notna(x) and pd.notna(y):  # Skip NaN values
                    unique_coords.add((float(x), float(y)))
        
    except ImportError:
        print("PyArrow not found. Falling back to pandas for Parquet reading.")
        
        # Read with pandas (less memory efficient but more widely available)
        df = pd.read_parquet(file_path)
        print(f"Found columns: {df.columns.tolist()}")
        
        # Find the correct column names for X and Y coordinates
        x_col = None
        y_col = None
        
        for col in df.columns:
            col_lower = col.lower()
            if col_lower == 'x' or col_lower == 'x_location':
                x_col = col
            elif col_lower == 'y' or col_lower == 'y_location':
                y_col = col
        
        if not x_col or not y_col:
            raise ValueError(f"Could not find X and Y columns in the dataframe: {df.columns.tolist()}")
            
        print(f"Using columns: {x_col} and {y_col} for coordinates")
        
        # Extract unique coordinates
        coords = df[[x_col, y_col]].dropna()
        unique_coords = set(zip(coords[x_col].astype(float), coords[y_col].astype(float)))
        total_rows = len(df)
    
    # Convert set to numpy array for plotting
    unique_coords_array = np.array(list(unique_coords))
    
    print(f"Found {len(unique_coords)} unique positions out of {total_rows} total transcripts")
    
    # Plot
    plt.figure(figsize=(10, 10))
    plt.scatter(unique_coords_array[:, 0], unique_coords_array[:, 1], s=1, alpha=0.5)
    plt.title('Unique Transcript Positions (Parquet)')
    plt.xlabel(f'{x_col} Position')
    plt.ylabel(f'{y_col} Position')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f'Plot saved as {output_path}')
    plt.show()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot unique transcript positions from a file')
    parser.add_argument('--input', '-i', type=str, 
                        default='Xenium_data/output-XETG00447__0005530__A1__20250120__214610/ficture/input/transcripts.tsv.gz',
                        help='Path to the input file')
    parser.add_argument('--output', '-o', type=str, default='transcript_positions.png',
                        help='Path to save the output plot')
    parser.add_argument('--format', '-f', type=str, choices=['tsv', 'csv', 'parquet'], default='tsv',
                        help='Format of the input file (tsv, csv, or parquet)')
    
    args = parser.parse_args()
    
    if args.format == 'tsv':
        plot_unique_transcript_positions(args.input, args.output)
    elif args.format == 'csv':
        plot_csv_transcript_positions(args.input, args.output)
    else:  # parquet
        plot_parquet_transcript_positions(args.input, args.output)
