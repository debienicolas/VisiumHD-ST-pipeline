import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
import glob
import argparse

import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def create_decoded_images_grid(base_dir):
    # Find all directories matching the pattern nF{n_factors}.d_{train_width}
    dir_pattern = re.compile(r'nF(\d+)\.d_(\d+)')
    
    # Store information about each directory
    dirs_info = []
    
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            match = dir_pattern.match(item)
            if match:
                n_factors = int(match.group(1))
                train_width = int(match.group(2))
                
                # Look for the decoded image in this directory
                # The pattern might need adjustment based on actual file naming
                image_pattern = os.path.join(item_path, "figure", f"{item}.decode.prj_*.png")
                image_files = glob.glob(image_pattern)
                
                # If no image found with that pattern, try looking for any PNG file
                if not image_files:
                    image_files = glob.glob(os.path.join(item_path, "*.png"))
                
                if image_files:
                    # Use the first matching image
                    dirs_info.append({
                        'n_factors': n_factors,
                        'train_width': train_width,
                        'image_path': image_files[0]
                    })
    
    if not dirs_info:
        print("No matching directories or images found.")
        return
    
    # Get unique values for n_factors and train_width
    n_factors_values = sorted(set(info['n_factors'] for info in dirs_info))
    train_width_values = sorted(set(info['train_width'] for info in dirs_info), reverse=True)
    
    # Create a grid for the images
    fig = plt.figure(figsize=(len(n_factors_values) * 4, len(train_width_values) * 4))
    gs = GridSpec(len(train_width_values), len(n_factors_values))
    
    # Create a mapping for quick lookup
    image_map = {(info['n_factors'], info['train_width']): info['image_path'] for info in dirs_info}
    
    # Plot each image in its corresponding position
    for i, train_width in enumerate(train_width_values):
        for j, n_factors in enumerate(n_factors_values):
            ax = fig.add_subplot(gs[i, j])
            
            # Check if we have an image for this combination
            if (n_factors, train_width) in image_map:
                try:
                    img = Image.open(image_map[(n_factors, train_width)])
                    ax.imshow(np.array(img))
                except Exception as e:
                    ax.text(0.5, 0.5, f"Error loading image:\n{str(e)}", 
                            ha='center', va='center', transform=ax.transAxes)
            else:
                ax.text(0.5, 0.5, "No image", ha='center', va='center', transform=ax.transAxes)
            
            ax.set_title(f"nF{n_factors}.d_{train_width}")
            ax.axis('off')
    
    plt.tight_layout()
    # Save as PNG
    plt.savefig(os.path.join(base_dir, "decoded_images_grid.png"), dpi=300)
    # Save as high-quality PDF
    plt.savefig(os.path.join(base_dir, "decoded_images_grid.pdf"), format='pdf', dpi=300)
    plt.close()
    
    print(f"Grid images saved to:")
    print(f"  - {os.path.join(base_dir, 'decoded_images_grid.png')} (PNG)")
    print(f"  - {os.path.join(base_dir, 'decoded_images_grid.pdf')} (PDF)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a grid of decoded images')
    parser.add_argument('--base_dir', '-b', type=str,
                        help='Path to the base directory containing the decoded images (analysis directory in ficture output)')
    args = parser.parse_args()
    
    if args.base_dir is None:
        args.base_dir = "Xenium_data/output-XETG00447__00051697__B3__20250120__214611/ficture/output/analysis"
    
    create_decoded_images_grid(args.base_dir)
    
    