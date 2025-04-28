

import scanpy as sc
import squidpy as sq
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib 


def main(output_file:str):
    
    adata = ad.read_h5ad(output_file)
    
    print(adata)
    # provide a custom color palette that doesn't have pink
    cmap = matplotlib.colormaps["viridis"]
    cmap = "viridis"
    dpi = 400
    
    sq.pl.spatial_scatter(adata, color="cnv_leiden", img_alpha=0.5, cmap=cmap)
    plt.savefig("infercnvpy_output.png", dpi=dpi)
    sq.pl.spatial_scatter(adata, color="factor", img_alpha=0.5, cmap=cmap)
    plt.savefig("infercnvpy_output_factor.png", dpi=dpi)









if __name__ == "__main__":
    infercnv_output_file = "results/human_breast_cancer_final/infercnvpy/infercnvpy_output.h5ad"
    main(infercnv_output_file)