# standard library imports
from pathlib import Path
import os
import glob
import argparse
import subprocess 
import json
import zlib
import base64

import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
#from MulticoreTSNE import MulticoreTSNE as TSNE
from loguru import logger as logging
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap
from pyscenic.utils import load_motifs
import operator as op
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
import matplotlib as mpl

RESOURCES_DIR = "resources_data/pyscenic_dbs"

def create_motifs_table(motifs_file:str, output_file:str):

    BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
    COLUMN_NAME_LOGO = "MotifLogo"
    COLUMN_NAME_MOTIF_ID = "MotifID"
    COLUMN_NAME_TARGETS = "TargetGenes"

    def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
        """
        :param df:
        :param base_url:
        """
        # Make sure the original dataframe is not altered.
        df = df.copy()
        
        # Add column with URLs to sequence logo.
        def create_url(motif_id):
            return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
        df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
        
        # Truncate TargetGenes.
        def truncate(col_val):
            return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
        df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
        
        MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
        pd.set_option('display.max_colwidth', 200)
        print(df.head())
        pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
        return df

    df_motifs = load_motifs(motifs_file)
    df = display_logos(df_motifs)
    df.to_csv(output_file, index=True)
    # save the html render of this table to a file
    df.to_html(str(output_file).replace(".csv", ".html"), escape=False)
    return 


def rss_panel(auc_mtx, anndata, output_dir, annotation_column:str="Factor", heatmap_top_n:int=3):
    
    rss_celltype = regulon_specificity_scores(auc_mtx, anndata.obs[annotation_column])
    cats = list(anndata.obs[annotation_column].unique())
    
    fig = plt.figure(figsize=(15, 8))
    for c,num in zip(cats, range(1,len(cats)+1)):
        x=rss_celltype.T[c]
        ax = fig.add_subplot(2,5,num)
        plot_rss(rss_celltype, c, top_n=5, max_n=None, ax=ax)
        ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
        for t in ax.texts:
            t.set_fontsize(12)
        ax.set_ylabel('')
        ax.set_xlabel('')
        adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )

    fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
    fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
    plt.tight_layout()
    plt.rcParams.update({
        'figure.autolayout': True,
            'figure.titlesize': 'large' ,
            'axes.labelsize': 'medium',
            'axes.titlesize':'large',
            'xtick.labelsize':'medium',
            'ytick.labelsize':'medium'
            })
    plt.savefig(f"{output_dir}/rss_top5_panel.pdf", dpi=600, bbox_inches = "tight")
    plt.show()
    plt.close()
    
    ### Select the top_n regulons from each cell type
    topreg = []
    for i,c in enumerate(cats):
        topreg.extend(
            list(rss_celltype.T[c].sort_values(ascending=False)[:5].index)
        )
    topreg = list(set(topreg))

    # generate a z-score for each regulon to enable comparison between regulons
    auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
    for col in list(auc_mtx.columns):
        auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
        
    
    # generate a heatmap
    def palplot(pal, names, colors=None, size=1):
        n = len(pal)
        f, ax = plt.subplots(1, 1, figsize=(n * size, size))
        ax.imshow(np.arange(n).reshape(1, n),
                cmap=mpl.colors.ListedColormap(list(pal)),
                interpolation="nearest", aspect="auto")
        ax.set_xticks(np.arange(n) - .5)
        ax.set_yticks([-.5, .5])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        colors = n * ['k'] if colors is None else colors
        for idx, (name, color) in enumerate(zip(names, colors)):
            ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
        return f

    colors = sns.color_palette('bright',n_colors=len(cats) )
    colorsd = dict( zip( cats, colors ))
    colormap = [ colorsd[x] for x in anndata.obs[annotation_column] ]
    
    sns.set()
    sns.set(font_scale=0.8)
    fig = palplot( colors, cats, size=1.0)
    plt.savefig(f"{output_dir}/rss_heatmap_top{heatmap_top_n}_legend.pdf", dpi=600, bbox_inches = "tight")
    
    sns.set(font_scale=1.2)
    g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
        yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
        cmap="YlGnBu", figsize=(21,16) )
    g.cax.set_visible(True)
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xlabel('')
    plt.savefig(f"{output_dir}/rss_heatmap_top{heatmap_top_n}.pdf", dpi=600, bbox_inches = "tight")


def project_regulons_to_tissue(input_anndata_path, auc_mtx, cellAnnot, output_dir="results/human_breast_cancer_final/pyscenic", top_n=3):
    """
    Project regulon activity to spatial coordinates using the original input AnnData
    
    Parameters:
    -----------
    input_anndata_path : str
        Path to the original input AnnData file that contains spatial coordinates
    auc_mtx : DataFrame
        AUC matrix with cells as rows and regulons as columns from pySCENIC
    cellAnnot : DataFrame
        Cell annotations including cell type info in 'Factor' column
    output_dir : str
        Directory to save output plots
    top_n : int
        Number of top regulons per cell type to visualize
    """
    import scanpy as sc
    import matplotlib.pyplot as plt
    from pyscenic.rss import regulon_specificity_scores
    import pandas as pd
    import os
    import numpy as np
    
    # Load the original AnnData with spatial coordinates
    print(f"Loading spatial data from {input_anndata_path}")
    spatial_adata = sc.read_h5ad(input_anndata_path)
    
    # Calculate regulon specificity scores
    print("Calculating regulon specificity scores...")
    rss_celltype = regulon_specificity_scores(auc_mtx, cellAnnot["Factor"])
    
    # Get cell type categories
    cats = sorted(list(set(cellAnnot['Factor'])))
    
    # Get top regulons for each cell type
    top_regulons = []
    for cell_type in cats:
        top_regulons.extend(
            list(rss_celltype.T[cell_type].sort_values(ascending=False)[:top_n].index)
        )
    top_regulons = list(set(top_regulons))
    print(f"Top regulons selected: {', '.join(top_regulons)}")
    
    # Ensure the AnnData object and auc_mtx share the same cell IDs
    # Convert both index types to string for safer comparison
    spatial_cells = set(spatial_adata.obs.index.astype(str))
    auc_cells = set(auc_mtx.index.astype(str))
    common_cells = list(spatial_cells.intersection(auc_cells))
    
    if len(common_cells) == 0:
        print("No common cells between spatial data and regulon data")
        print(f"Spatial data has {len(spatial_cells)} cells with IDs like: {list(spatial_cells)[:5]}")
        print(f"AUC matrix has {len(auc_cells)} cells with IDs like: {list(auc_cells)[:5]}")
        return None
    
    print(f"Found {len(common_cells)} common cells between spatial data and regulon data")
    
    # Create a mapping dictionary if cell IDs in spatial_adata are different format but same content
    if len(common_cells) < min(len(spatial_cells), len(auc_cells)) * 0.5:
        print("Warning: Few common cells found. Attempting to match cell IDs...")
        # Try to match by removing common prefixes/suffixes or converting formats
        # For example, converting "cell_001" to "001" or vice versa
        
        # Example matching strategy - modify based on your specific ID formats
        spatial_to_auc = {}
        for s_id in spatial_cells:
            # Try various transformations of the ID
            stripped_id = s_id.split('-')[-1] if '-' in s_id else s_id
            if stripped_id in auc_cells:
                spatial_to_auc[s_id] = stripped_id
        
        if spatial_to_auc:
            print(f"Found {len(spatial_to_auc)} matched cells after ID transformation")
            
            # Create a modified spatial object with matched cells
            spatial_regulon = spatial_adata[list(spatial_to_auc.keys())].copy()
            
            # Add regulon AUC scores to the spatial data
            for regulon in top_regulons:
                if regulon in auc_mtx.columns:
                    spatial_regulon.obs[f'AUC_{regulon}'] = [
                        auc_mtx.loc[spatial_to_auc[cell_id], regulon] 
                        if spatial_to_auc[cell_id] in auc_mtx.index else np.nan
                        for cell_id in spatial_regulon.obs.index
                    ]
        else:
            print("Could not match cell IDs. Using only common cells.")
            spatial_regulon = spatial_adata[common_cells].copy()
            
            # Add regulon AUC scores to the spatial data
            for regulon in top_regulons:
                if regulon in auc_mtx.columns:
                    spatial_regulon.obs[f'AUC_{regulon}'] = auc_mtx.loc[common_cells, regulon].values
    else:
        # Use the directly common cells
        spatial_regulon = spatial_adata[common_cells].copy()
        
        # Add regulon AUC scores to the spatial data
        for regulon in top_regulons:
            if regulon in auc_mtx.columns:
                spatial_regulon.obs[f'AUC_{regulon}'] = auc_mtx.loc[common_cells, regulon].values
    
    # Check if the spatial coordinates exist
    spatial_keys = []
    for key in ['spatial', 'X_spatial']:
        if key in spatial_regulon.obsm:
            spatial_keys.append(key)
    
    if not spatial_keys:
        print("Warning: No spatial coordinates found in the AnnData object.")
        print("Available obsm keys:", list(spatial_regulon.obsm.keys()))
        
        # Create a directory for UMAP-based visualization instead
        os.makedirs(f"{output_dir}/regulon_projections", exist_ok=True)
        
        # Make sure UMAP coordinates exist
        if 'X_umap' in spatial_regulon.obsm:
            for regulon in top_regulons:
                if f'AUC_{regulon}' in spatial_regulon.obs.columns:
                    plt.figure(figsize=(10, 8))
                    sc.pl.umap(
                        spatial_regulon, 
                        color=f'AUC_{regulon}',
                        title=f'Regulon activity: {regulon}',
                        show=False
                    )
                    plt.savefig(f"{output_dir}/regulon_projections/{regulon}_umap.pdf", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
        return spatial_regulon
    
    # Plot spatial distribution of top regulons
    os.makedirs(f"{output_dir}/spatial_regulons", exist_ok=True)
    
    for regulon in top_regulons:
        if f'AUC_{regulon}' in spatial_regulon.obs.columns:
            plt.figure(figsize=(10, 8))
            
            # Use the first available spatial key
            spatial_key = spatial_keys[0]
            
            try:
                sc.pl.spatial(
                    spatial_regulon, 
                    color=f'AUC_{regulon}',
                    size=1.5,
                    title=f'Regulon activity: {regulon}',
                    show=False,
                    basis=spatial_key.replace('X_', '')  # Remove X_ prefix if present
                )
                plt.savefig(f"{output_dir}/spatial_regulons/{regulon}_spatial.pdf", 
                          dpi=300, bbox_inches='tight')
            except Exception as e:
                print(f"Error plotting spatial data for {regulon}: {e}")
                # Fallback to tissue image if available
                try:
                    sc.pl.spatial(
                        spatial_regulon,
                        img_key="hires",
                        color=f'AUC_{regulon}',
                        size=1.5,
                        title=f'Regulon activity: {regulon}',
                        show=False
                    )
                    plt.savefig(f"{output_dir}/spatial_regulons/{regulon}_spatial_hires.pdf", 
                              dpi=300, bbox_inches='tight')
                except Exception as e2:
                    print(f"Also failed with hires image: {e2}")
            
            plt.close()
    
    return spatial_regulon


def run_pyscenic(input_anndata:str, output_dir:str, species:str):
    os.makedirs(output_dir, exist_ok=True)
    
    ### Set the ouput files ###
    output_dir = Path(output_dir)
    # path to unfiltered loom file (this will be created in the optional steps below)
    f_loom_path_unfilt = output_dir / "pbmc10k_unfiltered.loom" # test dataset, n=500 cells
    # # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
    f_loom_path_scenic = output_dir / "pbmc10k_filtered_scenic.loom"
    # path to anndata object, which will be updated to store Scanpy results as they are generated below
    f_anndata_path = output_dir / "anndata.h5ad"
    # path to pyscenic output
    f_pyscenic_output = output_dir / "pyscenic_output.loom"
    # loom output, generated from a combination of Scanpy and pySCENIC results:
    f_final_loom = output_dir / 'pbmc10k_scenic_integrated-output.loom'
    
    ### Scanpy output settings ###
    sc.settings.verbosity = 3 # errors, warnings, info, hints
    sc.logging.print_versions()
    sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)
    sc.settings.njobs = os.environ.get("SLURM_CPUS_PER_TASK", 4) # extract this from slurm job if running on dcc
    sc.settings.njobs = 4
    ### Load the data from filtered_feature_bc_matrix ###
    adata = sc.read_h5ad(input_anndata)
    adata.var_names_make_unique()
    logging.info(f"Loaded data from {input_anndata}")
    
    
    ### write to the unfiltered loom file ### 
    row_attrs = {"Gene": np.array(adata.var.index)}
    # nUMI is the total counts per cell
    col_attrs = { "CellID":  np.array(adata.obs.index) ,"nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,"nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten()}
    lp.create( str(f_loom_path_unfilt), adata.X.transpose(), row_attrs, col_attrs )
    logging.info(f"Transpose X shape: {adata.X.transpose().shape}")
    logging.info(f"Written to {f_loom_path_unfilt}")
    
    
    #### for computational reasons sample 2000 cells ####
    adata = adata[np.random.choice(adata.n_obs, size=1_000, replace=False), :]
    
    ### log some attributes ###
    nCountsPerGene = np.sum(adata.X, axis=0)
    nCellsPerGene = np.sum(adata.X>0, axis=1)
    logging.info(f"Number of counts (in the dataset units) per gene: {nCountsPerGene.min()} - {nCountsPerGene.max()}")
    logging.info(f"Number of cells in which each gene is detected: {nCellsPerGene.min()} - {nCellsPerGene.max()}")
     

    
    ### pre-processing of expression data ###
    
    # save a copy of the raw data
    adata.raw = adata
    # Total-count normalize (library-size correct) to 10,000 reads/cell
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # log transform the data.
    sc.pp.log1p(adata)
    # identify highly variable genes.
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    # keep only highly variable genes:
    adata = adata[:, adata.var['highly_variable']]
    # regress out total counts per cell and the percentage of mitochondrial genes expressed
    #sc.pp.regress_out(adata, ['n_counts', 'percent_mito'] ) #, n_jobs=args.threads)
    # scale each gene to unit variance, clip values exceeding SD 10.
    sc.pp.scale(adata, max_value=10)
    # update the anndata file:
    adata.write(f_anndata_path)
    logging.info(f"Finished pre-processing of expression data")
    logging.info(f"Written to {f_anndata_path}")
    
    
    ### PCA and clustering ### 
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    adata.write( f_anndata_path )
    
    # visualization of highly variable genes
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    # tSNE
    # use regular TSNE
    sc.tl.tsne(adata, n_jobs=20)
    adata.write(f_anndata_path)
    
    # clustering
    sc.tl.leiden(adata)
    
    # find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    print(pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10))
    adata.write( f_anndata_path )
    
    # output the basic filtered expression matrix to loom file
    row_attrs = {"Gene": np.array(adata.var_names)}
    col_attrs = { "CellID":  np.array(adata.obs.index) ,"nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,"nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten()}
    lp.create( str(f_loom_path_scenic), adata.X.transpose(), row_attrs, col_attrs )
        
    
    #### SCENIC STEPS #### 
    
    ##### Gene regulatory network inference, and generation of co-expression 
    if species == "mouse":
        f_tfs_path = RESOURCES_DIR + "/mouse/allTFs_mm.txt"
    elif species == "human":
        f_tfs_path = RESOURCES_DIR + "/human/allTFs_hg38.txt"
    else:
        raise ValueError(f"Species {species} not supported")
    num_workers = os.environ.get("SLURM_CPUS_PER_TASK", 4)
    num_workers = 10
    logging.info(f"Number of workers: {num_workers}")
    output_adjacency_matrix = output_dir / "adjacency_matrix.tsv"
    
    logging.info(f"Subprocess echo the path: {subprocess.getoutput('pwd')}")
    try:
        result = subprocess.run(
            f"pyscenic grn {str(f_loom_path_scenic)} {str(f_tfs_path)} -o {str(output_adjacency_matrix)} --num_workers {num_workers}",
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        logging.info("GRN inference completed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"GRN inference failed with exit code {e.returncode}")
        logging.error(f"Error output: {e.stderr}")
        raise
        
    
    # adjacency matrix can be read as pandas dataframe
    
    ##### Regulon prediction aka cisTarget from cli 
    if species == "mouse":
        f_db_glob = RESOURCES_DIR + "/mouse/*feather"
        f_motif_path = RESOURCES_DIR + "/mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
    elif species == "human":
        f_db_glob = RESOURCES_DIR + "/human/*feather"
        f_motif_path = RESOURCES_DIR + "/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    else:
        raise ValueError(f"Species {species} not supported")
    f_db_names = ' '.join(glob.glob(f_db_glob))
    try:
        result = subprocess.run([
            "pyscenic", "ctx",
            str(output_dir/'adjacency_matrix.tsv'),
            *f_db_names.split(),  # Split the feather files into separate arguments
            "--annotations_fname", str(f_motif_path),
            "--expression_mtx_fname", str(f_loom_path_scenic),
            "--output", str(output_dir/'regulon.csv'),
            "--mask_dropouts",
            "--num_workers", str(num_workers)
        ], check=True, capture_output=True, text=True)
        logging.info("Regulon prediction completed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"Regulon prediction failed with exit code {e.returncode}")
        logging.error(f"Error output: {e.stderr}")
        raise
    
    
    ##### cellular enrichment
    try:
        result = subprocess.run([
            "pyscenic", "aucell",
            str(f_loom_path_scenic),
            str(output_dir/'regulon.csv'),
            "--output", str(f_pyscenic_output),
            "--num_workers", str(num_workers)
        ], check=True, capture_output=True, text=True)
        logging.info("Cellular enrichment completed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"Cellular enrichment failed with exit code {e.returncode}")
        logging.error(f"Error output: {e.stderr}")
        raise
    
    
    
def postprocess_scenic(f_pyscenic_output:str, anndata_path:str, input_anndata:str, output_dir:str, species:str):
    
    f_final_loom = Path(output_dir) / 'pbmc10k_scenic_integrated-output.loom'
    output_dir = Path(output_dir)
    ### Scenic postprocessing ### 
    adata = sc.read_h5ad(anndata_path)
    logging.info(f"Loaded data from {anndata_path}")
    logging.info(f"Anndata shape:")
    logging.info(adata)
    
    # load the input anndata
    adata_input = sc.read_h5ad(input_anndata)
    logging.info(f"Loaded data from {input_anndata}")
    logging.info(f"Input anndata shape:")
    logging.info(adata_input)
    
    # collect SCENIC AUCell output
    lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
    auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()
    
    # UMAP
    runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
    dr_umap = runUmap( auc_mtx )
    pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( output_dir / "scenic_umap.txt", sep='\t')
    # tSNE
    # use regular TSNE
    from sklearn.manifold import TSNE
    tsne = TSNE( n_jobs=20 )
    dr_tsne = tsne.fit_transform( auc_mtx )
    pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( output_dir / "scenic_tsne.txt", sep='\t')
        
        
    ### Integrate the output
    # Here we combine the results from SCENIC and Scanpy analysis into SCope compatible loom file
    # scenic output
    lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
    meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
    #exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
    auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons = lf.ra.Regulons
    dr_umap = pd.read_csv( output_dir / "scenic_umap.txt", sep='\t', header=0, index_col=0 )
    dr_tsne = pd.read_csv( output_dir / "scenic_tsne.txt", sep='\t', header=0, index_col=0 )
    ###
    
    # Fix regulon objects to display properly in SCope (not sure if necessary)
    auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
    regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
    # regulon thresholds
    rt = meta['regulonThresholds']
    for i,x in enumerate(rt):
        tmp = x.get('regulon').replace("(","_(")
        x.update( {'regulon': tmp} )
        
    logging.info("Regulons:")
    logging.info(regulons)
    logging.info("AUC matrix:")
    logging.info(auc_mtx)
    
    
    #### Postprocessing + pyscenic plots ####
    
    # Create motifs table
    regulon_file = Path(output_dir) / "regulon.csv"
    create_motifs_table(regulon_file, Path(output_dir) / "motifs.csv")
    
    # Plot the Rss panel
    rss_panel(auc_mtx, adata_input, output_dir, annotation_column="Factor", heatmap_top_n=3)
        
    ### concatenate embeddings (tSNE and UMAP)
    # print the type of the index of dr_umap
    logging.info(f"Type of dr_umap index: {dr_umap.index.dtype}")
    logging.info(f"Type of dr_tsne index: {dr_tsne.index.dtype}")
    logging.info(f"Type of adata.obs index: {adata.obs.index.dtype}")
    
    # turn all the indices of dr_umap, dr_tsne, and adata.obs to string
    dr_umap.index = dr_umap.index.astype(str)
    dr_tsne.index = dr_tsne.index.astype(str)
    adata.obs.index = adata.obs.index.astype(str)
    
    tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])
    Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
    Embeddings_X = pd.concat( [
            pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
            pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
            dr_tsne['X'] ,
            dr_umap['X']
        ], sort=False, axis=1, join='outer' )
    Embeddings_X.columns = ['1','2','3','4']
    
    logging.info(f"Embeddings_X created")
    logging.info(f"Embeddings_X head:")
    logging.info(Embeddings_X.head())
    logging.info(f"Embeddings_X shape: {Embeddings_X.shape}")

    Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
    Embeddings_Y = pd.concat( [
            pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
            pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
            dr_tsne['Y'] ,
            dr_umap['Y']
        ], sort=False, axis=1, join='outer' )
    Embeddings_Y.columns = ['1','2','3','4']
    
    # Check for NaN values in AUC matrix and replace with zeros
    nan_count = auc_mtx.isna().sum().sum()
    if nan_count > 0:
        logging.info(f"Found {nan_count} NaN values in AUC matrix, replacing with zeros")
        auc_mtx = auc_mtx.fillna(0)
    
    # concatenate embeddings (tSNE and UMAP)
    # tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

    # Fix: Convert named arrays to regular arrays for embeddings
    # umap_x = pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0].values
    # umap_y = pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1].values
    # pca_x = pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0].values
    # pca_y = pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1].values
    # tsne_x = dr_tsne['X'].values
    # tsne_y = dr_tsne['Y'].values
    # umap_scenic_x = dr_umap['X'].values
    # umap_scenic_y = dr_umap['Y'].values

    # metadata
    ### metadata
    metaJson = {}

    metaJson['embeddings'] = [
        {
            "id": -1,
            "name": f"Scanpy t-SNE (highly variable genes)"
        },
        {
            "id": 1,
            "name": f"Scanpy UMAP  (highly variable genes)"
        },
        {
            "id": 2,
            "name": "Scanpy PC1/PC2"
        },
        {
            "id": 3,
            "name": "SCENIC AUC t-SNE"
        },
        {
            "id": 4,
            "name": "SCENIC AUC UMAP"
        },
    ]

    metaJson["clusterings"] = [{
                "id": 0,
                "group": "Scanpy",
                "name": "Scanpy leiden default resolution",
                "clusters": [],
            },
            {
                "id": 1,
                "group": "ficture",
                "name": "ficture factors",
                "clusters": [],
            }
            ]

    metaJson["metrics"] = [
            {
                "name": "nUMI"
            }, {
                "name": "nGene"
            }
    ]

    metaJson["annotations"] = [
        {
            "name": "Leiden_clusters_Scanpy",
            "values": list(set(adata.obs['leiden'].astype(str)))
        },
        {
            "name": "clusters_ficture",
            "values": list(set(adata.obs['factor'].astype(str)))
        }
    ]

    # SCENIC regulon thresholds:
    metaJson["regulonThresholds"] = rt

    for i in range(max(set([int(x) for x in adata.obs['leiden']])) + 1):
        clustDict = {}
        clustDict['id'] = i
        clustDict['description'] = f'Unannotated Cluster {i + 1}'
        metaJson['clusterings'][0]['clusters'].append(clustDict)
    clusterings = pd.DataFrame()
    clusterings["0"] = adata.obs['leiden'].values.astype(np.int64)
    
    # Fix: Create a mapping of factor labels to numeric indices
    unique_factors = adata.obs['factor'].unique()
    factor_to_index = {factor: idx for idx, factor in enumerate(unique_factors)}
    factor_index_to_name = {idx: factor for idx, factor in enumerate(unique_factors)}
    sorted_factor_index_to_name = sorted(factor_index_to_name.items())
    metaJson["annotations"].append({"name": "factor_index_to_name", "values": [x[1] for x in sorted_factor_index_to_name]})
    
    # Create the clusters metadata
    for i, factor in enumerate(unique_factors):
        clustDict = {}
        clustDict['id'] = i
        clustDict['description'] = f'Factor {factor}'
        metaJson['clusterings'][1]['clusters'].append(clustDict)
    
    # Map the string factors to numeric indices
    clusterings["1"] = adata.obs['factor'].map(factor_to_index).values.astype(np.int64)
    
    # # Fix: Convert clusterings to a regular array
    # leiden_clusters = adata.obs['leiden'].values.astype(np.int64)
    

    # # Fix: Convert regulon data to standard arrays
    # # Extract regulon data as a regular numpy array
    # regulon_genes = {}
    # for name in regulons.dtype.names:
    #     regulon_genes[name] = regulons[name]

    # # Convert AUC matrix to a standard numpy array
    # auc_values = auc_mtx.values
    # auc_regulon_names = list(auc_mtx.columns)

    # Assemble loom file row and column attributes 
    def dfToNamedMatrix(df):
        arr_ip = [tuple(i) for i in df.values]
        dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
        arr = np.array(arr_ip, dtype=dtyp)
        return arr
    
    col_attrs = {
        "CellID": np.array(adata.obs.index),
        "nUMI": np.array(adata.obs['n_counts'].values),
        #"nGene": np.array(adata.obs['n_genes'].values),
        "Leiden_clusters_Scanpy": np.array(adata.obs['leiden'].values),
        "Factor": np.array(adata.obs['factor'].map(factor_to_index).values),
        "ClusterID": np.array(adata.obs['leiden'].values),
        # Fix: Store embeddings as separate arrays
        "Embedding": dfToNamedMatrix(tsneDF),
        "Embeddings_X": dfToNamedMatrix(Embeddings_X),
        "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
        "RegulonsAUC": dfToNamedMatrix(auc_mtx),
        "Clusterings": dfToNamedMatrix(clusterings),
        "ClusterID": np.array(adata.obs['factor'].values)
    }
    logging.info(f"Columns attributes created")

    # # Add AUC values as separate columns
    # for i, name in enumerate(auc_regulon_names):
    #     # Replace NaN values with 0 for each column
    #     col_data = auc_values[:, i]
    #     col_data = np.nan_to_num(col_data, nan=0.0)
    #     col_attrs[f"AUC_{name}"] = col_data

    row_attrs = {
        "Gene": lf.ra.Gene,
        "Regulons": regulons
    }
    logging.info(f"Row attributes created")

    # Add regulon genes as separate row attributes
    # for name, genes in regulon_genes.items():
    #     row_attrs[f"Regulon_{name}"] = genes

    if species == "mouse":
        genome = "mm10"
    elif species == "human":
        genome = "hg38"
    else:
        raise ValueError(f"Species {species} not supported")
    
    attrs = {
        "title": "sampleTitle",
        "MetaData": json.dumps(metaJson),
        "Genome": genome,
        "SCopeTreeL1": "",
        "SCopeTreeL2": "",
        "SCopeTreeL3": ""
    }

    # compress the metadata field:
    attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')
    
    lp.create(
        filename = str(f_final_loom),
        layers=lf[:,:],
        row_attrs=row_attrs, 
        col_attrs=col_attrs, 
        file_attrs=attrs
    )
    lf.close() # close original pyscenic loom file
    
    logging.info(f"Written to {f_final_loom}")
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_anndata", type=str, default=None, help="Path to the filtered feature matrix")
    parser.add_argument("--output_dir", type=str, default=None, help="Path to the output directory")
    parser.add_argument("--species", type=str, default=None, help="Species")
    args = parser.parse_args()
    
    if args.input_anndata is None or args.output_dir is None:
        output_dir = "run_4"
        output_dir = "human_breast_cancer_final"
        args.input_anndata = f"results/{output_dir}/ficture/annotated_anndata.h5ad"
        args.output_dir = f"results/{output_dir}/pyscenic"
        args.species = "mouse"
        args.species = "human"
        
    logging.info(f"Running pyscenic with {args.input_anndata} and {args.output_dir}")
    run_pyscenic(args.input_anndata, args.output_dir, args.species)
    postprocess_scenic(Path(args.output_dir) / "pyscenic_output.loom", Path(args.output_dir) / "anndata.h5ad", args.input_anndata, args.output_dir, args.species)
    
    
    
    