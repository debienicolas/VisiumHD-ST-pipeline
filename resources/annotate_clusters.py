"""
Use this script for manual annotattion of the clusters.

Currenlty using claude 3.5 sonnet to annotate the clusters.
We can integrate an api call later, currenlty just 
"""
import argparse, os
import scanpy as sc
import pandas as pd
import anndata as ad
from dotenv import load_dotenv
from loguru import logger as logging
from pathlib import Path
import anthropic
import celltypist
import anndata as ad
import numpy as np
import scanpy as sc

# add parent directory to path
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def manual_annotation(input_file, ficture_output_folder, chosen_id, output_file):
    
    # load the anndata object
    adata = sc.read(input_file)
    logging.info(f"Loaded anndata object from {input_file}")
    
    # construct the cluster info path based on ficture_output_folder and chosen_id
    with open(chosen_id, 'r') as file:
        chosen_id = file.read()
    cluster_info_path = Path(ficture_output_folder) / "analysis" / chosen_id
    # find the file ending in .factor.info.tsv
    cluster_info_path = next(Path(cluster_info_path).glob('*.factor.info.tsv'))

    # load the anndata object
    adata = sc.read(input_file)
    print("adata", adata)
    print("adata.obs", adata.obs)
    print("adata.var", adata.var)
    #print("adata counts", adata.X.toarray().sum())
    # get the clusters
    factors = adata.obs["factor"]

    # take a , seperated string input from the user
    user_input = input("Enter the cluster info: ")
    # split the input into a list
    user_input = user_input.split(",")
    mapping = {i: user_input[i] for i in range(len(user_input))}

    # remap the factors to the cluster info
    adata.obs["factor"] = factors.map(mapping)
    
    #print("adata counts", adata.X.toarray().sum())

    # save the anndata object
    adata.write(output_file, compression="gzip")
    
    adata = ad.read_h5ad(output_file)
    #print("adata counts", adata.X.toarray().sum())
    
    # load the tsv file
    cluster_info = pd.read_csv(cluster_info_path, sep="\t")
    logging.info(f"Loaded cluster info from {cluster_info_path}")
    ## add functionality to add the annotation to the html and tsv files
    annotation_ordered = cluster_info["Factor"].map(mapping).to_list()
    # update the tsv file
    cluster_info["annotation"] = cluster_info["Factor"].map(mapping)
    cluster_info.to_csv(cluster_info_path, sep="\t", index=False)
    
    # update the html file
    from bs4 import BeautifulSoup
    cluster_info_html = cluster_info_path.with_suffix("").with_suffix(".info.html")
    with open(cluster_info_html, "r") as file:
        soup = BeautifulSoup(file, "html.parser")
        
    # Add new header column
    header_row = soup.find('thead').find('tr')
    new_header = soup.new_tag('th')
    new_header.string = 'Annotation'
    header_row.append(new_header)

    # Add new data cells to each row
    for row, data in zip(soup.find('tbody').find_all('tr'), annotation_ordered):
        new_cell = soup.new_tag('td')
        new_cell.string = str(data)  # Convert data to string before setting
        row.append(new_cell)

    # Save the updated HTML file
    with open(cluster_info_html, 'w') as file:
        file.write(soup.prettify())


    
    return None


def claude_annotation(input_file, ficture_output_folder, chosen_id, output_file):
    """
    Use this function to annotate the clusters using claude 3.5 sonnet.

    Args:
        input_file (str): The path to the input file: The clustered anndata object
        cluster_info_path (str): The path to the cluster info file: The file that contains the cluster info (.tsv)
        output_file (str): The path to the output file: The annotated anndata object
    """
    
    # load the environment variables: Claude api key should be in the .env file
    load_dotenv()
    
    # load the anndata object
    adata = ad.read_h5ad(input_file)
    logging.info(f"Loaded anndata object from {input_file}")
    logging.info(f"Anndata object: ")
    logging.info(adata)
    
    # construct the cluster info path based on ficture_output_folder and chosen_id
    with open(chosen_id, 'r') as file:
        chosen_id = file.read()
    cluster_info_path = Path(ficture_output_folder) / "analysis" / chosen_id
    # find the file ending in .factor.info.tsv
    cluster_info_path = next(Path(cluster_info_path).glob('*.factor.info.tsv'))
    
    # load the tsv file
    cluster_info = pd.read_csv(cluster_info_path, sep="\t")
    logging.info(f"Loaded cluster info from {cluster_info_path}")
    # turn this cluster info into a string
    cluster_info_string = cluster_info.to_string(index=False)
    logging.info(f"Cluster info string: {cluster_info_string}")
    
    # set up the claude prompt to annotate the clusters
    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))
    
    message = client.messages.create(
        model="claude-3-5-sonnet-20241022",
        max_tokens=4096,
        temperature=0,
        system="You are a world-class cell biology expert.",
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": f"""Below are the results of a spatial transcriptomics clustering analysis. \
                                The clusters are called factors. Please annotate the clusters based on the information provided. \
                                If clusters are the same type, add an incrementing number to the end of the cluster name. \
                                Only return a comma seperated list of the cluster names in ascending factor number order. \
                                The cluster info is:\
                                {cluster_info_string}
                                """
                    }
                ]
            }
        ]
    )
        
    annotated_clusters = message.content[0].text
    logging.info(f"Annotated clusters: {annotated_clusters}")
    
    # split the input into a list
    annotated_clusters = annotated_clusters.split(",")
    # remove any spaces from the annotated clusters
    annotated_clusters = [cluster.strip() for cluster in annotated_clusters]
    # replace any spaces in the annotated clusters with underscores
    annotated_clusters = [cluster.replace(" ", "_") for cluster in annotated_clusters]
    

    mapping = {i: annotated_clusters[i] for i in range(len(annotated_clusters))}

    # remap the factors to the cluster info
    factors = adata.obs["factor"]
    adata.obs["factor"] = factors.map(mapping)
    
    

    # save the anndata object
    adata.write(output_file, compression="gzip")
    
    
    ## add functionality to add the annotation to the html and tsv files
    annotation_ordered = cluster_info["Factor"].map(mapping).to_list()
    # update the tsv file
    cluster_info["annotation"] = cluster_info["Factor"].map(mapping)
    cluster_info.to_csv(cluster_info_path, sep="\t", index=False)
    
    # update the html file
    from bs4 import BeautifulSoup
    cluster_info_html = cluster_info_path.with_suffix("").with_suffix(".info.html")
    with open(cluster_info_html, "r") as file:
        soup = BeautifulSoup(file, "html.parser")
        
    # Add new header column
    header_row = soup.find('thead').find('tr')
    new_header = soup.new_tag('th')
    new_header.string = 'Annotation'
    header_row.append(new_header)

    # Add new data cells to each row
    for row, data in zip(soup.find('tbody').find_all('tr'), annotation_ordered):
        new_cell = soup.new_tag('td')
        new_cell.string = str(data)  # Convert data to string before setting
        row.append(new_cell)

    # Save the updated HTML file
    with open(cluster_info_html, 'w') as file:
        file.write(soup.prettify())

    return None


def atlas_annotation(input_file, ficture_output_folder, chosen_id, output_file, model_path:str="resources_data/celltypist/brca_atlas_model.pkl"):
    
    adata = ad.read_h5ad(input_file)
    
    # construct the cluster info path based on ficture_output_folder and chosen_id
    with open(chosen_id, 'r') as file:
        chosen_id = file.read()
    logging.info(f"Chosen id: {chosen_id}")
    cluster_info_path = Path(ficture_output_folder) / "analysis" / chosen_id
    # find the file ending in .factor.info.tsv
    cluster_info_path = next(Path(cluster_info_path).glob('*.factor.info.tsv'))
    
    # make sure that the counts are normalized to 10000 and log1p transformed
    if not np.allclose(1 + adata.X.toarray().sum(axis=1), 10000):
        logging.info("Normalizing the counts to 10000")
        sc.pp.normalize_total(adata, target_sum=10000)
        logging.info("Log1p transforming the counts")
        sc.pp.log1p(adata)
                
    predictions = celltypist.annotate(adata, model=model_path, majority_voting=True)
    result_adata = predictions.to_adata()
    
    print(result_adata)
    
    # the factor column contains the ficture factors (number)
    # create a original_factor column with this number and replace the factor column with the predicted labels
    
    cells_df = result_adata.obs
    
    #grouped_df = cells_df.groupby(["factor", "predicted_labels"])
    factor_grouped_df = cells_df.groupby("factor")
    mapping = {}
    # loof over factors and print the value counts of the predicted labels
    for factor, group in factor_grouped_df:
        logging.info(f"Factor: {factor}")
        logging.info(group["predicted_labels"].value_counts())
        
        majority_label = group["predicted_labels"].value_counts().idxmax()
        # if the majority label is already in the values of the mapping add the factor number to the end of the label
        mapping[factor] = f"{majority_label}_F{factor}"
        logging.info(f"Majority label: {majority_label}")
        
        # log the percentage of the majority label
        logging.info(f"Percentage of majority label: {group['predicted_labels'].value_counts()[majority_label] / len(group)}")

    
    result_adata.obs["original_factor"] = result_adata.obs["factor"]
    result_adata.obs["factor"] = result_adata.obs["factor"].map(mapping)
    
    result_adata.write(output_file, compression="gzip")
    
    
    logging.info(f"Unique factors:")
    logging.info(result_adata.obs["factor"].unique())
    
    return
    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=False, help="Path to the input file: The clustered anndata object")
    parser.add_argument("--output_file", type=str, required=False, help="Path to the output file: The annotated anndata object")
    parser.add_argument("--ficture_output_folder", type=str, required=False, help="Path to the ficture output folder: The folder that contains the ficture output")
    parser.add_argument("--chosen_id", type=str, required=False, help="Path to the chosen id file: The file that contains the chosen id")
    parser.add_argument("--man_annot", action="store_true", help="Boolean to indicate if manual annotation is used")
    args = parser.parse_args()
    
    if args.input_file is None or args.output_file is None:
        output_name = "human_DCIS_16um_correct"
        args.input_file = f"results/{output_name}/ficture/postprocessed_ficture/ficture_anndata.h5ad"
        args.output_file = f"results/{output_name}/ficture/annotated_anndata.h5ad"
        args.ficture_output_folder = f"results/{output_name}/ficture/ficture_output"
        args.chosen_id = f"results/{output_name}/ficture/chosen_id.txt"

    print("Manual annotation: ", args.man_annot)
    
    # This argument is ultimately set in the Snakefile
    if args.man_annot:
        manual_annotation(args.input_file,args.ficture_output_folder, args.chosen_id, args.output_file)
    else:
        #atlas_annotation(args.input_file, args.ficture_output_folder, args.chosen_id, args.output_file)
        claude_annotation(args.input_file, args.ficture_output_folder, args.chosen_id, args.output_file)



