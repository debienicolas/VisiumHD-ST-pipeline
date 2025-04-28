import argparse
import os
from loguru import logger
import re
import pandas as pd
from ficture.utils.utilt import make_mtx_from_dge # not to be confused with a relative import
import pickle
from sklearn.decomposition import LatentDirichletAllocation as LDA
from kneed import KneeLocator
import matplotlib.pyplot as plt



def ficture_coherence(input_path, output_path):
    """
    Check the perplexity scores for each factor and training width pair present in the input_path.
    Return a tsv file with the ranking.
    
    args:
        input_path: path to the ficture analysis folder.
        output_path: path to the output tsv file.
    """
    # declare results dataframe
    results = []
    
    # get all the folders in the input_path
    folders = [f for f in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, f))]
    logger.info(f"Found {len(folders)} folders in {input_path}")
    
    # use regex to get the n_factors and train_width from the folder names
    n_factors = list(set([re.search(r'nF(\d+)', f).group(1) for f in folders])) # type: ignore
    train_width = list(set([re.search(r'd_(\d+)', f).group(1) for f in folders])) # type: ignore
    logger.info(f"Found {len(n_factors)} n_factors and {len(train_width)} train_widths")
    
    # go to parent folder of input_path
    ficture_output_folder = os.path.dirname(input_path)
    
    # loop over all the training_widths
    for tw in train_width:
        logger.info(f"Processing training width {tw}")
        # load the sparse matrix for the current training_width
        _, _, mtx_org, _, _ = make_mtx_from_dge( # type: ignore
            os.path.join(ficture_output_folder, f"hexagon.d_{tw}.tsv.gz"),
            min_ct_per_feature = 20,
            min_ct_per_unit = 50,
            feature_list = None,
            unit = "random_index",
            key = "Count"
        ) 
        logger.info(f"Loaded sparse matrix for training width {tw}")
        for nF in n_factors:
            logger.info(f"Processing n_factors {nF}")
            # load the model
            with open(os.path.join(input_path, f"nF{nF}.d_{tw}", f"nF{nF}.d_{tw}.model.p"), "rb") as f:
                model = pickle.load(f)
            
            # calculat the perplexity score
            score = model.perplexity(mtx_org)
            logger.info(f"Calculated perplexity score for n_factors {nF} and training width {tw}")
            # append the results
            results.append({
                "id": f"nF{nF}.d_{tw}",
                "n_factors": nF,
                "train_width": tw,
                "perplexity": score
            })
            logger.info(f"Appended results for n_factors {nF} and training width {tw}")

        
    # create dataframe from results
    results_df = pd.DataFrame(results)
    
    # sort the dataframe by perplexity ascending
    results_df = results_df.sort_values(by="perplexity", ascending=True)
        
    # save the results to a tsv file
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Saved results to {output_path}")
    
    #### elbow plot ####
    # find the elbow point using the kneed package
    min_perplexity_tw = results_df.groupby('train_width')['perplexity'].min().idxmin()
    
    # if min_perplexity_tw has a single row than stop here
    if len(results_df[results_df['train_width'] == min_perplexity_tw]) == 1:
        logger.info(f"Only one row for train_width {min_perplexity_tw}, stopping here")
        # make the chosen_id.txt file
        with open(output_path.replace("coherence.tsv", "chosen_id.txt"), "w") as f:
            f.write(results_df[results_df['train_width'] == min_perplexity_tw]["id"].values[0])
        logger.info(f"Saved chosen_id to {output_path.replace('coherence.tsv', 'chosen_id.txt')}")
        return
    
    results_df = results_df[results_df['train_width'] == min_perplexity_tw]
    # turn the n_factors column into an integer
    results_df["n_factors"] = results_df["n_factors"].astype(int)
    results_df["perplexity"] = results_df["perplexity"].astype(float)
    # sort by n_factors ascending
    results_df = results_df.sort_values(by="n_factors", ascending=True)
    logger.info(f"results_df: {results_df}")
    logger.info("Finding elbow point")
    
    # Convert n_factors to integers before using them
    x = list(results_df["n_factors"].values)
    logger.info(f"x: {x}")
    y = list(results_df["perplexity"].values)
    logger.info(f"y: {y}")
    kl = KneeLocator(x, y, S=1.0, curve="convex", direction="decreasing")
    elbow_point = kl.elbow
    logger.info(f"Elbow point: {elbow_point}")
    
    # plot the elbow plot and save it to a png file
    logger.info("Plotting elbow plot")
    kl.plot_knee()
    plt.savefig(output_path.replace(".tsv", "_kneedle.png"))
    logger.info(f"Saved elbow plot to {output_path.replace('.tsv', '_kneedle.png')}")
    
    # write the elbow point to chosen_id.txt
    with open(output_path.replace("coherence.tsv", "chosen_id.txt"), "w") as f:
        f.write(results_df[results_df["n_factors"] == elbow_point]["id"].values[0])
    logger.info(f"Saved chosen_id to {output_path.replace('coherence.tsv', 'chosen_id.txt')}")
    
    return

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", type=str, default=None, required=False)
    parser.add_argument("--output_path", type=str, default=None, required=False)
    args = parser.parse_args()
    
    if args.input_path is None or args.output_path is None:
        results_dir = "human_breast_cancer"
        args.input_path = f"results/{results_dir}/ficture/ficture_output/analysis"
        args.output_path = f"results/{results_dir}/ficture/coherence.tsv"
    
    ficture_coherence(args.input_path, args.output_path)
    
    
    