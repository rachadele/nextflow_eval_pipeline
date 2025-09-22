#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
#import utils
#from utils import *
from pathlib import Path
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import seaborn as sns

import argparse
import os
import json

# Common optimizations for reading large files
read_opts = {
    "sep": "\t",
    "low_memory": False,          # Avoid dtype guessing in chunks
    "engine": "c",                # Faster C engine
    "na_filter": False,           # Skip NA parsing if not needed
    "memory_map": True            # Memory-map the file for faster IO
}

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["rachel_subclass", "rachel_class", "rachel_family"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff used for classification")
    parser.add_argument('--f1_results', type=str, nargs='+', help="Files containing F1 results")
    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    
def setup_plot(var, split):
    """Set up the basic plot structure."""
    plt.figure(figsize=(17, 8))
    plt.xlabel(split.split("_")[0].capitalize(), fontsize=25)
    var = var.replace("_", " ")
    #plt.ylabel(f"{var}", fontsize=25)
    plt.ylabel("Performance (weighted F1)", fontsize=25)
    #plt.title(f'Distribution of {var} across {split}', fontsize=25)


def add_violin_plot(df, var, split, facet):
    """Add a violin plot to the figure."""
    # remove extra weighted f1 values
    
    df = df.drop_duplicates(subset=[split, facet, var])
    sns.violinplot(
        data=df, 
        y=var, 
        x=split, 
        palette="Set2", 
        hue=facet, 
        split=False, 
        dodge=True
    )

def add_strip_plot(df, var, split, facet):
    """Add a strip plot to the figure."""
    # remove extra weighted f1 values
    # doesn't change overall values
    df = df.drop_duplicates(subset=[split, facet, var])
    df['match_region'] = df.apply(lambda row: row['query_region'] in row['ref_region'], axis=1)
    # Map match_region to colors before plotting
    df['color'] = df['match_region'].map({True: 'red', False: 'grey'})
    
    # Separate data into two groups based on 'match_region'
    mask = df['match_region']
    match_region_df = df[mask]
    non_match_region_df = df[~mask]
    
    # Create the strip plot for non-matching region data
    ax = sns.stripplot(
        data=non_match_region_df,
        y=var,
        x=split,
        hue=facet,
        dodge=True,          
        palette="Set2",      
        size=3,
        alpha=0.8,           
        jitter=True,
        marker="o",
        edgecolor='black',   
        linewidth=2
    )

    # Create the strip plot for matching region data with customized edgecolor
    sns.stripplot(
        data=match_region_df,
        y=var,
        x=split,
        hue=facet,
        dodge=True,          
        palette="Set2",      
        size=3,
        alpha=0.8,           
        jitter=True,
        marker="o",
        edgecolor='r',       # Red edge color for match region
        linewidth=2,         
        legend=None,         # Disable legend for second plot
        ax=ax                # Add to same axis
    )
     # Create custom legend handles for edge color
          

def add_acronym_legend(acronym_mapping, figure=None, x=1.05, y=0.5, title=None):
    if acronym_mapping:
        legend_text = f"{title}\n" + "\n".join([f"{k}: {v}" for k, v in acronym_mapping.items()])        
        figure = figure or plt.gcf()
        figure.text(
            x, y, legend_text,
            fontsize=14,
            verticalalignment='center',
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.1')
        )

def save_plot(var, split, facet, outdir):
    """Save the plot to the specified directory."""
    os.makedirs(outdir, exist_ok=True)
    var = var.replace(" ", "_")
    save_path = os.path.join(outdir, f"{var}_{split}_{facet}_distribution.png")
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()

def plot_distribution(df, var, outdir, split=None, facet="key", acronym_mapping=None):
    """
    Create a violin and strip plot for the given variable across groups.
    
    Parameters:
        df (pd.DataFrame): Data to plot.
        var (str): Column name for the variable to plot.
        outdir (str): Directory to save the plot.
        split (str): Column name to split the x-axis.
        facet (str): Column name for facet grouping (optional).
        acronym_mapping (dict): Mapping for acronyms to add as a legend (optional).
    """
    # sort df by split in alphabetical order
    df = df.sort_values(by=split)
    setup_plot(var, split)
    add_violin_plot(df, var, split, facet)
    add_strip_plot(df, var, split, facet)
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
    plt.xticks(rotation=90, ha="right", fontsize=25)
    plt.yticks(fontsize=25)
    red_patch = mlines.Line2D([], [], color='red', marker='o', markersize=7, label='Matching region')
    grey_patch = mlines.Line2D([], [], color='grey', marker='o', markersize=7, label='Non-Matching region')

   # Get current legend handles and labels (including the hue legend from violin plot)
    handles, labels = plt.gca().get_legend_handles_labels()
    
    # Add the custom patches to the legend
    handles.extend([red_patch, grey_patch])
    labels.extend(['Matching region', 'Non-Matching region'])
    
    # Update the legend with both the existing and custom legend items
    plt.legend(handles=handles, labels=labels, title="Match region", loc='upper left', bbox_to_anchor=(1, 1.02))
    #sns.move_legend(plt, bbox_to_anchor=(1, 1.02), loc='upper left')
    if acronym_mapping:
        add_acronym_legend(acronym_mapping, title=split.split("_")[0].capitalize())
    plt.tight_layout()
    save_plot(var, split, facet, outdir)
    
def make_acronym(ref_name):
    # Split on "_" and replace with spaces
    words = ref_name.split("_")
    # Create acronym from the first letter of each word
    acronym = "".join(word[0].upper() for word in words if word)
    return acronym
    
def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set organism and census_version from arguments
    ref_keys=args.ref_keys
    f1_results = args.f1_results
    cutoff=args.cutoff

    all_f1_scores = {}
    
    f1_df = pd.DataFrame()
    for file in f1_results:
        temp_df = pd.read_csv(file, **read_opts)
        f1_df = pd.concat([temp_df, f1_df], ignore_index=True)
    #need to flatten all_f1_scores into data frame
    f1_df["reference_acronym"] = f1_df["reference"].apply(make_acronym)
    f1_df["reference"] = f1_df["reference"].str.replace("_", " ")
 
    acronym_mapping_ref = f1_df[["reference", "reference_acronym"]].drop_duplicates().set_index("reference")["reference_acronym"].to_dict()
    
    plot_distribution(f1_df, var="weighted_f1",outdir="dists", split="reference_acronym", facet="key", 
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="weighted_f1",outdir="dists", split="study", facet="key", 
                      acronym_mapping = None)
    
    # plot NMI and ARI
    plot_distribution(f1_df, var="nmi", outdir="dists", split="reference_acronym", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="nmi", outdir="dists", split="study", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="ari", outdir="dists", split="reference_acronym", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="ari", outdir="dists", split="study", facet="key",
                      acronym_mapping = acronym_mapping_ref)


    # plot macro f1
    plot_distribution(f1_df, var="macro_f1", outdir="dists", split="reference_acronym", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="macro_f1", outdir="dists", split="study", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    
    
    # plot accuracy
    plot_distribution(f1_df, var="overall_accuracy", outdir="dists", split="reference_acronym", facet="key",
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="overall_accuracy", outdir="dists", split="study", facet="key",
                      acronym_mapping = acronym_mapping_ref)

if __name__ == "__main__":
    main()
    
