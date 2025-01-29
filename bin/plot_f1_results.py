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
#import adata_functions
#from adata_functions import *
from pathlib import Path
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import seaborn as sns

import argparse
import os
import json

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
    setup_plot(var, split)
    add_violin_plot(df, var, split, facet)
    add_strip_plot(df, var, split, facet)
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
    plt.xticks(rotation=90, ha="right", fontsize=25)
    plt.yticks(fontsize=25)
    red_patch = mlines.Line2D([], [], color='red', marker='o', markersize=7, label='Matching region')
    grey_patch = mlines.Line2D([], [], color='grey', marker='o', markersize=7, label='Non-Matching region')

    # Add the custom legend to the plot
    plt.legend(handles=[red_patch, grey_patch], title="Match region", loc='upper left', bbox_to_anchor=(1, 1.02))

    # Move the legend to the desired location
    #sns.move_legend(plt, bbox_to_anchor=(1, 1.02), loc='upper left')

    add_acronym_legend(acronym_mapping, title=split.split("_")[0].capitalize())
    plt.tight_layout()
    save_plot(var, split, facet, outdir)
    
    
def calculate_global_f1_range(all_f1_scores):
    """Determine the global F1 score range."""
    all_scores = pd.concat([df['f1_score'] for df in all_f1_scores.values()])
    return all_scores.min(), all_scores.max()

def prepare_queries_and_keys(all_f1_scores):
    """Extract unique queries and keys from the F1 scores dictionary."""
    queries = set()
    keys = list(all_f1_scores.keys())
    for key, df in all_f1_scores.items():
        queries.update(df['query'].unique())
    return sorted(queries), keys

def plot_heatmap(data, ax, vmin, vmax, show_cbar, mask=None, title=None, 
                 x_fontsize=20, y_fontsize=20, leftmost=False):
    sns.heatmap(
        data,
        annot=True,
        cmap='YlOrRd',
        cbar=show_cbar,
        cbar_kws={'label': 'F1 Score'} if show_cbar else None,
        mask=mask,
        ax=ax,
        linewidths=0.5,
        annot_kws={"size": 8},
        vmin=vmin,
        vmax=vmax,
    )
    if title:
        ax.set_title(title, fontsize=25)
    
    ax.set_xlabel('', fontsize=x_fontsize)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=x_fontsize)
    
    return ax



def save_heatmap(fig, outpath, filename):
    """Save the heatmap figure."""
    fig.savefig(os.path.join(outpath, filename), bbox_inches='tight')
    plt.close(fig)

def plot_label_f1_heatmaps(all_f1_scores, threshold, outpath, widths=[1, 0.8, 0.5], acronym_mapping=None, fontsize=20):
    """Main function to plot label-level F1 heatmaps."""
    sns.set(style="whitegrid")
    os.makedirs(outpath, exist_ok=True)

    vmin, vmax = calculate_global_f1_range(all_f1_scores)
    queries, keys = prepare_queries_and_keys(all_f1_scores)

    if widths is None:
        widths = [1] * len(keys)

    for query in queries:
        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(keys),
            figsize=(sum(widths) * 10, 8),
            gridspec_kw={'width_ratios': widths},
            constrained_layout=True
        )
        fig.suptitle(f'Class-level F1 for Query: {query}\nThreshold = {threshold:.2f}', fontsize=fontsize+5, y=1.1)

        for i, key in enumerate(keys):
            if key not in all_f1_scores:
                continue

            query_df = all_f1_scores[key][all_f1_scores[key]['query'] == query]
            pivot_df = query_df.pivot_table(index='reference_acronym', columns='label', values='f1_score')
            mask = pivot_df.isnull() | (pivot_df == "nan")

            plot_heatmap(pivot_df, axes[i], vmin, vmax, show_cbar=(i == len(keys) - 1), 
                         mask=mask, title=key)
            if i == 0:
                axes[i].set_ylabel('Reference', fontsize=fontsize)
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=fontsize, 
                                        rotation=0)
                if acronym_mapping:
                    legend_text = "\n".join([f"{k}: {v}" for k, v in acronym_mapping.items()])
                    axes[i].text(
                        0.1, 0.5, legend_text,  # Position text to the left of the subplot
                        fontsize=14,
                        verticalalignment='center',
                        horizontalalignment='right',
                        transform=axes[i].transAxes,  # Use subplot's coordinate system
                        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.1')
                    )
            else:
                axes[i].set_ylabel("", fontsize=11)
                axes[i].set_yticks([]) # Hide y-axis labels
        save_heatmap(fig, outpath, f'f1_heatmaps_{query}_threshold_{threshold:.2f}.png')
        
        
        

def plot_f1_heatmap_for_level(data, vmin, vmax, level, threshold, outpath, 
                              acronym_mapping=None, fontsize=25):
    """Plot a single level heatmap."""
    pivot_f1 = data.pivot_table(index='reference_acronym', columns='query', values='weighted_f1')

    plt.figure(figsize=(20, 15))
    sns.heatmap(
        pivot_f1,
        annot=True,
        cmap='YlOrRd',
        cbar_kws={'label': 'Weighted F1 Score'},
        fmt='.3f',
        annot_kws={"size": 15},
        vmin=vmin,
        vmax=vmax
    )
    plt.title(f'Weighted F1 Scores for {level}\nThreshold = {threshold:.2f}', fontsize=fontsize)
    plt.ylabel('Reference', fontsize=fontsize)
    plt.xlabel('Query', fontsize=fontsize)
    plt.xticks(rotation=90, ha='right', fontsize=fontsize)
    plt.yticks(fontsize=fontsize, rotation=45)

    if acronym_mapping:
        legend_text = "\n".join([f"{k}: {v}" for k, v in acronym_mapping.items()])
        plt.gcf().text(
            0.85, 0.5, legend_text,
            fontsize=15,
            verticalalignment='center',
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.1')
        )

    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f'{level}_weighted_f1_heatmap_threshold_{threshold:.2f}.png'))
    plt.close()

def plot_f1_heatmaps_by_level(weighted_f1_data, threshold, outpath, ref_keys, acronym_mapping=None):
    """Main function to plot F1 heatmaps by level."""
    sns.set(style="whitegrid")
    os.makedirs(outpath, exist_ok=True)

    vmin, vmax = weighted_f1_data['weighted_f1'].min(), weighted_f1_data['weighted_f1'].max()

    for level in ref_keys:
        level_data = weighted_f1_data[weighted_f1_data['key'] == level]
        plot_f1_heatmap_for_level(level_data, vmin, vmax, level, threshold, outpath, acronym_mapping)


    
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
        temp_df = pd.read_csv(file,sep="\t")
        f1_df = pd.concat([temp_df, f1_df], ignore_index=True)
    #need to flatten all_f1_scores into data frame
    f1_df["reference_acronym"] = f1_df["reference"].apply(make_acronym)
    f1_df["query_acronym"] = f1_df["query"].apply(make_acronym)
    f1_df["reference"] = f1_df["reference"].str.replace("_", " ")
 
    acronym_mapping_ref = f1_df[["reference", "reference_acronym"]].drop_duplicates().set_index("reference")["reference_acronym"].to_dict()
    acronym_mapping_query = f1_df[["query", "query_acronym"]].drop_duplicates().set_index("query")["query_acronym"].to_dict()
    
    #plot_distribution(f1_df,var="f1_score",outdir="dists",split="label",facet="reference")
    plot_distribution(f1_df, var="weighted_f1",outdir="dists", split="reference_acronym", facet="key", 
                      acronym_mapping = acronym_mapping_ref)
    plot_distribution(f1_df, var="weighted_f1",outdir="dists", split="query_acronym", facet="key", 
                      acronym_mapping = acronym_mapping_query)

    for file in f1_results:
        # Get the full path of the file
        #file_path = os.path.join(f1_results_dir, file)
        
        # Read each file into a DataFrame
        x = pd.read_csv(file, sep="\t")
        
        # Get the unique keys from the 'key' column
        keys = x["key"].unique()
        
        # Loop through each unique key
        for key in keys:
            # Filter rows where the 'key' matches
            subset = x[x["key"] == key]
            subset["query"] = subset["query"].str.replace("_processed_", "")
            subset["query_acronym"] = subset["query"].apply(make_acronym)
            subset["reference_acronym"] = subset["reference"].apply(make_acronym)
            subset["reference"] = subset["reference"].str.replace("_", " ")
            # If this key doesn't exist in the dictionary, initialize an empty DataFrame
            if key not in all_f1_scores:
                all_f1_scores[key] = subset
            else:
                # If the key already exists, concatenate the new subset to the existing DataFrame
                all_f1_scores[key] = pd.concat([all_f1_scores[key], subset], ignore_index=True)
    
    plot_label_f1_heatmaps(all_f1_scores, threshold=cutoff, outpath="f1_plots", widths=[1,0.8,0.5], 
                           acronym_mapping=None)
    
    final_f1_data = pd.DataFrame()
    for key, df in all_f1_scores.items():
        macro = df.drop(columns=['label', 'f1_score'])
        macro["key"] = key
        final_f1_data = pd.concat([final_f1_data, macro], ignore_index=True)
    weighted_f1_data = final_f1_data[['reference', 'reference_acronym','key', 'query', 'weighted_f1']]
 
    plot_f1_heatmaps_by_level(weighted_f1_data, threshold=cutoff, outpath="f1_plots", 
                              ref_keys=ref_keys, acronym_mapping=None)
    


if __name__ == "__main__":
    main()
    
