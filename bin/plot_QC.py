#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
from scipy.sparse import csr_matrix
import warnings

import adata_functions
from adata_functions import *

import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
from scipy.stats import median_abs_deviation
from adata_functions import *

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
   # parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mmus/7f/b447fabc3ba6bd819a1e103ef04f26/GSE124952_483958_processed.h5ad")
    parser.add_argument('--predicted_meta', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mus_musculus/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi/predicted_meta/GSE124952_483958_whole_cortex.predictions.0.0.tsv")
    parser.add_argument('--markers_file', type=str, default=""),
    parser.add_argument('--gene_mapping', default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv", type=str),
    parser.add_argument('--ref_name', type=str),
    parser.add_argument('--ref_keys', type = str, nargs="+", default = ["subclass","class","family","global"]),
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
    parser.add_argument('--color_mapping', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/color_mapping.tsv")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


    
def read_adata(query_path, gene_mapping, predicted_meta):
    
    query = sc.read_h5ad(query_path)
    if "feature_name" not in query.var.columns:
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)

    query.obs["full_barcode"] = query.obs["sample_id"].astype(str) + "_" + query.obs["cell_id"].astype(str)
    predicted_meta["full_barcode"] = predicted_meta["sample_id"].astype(str) + "_" + predicted_meta["cell_id"].astype(str)
    
    query.obs = query.obs.merge(predicted_meta, left_on="full_barcode", right_on="full_barcode", how="left", suffixes=("", "_y"))
    columns_to_drop = [col for col in query.obs.columns if col.endswith("_y")]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query


def is_correct(adata, ref_keys, mapping_df):
    adata.obs = map_valid_labels(adata.obs, ref_keys = ref_keys, mapping_df = mapping_df) 
    # change to string type
    adata.obs["correct"] = adata.obs["predicted_subclass"].astype(str) == adata.obs["subclass"].astype(str)
    return adata
    
    
def process_adata(query):
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query)
    
    return query


def plot_jointplots(query, query_name):

    # First jointplot
    plot1 = sns.jointplot(
        data=query.obs,
        x="log1p_total_counts",
        y="log1p_n_genes_by_counts",
        hue="counts_outlier",
        kind="scatter"
    )
    plot1.savefig(f"{query_name}_genes_vs_counts.png")

    # Second jointplot
    plot2 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts_mito",
        hue="outlier_mito",
        kind="scatter"
    )
    plot2.savefig(f"{query_name}_mito_counts.png")

    # Third jointplot
    plot3 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts_ribo",
        hue="outlier_ribo",
        kind="scatter"
    )
    plot3.savefig(f"{query_name}_ribo_counts.png")
    
    plot4 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts_hb",
        hue="outlier_hb",
        kind="scatter"
    )
    plot4.savefig(f"{query_name}_hb_counts.png")

    
def plot_umap_qc(query, query_name, subclass_colors):
    colors = ["subclass", "predicted_subclass", "correct", "pct_counts_mito", "pct_counts_ribo", "pct_counts_hb","predicted_doublet"]

    for color in colors:
        is_categorical = query.obs[color].dtype.name == "category" or query.obs[color].dtype == object

        sc.pl.umap(
            query,
            color=color,
            use_raw=False,
            save=f"_{query_name}_{color}_qc_umap.png",
            show=False,
            title=f"{color} - {query_name}",
            palette=subclass_colors if is_categorical else None,
            color_map="viridis" if not is_categorical else None,
        )

def make_stable_colors(color_mapping_df):
    
    all_subclasses = sorted(color_mapping_df["subclass"])
    # i need to hardcode a separate color palette based on the mmus mapping file
    # Generate unique colors for each subclass
    color_palette = sns.color_palette("husl", n_colors=len(all_subclasses))
    subclass_colors = dict(zip(all_subclasses, color_palette))
    return subclass_colors
    
    
def plot_marker_genes(query, marker_file):
    if not os.path.exists(marker_file):
        raise ValueError(f"Marker file {marker_file} does not exist.")
    if "feature_name" not in query.var.columns:
        raise ValueError("feature_name column not found in query.var")
    query.var.set_index("feature_name", inplace=True)
    sc.pl.umap(
        query,
        color = ["Olig1", "Cldn5", "Aif1"],
        use_raw=False,
        #save=f"qc_metrics_{query_name}_{sample_id}.png",
        show=True,
       # title=f"QC Metrics for {query_name} {sample_id}",
        ncols=1,
        color_map="viridis")
    
        
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()
    color_mapping = args.color_mapping
    ref_keys = args.ref_keys
    mapping_file = args.mapping_file
    # Load the mapping file
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=0)
    color_mapping = pd.read_csv(color_mapping, sep="\t", header=0)
    subclass_colors = make_stable_colors(color_mapping)
     
    # Set variables from arguments
    query_path = args.query_path
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    predicted_meta = args.predicted_meta
    #markers_file = args.markers_file
    gene_mapping_path = args.gene_mapping
    gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)
    # Load query and reference datasets
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    predicted_meta = pd.read_csv(predicted_meta, sep="\t", header=0)
   # markers = pd.read_csv(markers_file, sep="\t", header=0)
    
    query = read_adata(query_path, gene_mapping, predicted_meta)
    query = process_adata(query)
    #query = get_qc_metrics(query, query_name)
    # for sample in query.obs["sample_id"].unique():
    #    query = query[query.obs["sample_id"] == sample]
    query = is_correct(query, ref_keys, mapping_df)
    plot_jointplots(query, query_name)
    plot_umap_qc(query, query_name, subclass_colors)
    

if __name__ == "__main__":
    main() 