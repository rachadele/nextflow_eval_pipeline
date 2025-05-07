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
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mmus/03/3540ff0c3d691ffe226be26dc011c1/GSE247339.2_1051967_processed.h5ad")
    parser.add_argument('--predicted_meta', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mus_musculus/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi/predicted_meta/GSE247339.2_1051967/whole_cortex/GSE247339.2_1051967_whole_cortex.predictions.0.0.tsv")
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


def is_outlier(query, metric: str, nmads=3):
    M = query.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier



def mad(var, scale='normal'):
    """Median Absolute Deviation. Set scale='normal' for consistency with R's default."""
    med = np.median(var)
    mad = np.median(np.abs(var - med))
    if scale == 'normal':
        return mad * 1.4826  # for normally distributed data
    return mad


def get_lm(query, nmads=5, scale="normal"):
    # Assume dataset is an AnnData object
    # Fit linear model: log10(n_genes_per_cell) ~ log10(counts_per_cell)
    lm_model = ols(formula='log1p_n_genes_by_counts ~ log1p_total_counts', data=query.obs).fit()
    # Calculate residuals
    residuals = lm_model.resid
    # If data is normally distributed, this is similar to std 
    mad_residuals = mad(residuals, scale=scale)
    # Intercept adjustment (upper bound)
    intercept_adjustment = np.median(residuals) + nmads * mad_residuals
    return {
        "model": lm_model,
        "intercept_adjustment": intercept_adjustment
    }

def is_correct(adata, ref_keys, mapping_df):
    adata.obs = map_valid_labels(adata.obs, ref_keys = ref_keys, mapping_df = mapping_df) 
    # change to string type
    adata.obs["correct"] = adata.obs["predicted_subclass"].astype(str) == adata.obs["subclass"].astype(str)
    return adata
    
    
def qc_preprocess(query):
    sc.pp.scrublet(query, batch_key="sample_id")
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query)
    
    return query

def get_qc_metrics(query, nmads):
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20], use_raw=True)

    metrics = {
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb",
    }
    
    for metric, col_name in metrics.items():
        query.obs[col_name] = is_outlier(query, metric, nmads)

    lm_dict = get_lm(query, nmads=nmads)
    intercept = lm_dict["model"].params[0]
    slope = lm_dict["model"].params[1]
    

    query.obs["counts_outlier"] = (
        query.obs["log1p_n_genes_by_counts"] < (query.obs["log1p_total_counts"] * slope + (intercept - lm_dict["intercept_adjustment"]))
        ) | (
        query.obs["log1p_n_genes_by_counts"] > (query.obs["log1p_total_counts"] * slope + (intercept + lm_dict["intercept_adjustment"]))
        )

    query.obs["total_outlier"] = (
        query.obs["counts_outlier"] | query.obs["outlier_mito"] | query.obs["outlier_ribo"] | query.obs["outlier_hb"] | query.obs["predicted_doublet"]
    )
    
    query.obs["non_outlier"] = ~query.obs["total_outlier"]

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
        y="pct_counts_mito",
        hue="outlier_mito",
        kind="scatter"
    )
    plot2.savefig(f"{query_name}_mito_counts.png")

    # Third jointplot
    plot3 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="pct_counts_ribo",
        hue="outlier_ribo",
        kind="scatter"
    )
    plot3.savefig(f"{query_name}_ribo_counts.png")
    
    plot4 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="pct_counts_hb",
        hue="outlier_hb",
        kind="scatter"
    )
    plot4.savefig(f"{query_name}_hb_counts.png")

    

def plot_jointplots_R(query, study_name, sample_name):
    os.makedirs(study_name, exist_ok=True)
    # Save query.obs to CSV
    query.obs.to_csv(f"{study_name}/{sample_name}_obs.tsv", sep="\t", index=False)
    tsv_path = os.path.abspath(f"{study_name}/{sample_name}_obs.tsv")

    # get path to Rscript
    rscript_path = os.path.join(os.path.dirname(__file__), "plot_jointplots.R")
    subprocess.run([
        "Rscript", rscript_path,
        tsv_path, study_name, sample_name
    ])

        

def plot_umap_qc(query, study_name, sample_name):
    colors = ["outlier_hb", "outlier_ribo", "outlier_mito","predicted_doublet","counts_outlier","total_outlier"]

    output_dir = os.path.join(study_name, sample_name)
    os.makedirs(output_dir, exist_ok=True)

    sc.pl.umap(
        query,
        color=colors,
        use_raw=False,
        save=None,
        show=False,
       # title=f"Sample {sample_name}",
        ncols=2)
        # Manually save the plot
    plt.savefig(os.path.join(output_dir, "umap_mqc.png"), dpi=150, bbox_inches='tight')
    plt.close()
            

def read_markers(markers_file, organism):
    df = pd.read_csv(markers_file, sep=None, header=0)
    
    # Split markers column into list
    df['markers'] = df['markers'].str.split(',\s*', regex=True)

    # Build nested dict: family > class > cell_type
    nested_dict = defaultdict(lambda: defaultdict(dict))
    
    for _, row in df.iterrows():
        fam = row['family']
        cls = row['class']
        cell = row['cell_type']
        markers = row['markers']
        nested_dict[fam][cls][cell] = markers
    
    if organism == "mus_musculus":
        for fam in nested_dict:
            for cls in nested_dict[fam]:
                for cell in nested_dict[fam][cls]:
                    nested_dict[fam][cls][cell] = [
                        x.lower().capitalize() for x in nested_dict[fam][cls][cell]
                    ]

    return nested_dict


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

def plot_ct_composition(query, subclass_colors, query_name):

    subclass_counts = query.obs["subclass"].value_counts()
    predicted_subclass_counts = query.obs["predicted_subclass"].value_counts()
    # Union of all subclasses
    all_subclasses = subclass_counts.index.union(predicted_subclass_counts.index)
    subclass_counts = subclass_counts.reindex(all_subclasses, fill_value=0)
    predicted_subclass_counts = predicted_subclass_counts.reindex(all_subclasses, fill_value=0)

   # Combine into DataFrame
    df_counts = pd.DataFrame({
        "True": subclass_counts,
        "Predicted": predicted_subclass_counts
    })

    # Crosstab of predicted vs true
    crosstab = pd.crosstab(query.obs["predicted_subclass"], query.obs["subclass"])
    crosstab = crosstab.reindex(index=all_subclasses, columns=all_subclasses, fill_value=0)

    df_counts = df_counts.T
    # First figure: stacked barh of total counts
    fig1, ax1 = plt.subplots(figsize=(10, 15))
    df_counts.plot(kind="barh", stacked=True, ax=ax1, color=subclass_colors)
    ax1.set_xlabel("Cell Count")
    ax1.set_ylabel("Subclass")
    ax1.set_title("True vs Predicted Subclass Composition")
    ax1.legend(
        title="Label Type",
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0.
    )
    plt.tight_layout()
    plt.savefig(f"{query_name}_subclass_totals.png", dpi=300, bbox_inches="tight")
    plt.close(fig1)

    # Second figure: breakdown of true subclass by predicted subclass
    fig2, ax2 = plt.subplots(figsize=(10, 15))
    crosstab.plot(kind="barh", stacked=True, ax=ax2, color=subclass_colors)
    ax2.set_xlabel("Cell Count")
    ax2.set_ylabel("Predicted Subclass")
    ax2.set_title("Predicted Subclass Composition by True Subclass")
    plt.tight_layout()
    plt.savefig(f"{query_name}_subclass_breakdown.png", dpi=300, bbox_inches="tight")
    plt.close(fig2) 

    
    
def plot_misclassified(query, subclass_colors, query_name):
    query_obs = query.obs.copy()

    # Add classification result to obs
    query_obs["classification"] = query_obs["correct"].map({True: "Correct", False: "Incorrect"})

    # Plot distribution of each QC metric grouped by subclass and classification
    metrics = ["pct_counts_mito", "pct_counts_ribo", "pct_counts_hb", "predicted_doublet"]
    for metric in metrics:
        plt.figure(figsize=(10, 5))
        sns.boxplot(data=query_obs, x="subclass", y=metric, hue="classification", palette="Set2")
        plt.title(f"{query_name}: {metric} by subclass (Correct vs Incorrect)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(f"{query_name}_{metric}_misclassified.png", dpi=300)
        plt.close()

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
    query = qc_preprocess(query)
    query = get_qc_metrics(query, query_name)
    # for sample in query.obs["sample_id"].unique():
    #    query = query[query.obs["sample_id"] == sample]
    query = is_correct(query, ref_keys, mapping_df)
   # plot_jointplots(query, query_name)
    plot_ct_composition(query, subclass_colors, query_name)
    plot_misclassified(query, subclass_colors, query_name)
   # plot_umap_qc(query, query_name, subclass_colors)
    

if __name__ == "__main__":
    main() 