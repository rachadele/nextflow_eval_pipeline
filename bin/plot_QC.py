#!/user/bin/python3
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
import warnings
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
import adata_functions
from adata_functions import *
from PIL import Image
import io
import os
import math
from scipy.stats import median_abs_deviation
from statsmodels.formula.api import ols


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/hsap/fc/928406f842fb7357c49eefb0328834/rosmap_R9426782_processed.h5ad")
    parser.add_argument('--predicted_meta', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/hsap/fc/928406f842fb7357c49eefb0328834/rosmap_R9426782_Dissection_Angular_gyrus_AnG.predictions.0.0.tsv")
    parser.add_argument('--markers_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/cell_type_markers.tsv")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
    parser.add_argument('--nmads',type=int, default=5)
    parser.add_argument('--ref_keys', nargs="+", type=str, default = ["subclass","class","family"])
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    #parser.add_argument('--sample_key', type=str, default="sample_id")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    
def is_correct(adata, ref_keys, mapping_df):
    adata.obs = map_valid_labels(adata.obs, ref_keys = ref_keys, mapping_df = mapping_df) 
    # change to string type
    adata.obs["correct"] = adata.obs["predicted_subclass"].astype(str) == adata.obs["subclass"].astype(str)
    return adata

def read_query(query_path, gene_mapping, predicted_meta):
    
    query = sc.read_h5ad(query_path)
    if "feature_name" not in query.var.columns:
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)

    if "sample_id" in query.obs.columns and "cell_id" in query.obs.columns:
        query.obs["full_barcode"] = query.obs["sample_id"].astype(str) + "_" + query.obs["cell_id"].astype(str)
        predicted_meta["full_barcode"] = predicted_meta["sample_id"].astype(str) + "_" + predicted_meta["cell_id"].astype(str)
        query.obs = query.obs.merge(predicted_meta, left_on="full_barcode", right_on="full_barcode", how="left", suffixes=("", "_y"))
    
    else:
        # they should be in the same order
        query.obs = query.obs.merge(predicted_meta, left_index=True, right_index=True, how="left", suffixes=("", "_y"))
        
    columns_to_drop = [col for col in query.obs.columns if col.endswith("_y")]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query


def is_outlier(query, metric: str, nmads=3):
    M = query.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def qc_preprocess(query):
    # check if any sample_id has fewer than 30 associated cells
    #sample_counts = query.obs["sample_id"].value_counts()
    #if (sample_counts < 30).any():
        #batch_key=None
    #else:
        #batch_key="sample_id"
   # sc.pp.scrublet(query, batch_key=batch_key)
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=0.3)
    
    return query

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
    mad_residuals = median_abs_deviation(residuals, scale=scale)
    # Intercept adjustment (add for upper bound, subtract for lower bound)
    intercept_adjustment = np.median(residuals) + nmads * mad_residuals
    return {
        "model": lm_model,
        "intercept_adjustment": intercept_adjustment
    }
    
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
        "log1p_total_counts": "umi_outlier",
        "log1p_n_genes_by_counts": "genes_outlier",
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
        ) | (
        query.obs["umi_outlier"] ) | (query.obs["genes_outlier"])
        

    query.obs["total_outlier"] = (
        query.obs["counts_outlier"] | query.obs["outlier_mito"] | query.obs["outlier_ribo"] | query.obs["outlier_hb"]
    )
    
    query.obs["non_outlier"] = ~query.obs["total_outlier"]

    return query

def map_celltype_hierarchy(query, markers_file):
    # Load the markers table
    df = pd.read_csv(markers_file, sep=None, header=0)
    df.drop(columns="markers", inplace=True)
    query.obs = query.obs.merge(df, left_on="subclass", right_on="subclass", how="left", suffixes=("", "_y"))
    return query

def get_gene_to_celltype_map(markers_file, organism="mus_musculus"):
    # Read the marker file
    df = pd.read_csv(markers_file, sep="\t")
    df = df[df["markers"].notnull()]
    gene_to_celltype = {}

    for _, row in df.iterrows():
        subclass = row["cell_type"]
        genes = [gene.strip() for gene in row["markers"].split(",")]
        for gene in genes:
            if organism == "mus_musculus":
                gene = gene.lower().capitalize()
                # Handle multiple cell types mapping to the same gene
            if gene not in gene_to_celltype:
                gene_to_celltype[gene] = []
            gene_to_celltype[gene].append(subclass)

    # Join multiple cell types into one label if needed
    gene_ct_dict = {
        gene: f"{gene}: {'_'.join(set(celltypes))}"
        for gene, celltypes in gene_to_celltype.items()
    }
    return gene_ct_dict

def make_celltype_matrices(query, markers_file, organism="mus_musculus", outdir=""):
    # Drop vars with NaN feature names
    # Map cell type hierarchy
   # query = map_celltype_hierarchy(query, markers_file=markers_file)
    
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]
    
    #Make raw index match processed var index
    query.raw.var.index = query.raw.var["feature_name"]

    # Read marker genes
    gene_ct_dict = get_gene_to_celltype_map(markers_file, organism=organism)
    # Collect all unique markers across all families/classes/cell types
    all_markers = list(gene_ct_dict.keys())
    valid_markers = [gene for gene in all_markers if gene in query.var_names]

    # Filter raw expression matrix to match query.var_names
    expr_matrix = query.raw.X.toarray()
    expr_matrix = pd.DataFrame(expr_matrix, index=query.obs.index, columns=query.raw.var.index)
    
    avg_expr = expr_matrix.groupby(query.obs["subclass"]).mean()
    avg_expr = avg_expr.loc[:, valid_markers]
    
    # Scale expression across genes
    scaled_expr = (avg_expr - avg_expr.mean()) / avg_expr.std()
    scaled_exp = scaled_expr.loc[:, valid_markers]
    scaled_expr.fillna(0, inplace=True)

    # Rename columns: gene -> gene (celltype)
    scaled_expr.rename(columns=gene_ct_dict, inplace=True)

    # Save matrix
    scaled_expr.to_csv(f"{outdir}/heatmap_mqc.tsv", sep="\t")

 

def plot_joint_umap(query, outdir):
    x_metric = "log1p_n_genes_by_counts"
    metrics = {
        "log1p_total_counts": "counts_outlier",
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb",
    }
    
    data = query.obs
    images = []
    for yval, hue in metrics.items():
        fig_joint = sns.jointplot(
            data=data,
            x=x_metric,
            y=yval,
            hue=hue,
            kind="scatter"
        )
        
        
        umap_fig = sc.pl.umap(
        query,
        color=hue,
        use_raw=False,
        save=None,
        show=False,
        title=f"{hue}",
        ncols=1,
        legend_loc="upper right",
        return_fig=True
        ) 


        joint_buf = io.BytesIO()
        fig_joint.savefig(joint_buf, format="png", bbox_inches='tight')
        plt.close(fig_joint.fig) 
        
        umap_buf = io.BytesIO()
        umap_fig.savefig(umap_buf, format="png", bbox_inches='tight')
        plt.close(umap_fig)
        
        joint_buf.seek(0)
        images.append(Image.open(joint_buf))
        
        umap_buf.seek(0)
        images.append(Image.open(umap_buf))
    
    scale = 0.5  # Resize to 50%
    resized_images = [img.resize((int(img.width * scale), int(img.height * scale))) for img in images]

    # Use resized dimensions
    img_width, img_height = resized_images[0].size
    grid_cols = 2
    grid_rows = math.ceil(len(resized_images) / grid_cols)

    combined_img = Image.new("RGB", (grid_cols * img_width, grid_rows * img_height), "white")

    for idx, img in enumerate(resized_images):
        row = idx // grid_cols
        col = idx % grid_cols
        x_offset = col * img_width
        y_offset = row * img_height
        combined_img.paste(img, (x_offset, y_offset))

    os.makedirs(outdir, exist_ok=True)
    out_path = f"{outdir}/outliers_mqc.png"
    combined_img.save(out_path)

def plot_ct_umap(query, outdir):
    colors = ["predicted_subclass", "subclass", "correct"]
    
    sc.pl.umap(
        query,
        color=colors,
        use_raw=False,
        save=None,
        show=False,
        ncols=1
       # legend_loc="upper right",
    )
    os.makedirs(outdir, exist_ok=True)
    out_path = f"{outdir}/celltype_umap_mqc.png"
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()

def main():
    # Parse command line arguments
    args = parse_arguments()
    # Set variables from arguments
    query_path = args.query_path
    predicted_meta = args.predicted_meta
    markers_file = args.markers_file
    gene_mapping_path = args.gene_mapping 
    organism = args.organism
    mapping_file = args.mapping_file
    ref_keys = args.ref_keys
    #sample_key = args.sample_key
    # Load the mapping file
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=0)
    
    gene_mapping = pd.read_csv(gene_mapping_path, sep=None, header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True) 

    # Load query and reference datasets
    study_name = os.path.basename(query_path).replace("_processed.h5ad", "")
    
    os.makedirs(study_name, exist_ok=True)
    outdir = study_name
    
    assigned_celltypes = pd.read_csv(predicted_meta, sep=None, header=0)
     
    query = read_query(query_path, gene_mapping, predicted_meta=assigned_celltypes)

    query.obs.index = query.obs["index"]
    query.raw = query.copy()
    query = qc_preprocess(query)
    query = is_correct(query, ref_keys, mapping_df)
    #plot_markers(query, markers_file, organism=organism)
    make_celltype_matrices(query, markers_file, organism=organism, outdir=outdir)
    query = get_qc_metrics(query, nmads=args.nmads) 
    
    plot_joint_umap(query, outdir=outdir)
    plot_ct_umap(query, outdir=outdir)
    # Count occurrences
    celltype_counts_correct = (
        query.obs
        .groupby(["subclass", "correct"])
        .size()                             # count cells per (sample, subclass)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    )
    celltype_counts_correct.to_csv(os.path.join(outdir,"celltype_counts_correct_mqc.tsv"), sep="\t", index=False)

    ## make a table of counts by outliers
    # count all combinations + non-outliers
    celltype_outlier_counts = (
        query.obs
        .groupby("subclass")[["counts_outlier", "outlier_mito", "outlier_ribo", "outlier_hb", "non_outlier"]]
        .sum()
        .astype(int)
    )
    celltype_outlier_counts.to_csv(os.path.join(outdir, "celltype_outlier_counts_mqc.tsv"), sep="\t", index=True)
 
    # correct by outlier composition
    correct_outlier_counts = (
        query.obs
        .groupby(["correct"])[["counts_outlier", "outlier_mito", "outlier_ribo", "outlier_hb", "non_outlier"]]
        .sum()
        .astype(int)
    )
    correct_outlier_counts.to_csv(os.path.join(outdir,"correct_outlier_counts_mqc.tsv"), sep="\t", index=True)
    
    predicted_vs_actual_counts = (
        query.obs
        .groupby(["subclass", "predicted_subclass"])
        .size()                             # count cells per (sample, subclass)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    ) 
    predicted_vs_actual_counts.to_csv(os.path.join(outdir,"actual_vs_predicted_counts_mqc.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
 
    
    

   