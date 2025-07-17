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
#get original files
# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--study_name', type=str, default="PTSDBrainomics", help='Name of the study')
    parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--predicted_meta_files', type=str, nargs="+", default = ["/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/homo_sapiens/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi/PTSDBrainomics_1203506/Dissection_Angular_gyrus_AnG/predicted_meta/PTSDBrainomics_1203506_Dissection_Angular_gyrus_AnG.predictions.0.0.tsv",
                                                                         "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/homo_sapiens/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi/PTSDBrainomics_1203507/Dissection_Angular_gyrus_AnG/predicted_meta/PTSDBrainomics_1203507_Dissection_Angular_gyrus_AnG.predictions.0.0.tsv",
                                                                         "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/homo_sapiens/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi/PTSDBrainomics_1203508/Dissection_Angular_gyrus_AnG/predicted_meta/PTSDBrainomics_1203508_Dissection_Angular_gyrus_AnG.predictions.0.0.tsv"])
    parser.add_argument('--query_paths', type=str, nargs="+", default=["/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/sample_subsets/PTSDBrainomics_1203506.h5ad",
                                                                                "/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/sample_subsets/PTSDBrainomics_1203507.h5ad",
                                                                                "/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/sample_subsets/PTSDBrainomics_1203508.h5ad"])
    parser.add_argument('--markers_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/cell_type_markers.tsv")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
    parser.add_argument('--nmads',type=int, default=5)
    parser.add_argument('--ref_keys', nargs="+", type=str, default = ["subclass","class","family"])
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    #parser.add_argument('--sample_key', type=str, default="sample_id")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    

def is_correct(adata, level="subclass"):
  # change to string type
  adata.obs["correct_"+level] = adata.obs["predicted_"+level].astype(str) == adata.obs[level].astype(str)
  return adata

def read_query(query_path, gene_mapping, predicted_meta):

    query = sc.read_h5ad(query_path)
    # check if query is empty
    if query.n_obs == 0:
        raise ValueError(f"Query dataset {query_path} is empty after reading.")
    if "sample_id" in query.obs.columns and "cell_id" in query.obs.columns:
        query.obs["full_barcode"] = query.obs["sample_id"].astype(str) + "_" + query.obs["cell_id"].astype(str)
        predicted_meta["full_barcode"] = predicted_meta["sample_id"].astype(str) + "_" + predicted_meta["cell_id"].astype(str)
        query.obs = query.obs.merge(predicted_meta, left_on="full_barcode", right_on="full_barcode", how="left", suffixes=("", "_y"))
    
    else:
        query.obs.index = query.obs.index.astype(int)
        query.obs = query.obs.merge(predicted_meta, left_index=True, right_index=True, how="left", suffixes=("", "_y"))

        
    columns_to_drop = [col for col in query.obs.columns if col.endswith("_y")]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query


def qc_preprocess(query):
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=0.3)
    
    return query



def plot_joint_umap(query, outdir):
    x_metric = "log1p_n_genes_by_counts"
    metrics = {
        "log1p_total_counts": ["counts_outlier", "umi_outlier", "genes_outlier"],
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb",
    }
    
    data = query.obs
    images = []
    for yval, hues in metrics.items():
        for hue in (hues if isinstance(hues, list) else [hues]):
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
    
    all_subclasses = sorted(set(query.obs["subclass"].unique()) | 
                            set(query.obs["predicted_subclass"].unique()))

    # Generate unique colors for each subclass
    color_palette = sns.color_palette("husl", n_colors=len(all_subclasses))
    subclass_colors = dict(zip(all_subclasses, color_palette))
    os.makedirs(outdir, exist_ok=True)

    for col in ["predicted_subclass", "subclass", "correct_subclass", "sample_id", "predicted_doublet"]:
        sc.pl.umap(
            query,
            color=col,
            use_raw=False,
            palette=subclass_colors if col in ["subclass", "predicted_subclass"] else None,
            save=None,
            show=False,
        )

        out_path = os.path.join(outdir, f"{col}_umap_mqc.png") 
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()



def plot_upset_by_group(obs, outlier_cols, group_col, outdir):
    os.makedirs(outdir, exist_ok=True)
    obs = obs.copy()
    obs["membership"] = obs[outlier_cols].apply(lambda row: tuple(c for c in outlier_cols if row[c]), axis=1)

    if group_col:
        for group in sorted(obs[group_col].unique()):
            counts = obs[obs[group_col] == group]["membership"].value_counts()
            if counts.empty:
                continue
    else:   
        group_col = "study"
        group = "all"
        counts = obs["membership"].value_counts()
        if counts.empty:
            return
    data = from_memberships(counts.index, data=counts.values)
    plt.figure(figsize=(8, 4))
    UpSet(data, show_counts=True).plot()
    plt.suptitle(f"{group_col} = {group}")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{group_col}_{group}_upset_mqc.png"))
    plt.close()


def main():
    # Parse command line arguments
    args = parse_arguments()
    # Set variables from arguments
    query_paths = args.query_paths
    predicted_meta_files = args.predicted_meta_files
    # print to std err
    if not query_paths or not predicted_meta_files:
        raise ValueError("Please provide at least one query path and one predicted meta file.")
    if len(query_paths) != len(predicted_meta_files):
        raise ValueError("Number of query paths must match number of predicted meta files")
    
    markers_file = args.markers_file
    gene_mapping_path = args.gene_mapping 
    organism = args.organism
    mapping_file = args.mapping_file
    ref_keys = args.ref_keys
    study_name = args.study_name
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
    os.makedirs(study_name, exist_ok=True) 
    # join query paths and predicted meta files by sample_id (split by "_")
    if len(query_paths) != len(predicted_meta_files):
        raise ValueError("Number of query paths must match number of predicted meta files")
   
    predicted_meta = {}
    for file in predicted_meta_files:
        sample_id = os.path.basename(file).split("_")[1]
        predicted_meta[sample_id] = pd.read_csv(file, sep=None, header=0)
        
    all_query_samples = {}
    for query_path in query_paths:
        # Extract sample name from query path
        sample_id = os.path.basename(query_path).split("_")[1]
       # sample_id = os.path.basename(query_path).split("_")[1].split(".")[0]

        # print to std err
        if not sample_id:
            raise ValueError(f"Sample ID could not be extracted from query path: {query_path}")
        assigned_celltypes = predicted_meta.get(sample_id, []) 
        if not predicted_meta:
            raise ValueError(f"No predicted meta file found for sample {sample_id}")
        query = read_query(query_path, gene_mapping, predicted_meta=assigned_celltypes)
        for key in ref_keys:
            query = is_correct(query, level=key)
        all_query_samples[str(sample_id)] = query


    # if length of all_query_samples is 0, raise error
    if len(all_query_samples.values()) == 0:
        raise ValueError("No valid query samples found. Please check the query paths and predicted meta files.")
    # merge
    query = ad.concat(all_query_samples, label="sample_id",axis=0)
    if "feature_name" not in query.var.columns:
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)

    # print to std err
    #query.obs.index = query.obs["cell_id"]
    query.raw = query.copy()
    

    query=qc_preprocess(query)
    
    query.obs["non_outlier"] = ~(
        query.obs["counts_outlier"] |
        query.obs["outlier_mito"] |
        query.obs["outlier_ribo"] |
        query.obs["outlier_hb"] |
        query.obs["predicted_doublet"]
    )
#plot_markers(query, markers_file, organism=organism)
    make_celltype_matrices(query, markers_file, organism=organism, outdir=study_name)
    
    
    for sample_id in query.obs["sample_id"].unique():
        sample_query = query[query.obs["sample_id"] == sample_id]
        
       # sample_query = get_qc_metrics(sample_query, nmads=args.nmads)
        outdir = os.path.join(study_name, sample_id)
        os.makedirs(outdir, exist_ok=True)
        plot_joint_umap(sample_query, outdir=outdir)
   
    
    plot_ct_umap(query, outdir=study_name)
    # Count occurrences
    celltype_counts_correct = (
        query.obs
        .groupby(["subclass", "correct_subclass"])
        .size()                             # count cells per (sample, subclass)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    )
    celltype_counts_correct.to_csv(os.path.join(study_name,"celltype_counts_correct_mqc.tsv"), sep="\t", index=False)

    predicted_vs_actual_counts = (
        query.obs
        .groupby(["subclass", "predicted_subclass"])
        .size()                             # count cells per (sample, subclass)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    ) 
    predicted_vs_actual_counts.to_csv(os.path.join(study_name,"actual_vs_predicted_counts_mqc.tsv"), sep="\t", index=False)
    
    
    sample_correct_counts = (
        query.obs.groupby(["sample_id","correct_subclass"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    sample_correct_counts.to_csv(os.path.join(study_name,"sample_correct_counts_mqc.tsv"), sep="\t", index=False)
    
    
   outlier_cols = [
        "counts_outlier", 
        "umi_outlier", 
        "genes_outlier",
        "outlier_mito", 
        "outlier_ribo", 
        "outlier_hb", 
        "predicted_doublet"
    ]
    plot_upset_by_group(query_combined.obs, outlier_cols, None, study_name)
 
if __name__ == "__main__":
    main()
 
    
    

   