#!/usr/bin/python3

from pathlib import Path
import random
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
import cellxgene_census.experimental
from sklearn.ensemble import RandomForestClassifier
import adata_functions
from adata_functions import *
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import yaml
import adata_functions
from adata_functions import *


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Download model file based on organism, census version, and tree file."
    )
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--ref_collections', type=str, nargs='+', default=[
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Adult mouse cortical cell taxonomy revealed by single cell transcriptomics"
    ])
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--outdir', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations")
    parser.add_argument('--organ', type=str, default="brain")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def extract_obs(filtered_ids, organism=None, census=None, 
    obs_filter=None, cell_columns=None, dataset_info=None, seed=42):
    # Assuming get_seurat is defined to return an AnnData object
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_value_filter=obs_filter,  # Ensure this is constructed correctly
        obs_column_names=cell_columns,
        obs_coords=filtered_ids,
        var_value_filter = "nnz > 10")
   #sc.pp.filter_genes(adata, min_cells=3)
    
    newmeta = adata.obs.merge(dataset_info, on="dataset_id", suffixes=(None,"y"))
    adata.obs = newmeta
    # Assuming relabel_wrapper is defined
    #adata = relabel(adata, relabel_path=relabel_path, join_key="cell_type", sep='\t')
    # Convert all columns in adata.obs to factors (categorical type in pandas)
    return adata.obs

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set variables from arguments
    organism = args.organism
    census_version = args.census_version
    # need to not subsample ref
    ref_collections = args.ref_collections
    SEED = args.seed
    outdir = args.outdir
    organ=args.organ 
    os.makedirs(outdir, exist_ok=True)

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    brain_obs = get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
    
    brain_obs = brain_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    brain_obs.drop(columns=['soma_joinid_y'], inplace=True)
    brain_obs_filtered = brain_obs[brain_obs['collection_name'].isin(ref_collections)] 

    organism_name_mapping = {
        "homo_sapiens": "Homo sapiens",
        "mus_musculus": "Mus musculus"
    }
    organism = organism_name_mapping.get(organism, organism)

    cell_columns = [
        "dataset_id",
        "soma_joinid"
    ]
   
    filtered_ids=brain_obs_filtered["soma_joinid"].values
    all_collection_obs = extract_obs(filtered_ids, organism=organism, census=census, cell_columns=cell_columns, 
                                     dataset_info=dataset_info, seed=42)
    

    for dataset_id in all_collection_obs["dataset_id"].unique():
        dataset_title = all_collection_obs.loc[all_collection_obs["dataset_id"] == dataset_id, "dataset_title"].values[0]
        new_title = (
            dataset_title.replace(" ", "_")
            .replace("\\/", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("\\", "")
            .replace("'", "")
            .replace(":", "")
            .replace(";", "")
            .replace("&", "")
        )
    #check if file exists
    # if f"original_{new_title}.h5ad") exists, skip
    # else download
        if not os.path.exists(os.path.join(outdir,f"original_{new_title}.h5ad")): 
            cellxgene_census.download_source_h5ad(
                dataset_id=dataset_id, 
                to_path=os.path.join(outdir,f"original_{new_title}.h5ad"), 
                census_version=census_version
            )

    #for dataset_title in all_collection_obs["dataset_title"].unique():
        og = sc.read_h5ad(os.path.join(outdir,f"original_{new_title}.h5ad"))
        #write to file
        og.obs.to_csv(os.path.join(outdir,f"original_{new_title}.obs.tsv"), sep="\t")


if __name__ == "__main__":
    main()
