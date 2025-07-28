#!/user/bin/python3

import warnings
warnings.filterwarnings("ignore")
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
import cellxgene_census
import cellxgene_census.experimental
import scvi
from sklearn.ensemble import RandomForestClassifier
import utils
from utils import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
import yaml


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2025-01-30', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--subsample_ref', type=int, default=50)
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "Transcriptomic cytoarchitecture reveals principles of human neocortex organization",
        "SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas",
        "Molecular and cellular evolution of the primate dorsolateral prefrontal cortex",
        "HVS: Human variation study"
    ]) 
    parser.add_argument('--split_column', type=str, default="dataset_id")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family"])
    parser.add_argument('--seed', type=int, default=42)
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        
     
def create_ref_region_yaml(refs, outdir):
    ref_names = list(refs.keys())
    ref_region_yaml = {}
    for ref_name in ref_names:
        unique_regions = refs[ref_name].obs["tissue"].unique().tolist()
        
        # If there are multiple region types, handle them by including them as a list
        if len(unique_regions) == 1:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = unique_regions[0]
        else:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = "multiple regions"
    
    with open(os.path.join(outdir, "ref_region.yaml"), 'w') as file:
        yaml.dump(ref_region_yaml, file)

         
def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set organism and census_version from arguments
    organism = args.organism
    census_version = args.census_version
    subsample_ref = args.subsample_ref
    relabel_path = args.relabel_path
    ref_collections = args.ref_collections
    split_column = args.split_column
    ref_keys = args.ref_keys
    SEED = args.seed

    if organism == "mus_musculus":
        original_celltypes = get_original_celltypes(columns_file=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}/original_celltype_columns.tsv",
                                                    author_annotations_path=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}") 
    else:
        original_celltypes = None

    refs = adata_functions.get_census(
        organism=organism,
        organ="brain",
        subsample=subsample_ref,
        split_column=split_column,
        census_version=census_version,
        relabel_path=relabel_path,
        ref_collections=ref_collections,
        seed=SEED,
        ref_keys=ref_keys,
        original_celltypes = original_celltypes
    )
    
    refs.pop('All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2', None)
    refs.pop('Microglia - 24 months old wild-type and Rag1-KO', None)
    refs.pop('Single-cell of aged oligodendrocytes', None)

    print("finished fetching anndata")
    outdir = "refs"
    os.makedirs(outdir, exist_ok=True)

    for ref_name, ref in refs.items():
        if ref.shape[0] < 50:
            continue

        new_dataset_title = (
            ref_name.replace(" ", "_")
            .replace("\\/", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("\\", "")
            .replace("'", "")
            .replace(":", "")
            .replace(";", "")
            .replace("&", "")
        )
        #if organism == "mus_musculus": 
            #ref.obs["new_dataset_title"] = ref.obs["dataset_title"].apply(lambda x:x.replace(" ", "_")
                                                                            #.replace("\\/", "_")
                                                                            #.replace("(", "")
                                                                            #.replace(")", "")
                                                                            #.replace("\\", "")
                                                                            #.replace("'", "")
                                                                            #.replace(":", "")
                                                                            #.replace(";", "")
                                                                            #.replace("&", "")
                                                                        #)
            #if new_dataset_title in original_celltypes.keys():
                #og_obs = original_celltypes[new_dataset_title]
                #ref.obs["new_observation_joinid"] = ref.obs["new_dataset_title"].astype(str) + "_" + ref.obs["observation_joinid"].astype(str)
                ## merge left, only keep leftmost observation_joinid
                ## instead of merging, create a dict and map using new_observation_joinid
                #mapping = dict(zip(og_obs["new_observation_joinid"], og_obs["author_cell_type"]))
                #ref.obs["author_cell_type"] = ref.obs["new_observation_joinid"].map(mapping)
                #ref.obs[["subclass","cell_type","author_cell_type","dataset_title"]].value_counts().reset_index().to_csv(f"refs/{ref_name}_new_celltypes.tsv", sep="\t", index=False)

                #refs[ref_name] = ref
            
        ref.write(os.path.join(outdir, f"{new_dataset_title}.h5ad"))
        ref.obs.to_csv(os.path.join(outdir, f"{new_dataset_title}.obs.tsv"), sep="\t")

  #  if organism == "mus_musculus":

   #    replace_ambiguous_cells(refs, ambiguous_celltypes)
        
    for ref_name, ref in refs.items():
        sc.pp.neighbors(ref, use_rep="scvi")
        sc.tl.umap(ref)
        colors = ["dataset_title","collection_name","subclass","class"]
        if "author_cell_type" in ref.obs.columns:
            colors = colors + ["author_cell_type"]
        if "new_cell_type" in ref.obs.columns:
            colors = colors + ["new_cell_type"]
        sc.pl.umap(ref, color=colors, ncols=1, save=f"{ref_name}_umap.png")
        


    create_ref_region_yaml(refs, outdir)

      
if __name__ == "__main__":
    main()