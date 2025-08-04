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
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2025-01-30', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Adult mouse cortical cell taxonomy revealed by single cell transcriptomics",
        "Tabula Muris Senis",
        "Single-cell transcriptomics characterization of oligodendrocytes and microglia in white matter aging"
    ]) 
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--organ', type=str, default="brain")
    parser.add_argument('--assay', type=str, nargs = "+", help="Assays to subset from referenc (unnecessary)", default=None)
    parser.add_argument('--tissue', type=str, nargs="+", default = "cortex", help = "tissues to pull from (different from organ, this can select for more specific brain regions)")
    parser.add_argument('--subsample', type=int, help="Number of cells per cell type to subsample from reference", default=50)
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/author_cell_annotations/rename_cells_mmus_author.tsv")
    parser.add_argument('--original_celltype_columns', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2025-01-30/original_celltype_columns.tsv")
    parser.add_argument('--author_annotations_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2025-01-30")
    parser.add_argument('--split_column', type=str, default="dataset_id")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family","global"])
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
    tissue = args.tissue
    assay = args.assay
    
    
    
    
    if organism == "mus_musculus":
        original_celltypes = get_original_celltypes(columns_file=args.original_celltype_columns,
                                                    author_annotations_path= args.author_annotations_path) 
    else:
        original_celltypes = None

    refs = utils.get_census(
        organism=organism,
        organ=organ,
        tissue=tissue,
        assay=assay,
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
                
        ref.write(os.path.join(outdir, f"{new_dataset_title}.h5ad"))
        ref.obs.to_csv(os.path.join(outdir, f"{new_dataset_title}.obs.tsv"), sep="\t") 
        
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