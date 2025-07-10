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
import adata_functions
from adata_functions import *
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
        "Molecular and spatial signatures of mouse brain aging at single-cell resolution"
    ]) 
    parser.add_argument('--organ', type=str, default="brain")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        
def main():
  # Parse command line arguments
  args = parse_arguments()
  # Set organism and census_version from arguments
  organism = args.organism
  census_version = args.census_version
  ref_collections=args.ref_collections
  os.makedirs(f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_info/{census_version}", exist_ok=True)
  
  #ref_collections = args.ref_collections
  organ = args.organ
  if organism == "mus_musculus":
    original_celltypes = get_original_celltypes(columns_file=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}/original_celltype_columns.tsv",
                                                      author_annotations_path=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}") 
  else:
    original_celltypes = None
          
  census = cellxgene_census.open_soma(census_version=census_version)
  dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
  dataset_info.to_csv(f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_info/{census_version}/{organism}_dataset_info.tsv", sep="\t", index=False)
  cellxgene_obs = adata_functions.get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
  
  cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
  cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
  # add author_cell_type to obs
  # this will enable proper relabeling and subsampling
  # need to add it back in after getting ids
  cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)] 

  if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
      cellxgene_obs_filtered = map_author_labels(cellxgene_obs_filtered, original_celltypes)
      cell_type_info = cellxgene_obs_filtered[["author_cell_type", "cell_type", "dataset_title", "collection_name","cell_type_ontology_term_id"]].value_counts().reset_index()
      cell_type_info.to_csv(f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_info/{census_version}/{organism}_cell_type_info.tsv", sep="\t", index=False)
  else:
      cell_type_info = cellxgene_obs_filtered[["cell_type", "dataset_title", "collection_name","cell_type_ontology_term_id"]].value_counts().reset_index()
      cell_type_info.to_csv(f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_info/{census_version}/{organism}_cell_type_info.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()