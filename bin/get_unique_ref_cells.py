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
dimport seaborn as sns
import json
import argparse
import os
import json
import yaml


# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
  #parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
      #"Transcriptomic cytoarchitecture reveals principles of human neocortex organization",
      #"SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas",
      #"Molecular and cellular evolution of the primate dorsolateral prefrontal cortex"
  #]),
  parser.add_argument("--organ", type=str, default="brain")
  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args
        
def main():
  # Parse command line arguments
  args = parse_arguments()

  # Set organism and census_version from arguments
  organism = args.organism
  census_version = args.census_version
  #ref_collections = args.ref_collections
  organ = args.organ

  census = cellxgene_census.open_soma(census_version=census_version)
  dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
  cellxgene_obs = adata_functions.get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
  
  cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
  cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
  #cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)] 

  unique_ref_cells = cellxgene_obs[["cell_type","cell_type_ontology_term_id","collection_name"]].value_counts().reset_index()
  unique_ref_cells.to_csv(f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_info/{census_version}/unique_cells_{organism}.tsv", sep="\t", index=False)
  
if __name__ == "__main__":
    main()