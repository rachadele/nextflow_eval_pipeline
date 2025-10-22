
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
from types import SimpleNamespace
import yaml
from sklearn.metrics import precision_recall_curve, average_precision_score, PrecisionRecallDisplay
from collections import defaultdict

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/2025-01-30/mus_musculus/100/dataset_id/SCT/gap_false/ref_50_cutoff_0/scvi/GSE152715.2/whole_cortex/GSE152715.2_1052248_GSM4624685/probs/GSE152715.2_1052248_GSM4624685.obs.relabel.tsv")
    parser.add_argument('--ref_name', type=str, default="whole_cortex") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family","global"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/2025-01-30/mus_musculus/100/dataset_id/SCT/gap_false/ref_50_cutoff_0/scvi/GSE152715.2/whole_cortex/GSE152715.2_1052248_GSM4624685/probs/probs/GSE152715.2_1052248_GSM4624685_whole_cortex.prob.df.tsv")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_region_mapping', type=str, default="")
    parser.add_argument('--study_name', type=str, default="GSE152715.2")
    parser.add_argument('--use_gap', action='store_true', help="Use gap analysis for classification")
    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def get_unique_value(df, column, default=None):
    if column in df.columns:
    # check how many unique values there are
        if len(df[column].unique()) == 1:
            return df[column].unique()[0] 
        else:
            return default

    
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()
    query_path = args.query_path
    ref_name = args.ref_name
    ref_keys = args.ref_keys
    cutoff = args.cutoff
    ref_region_mapping = args.ref_region_mapping
    study_name = args.study_name
    if args.use_gap:
        use_gap = True
    else:
        use_gap = False 
    
    # Load data
    ref_region_mapping = yaml.load(open(ref_region_mapping), Loader=yaml.FullLoader)
    ref_region=ref_region_mapping[ref_name]
    
    prob_df = pd.read_csv(args.probs, sep="\t")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
    query_name = os.path.basename(query_path).replace(".obs.relabel.tsv", "")
    query = pd.read_csv(query_path, sep="\t")
    #for factor in factors:
    

    query_region = get_unique_value(query, 'region')
    disease = get_unique_value(query, 'disease')
    sex = get_unique_value(query, 'sex')
    dev_stage = get_unique_value(query, 'dev_stage')
    treatment = get_unique_value(query, 'treatment')
    genotype = get_unique_value(query, 'genotype')
    strain = get_unique_value(query, 'strain')
    age = get_unique_value(query, 'age')

    os.makedirs("pr_curves", exist_ok=True) 
    
    # Classify cells and evaluate
    query = classify_cells(query=query, ref_keys=ref_keys, cutoff=cutoff, probabilities=prob_df, mapping_df=mapping_df, use_gap=use_gap)

    outdir = os.path.join("predicted_meta")
    os.makedirs(outdir, exist_ok=True)

    # map valid labels for given query granularity and evaluate
    query = map_valid_labels(query, ref_keys, mapping_df)  
    class_metrics = evaluate_sample_predictions(query, ref_keys, mapping_df)
    
    query.to_csv(os.path.join(outdir,f"{query_name}_{ref_name}.predictions.{cutoff}.tsv"), index=False, sep="\t")

    #class_metrics = update_classification_report(class_metrics, ref_keys)

    # Plot confusion matrices
    for key in ref_keys:
        outdir = os.path.join("confusion")
        os.makedirs(outdir, exist_ok=True)
        plot_confusion_matrix(query_name, ref_name, key, class_metrics[key]["confusion"], output_dir=outdir)

    # Collect F1 scores
    f1_data = []
    for key in ref_keys:
        label_metrics = class_metrics[key]["label_metrics"]
        weighted_metrics = class_metrics[key]["weighted_metrics"]
        macro_metrics = class_metrics[key]["macro_metrics"]
        nmi = class_metrics[key]["nmi"]
        ari = class_metrics[key]["ari"]
        overall_accuracy = class_metrics[key]["overall_accuracy"]
        
        for label, metrics in label_metrics.items():
            #if label not in ["macro avg", "micro avg", "weighted avg", "accuracy"]:
                f1_data.append({
                    'query': query_name,
                    'study': study_name,
                    'reference': ref_name,
                    'label': label,
                    'f1_score': metrics['f1_score'],
                    'accuracy': metrics['accuracy'],
                    'precision': metrics['precision'],
                    'recall': metrics['recall'],
                    'support': metrics['support'],
                    'weighted_f1': weighted_metrics.get('f1_score', None),
                    'weighted_precision': weighted_metrics.get('precision', None),
                    'weighted_recall': weighted_metrics.get('recall', None),
                    # add macro averages
                    'macro_f1': macro_metrics.get('f1_score', None),
                    'macro_precision': macro_metrics.get('precision', None),
                    'macro_recall': macro_metrics.get('recall', None),
                    # add micro averages
                    'micro_f1': class_metrics[key]["micro_metrics"].get('f1_score', None),
                    'micro_precision': class_metrics[key]["micro_metrics"].get('precision', None),
                    'micro_recall': class_metrics[key]["micro_metrics"].get('recall', None),
                    'nmi': nmi,
                    'ari': ari,
                    'overall_accuracy': overall_accuracy,
                    'key': key,
                    'cutoff': cutoff,
                    'ref_region': ref_region
                }
        )

    # Save F1 scores to a file
    df = pd.DataFrame(f1_data)
    
    fields_dict = {'disease': disease, 'sex': sex, 'dev_stage': dev_stage, 
                   'query_region': query_region, 
                   'treatment': treatment, 'genotype': genotype, 'strain': strain, 'age': age}   
    for field, value in fields_dict.items():
        df[field] = value if value is not None else np.nan


    outdir = "label_transfer_metrics"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.summary.scores.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
    

