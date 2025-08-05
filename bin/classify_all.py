
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
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/hsap/70/e01ffd3e44041673370e367f2157ac/lim_H5109Cin.obs.relabel.tsv")
    parser.add_argument('--ref_name', type=str, default="Human_Multiple_Cortical_Areas_SMART-seq") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/hsap/70/e01ffd3e44041673370e367f2157ac/lim_H5109Cin_Human_Multiple_Cortical_Areas_SMART-seq_prediction_scores_seurat.tsv")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_region_mapping', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/hsap/06/03e7ac72a1ef3b67ce6a357eebd8c3/refs/ref_region.yaml")
    parser.add_argument('--study_name', type=str, default="lim")
    
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
    
    #pr_metrics = pr_analysis(prob_df, query, key=ref_keys[0], mapping_df = mapping_df)
    #if not pr_metrics:
        #print("No PR metrics to plot")
        #fig, ax = plt.subplots(figsize=(5, 5))  # Create a small dummy figure
        #fig.savefig(os.path.join("pr_curves", f"{query_name}_{ref_name}_pr_curve.png"), bbox_inches='tight')
        #plt.close(fig)  # Close to free memory
    #else:    
        ## if pr metrics is not empty:
        #_, ax = plt.subplots(figsize=(10, 15))
        #plt.rcParams.update({'font.size': 20})
        
        #f_scores = np.linspace(0.2, 0.8, num=4)
        #lines, labels = [], []
        #for f_score in f_scores:
            #x = np.linspace(0.01, 1)
            #y = f_score * x / (2 * x - f_score)
            #(l,) = plt.plot(x[y >= 0], y[y >= 0], color="gray", alpha=0.2)
            #plt.annotate("f1={0:0.1f}".format(f_score), xy=(0.9, y[45] + 0.02))
            
        #class_labels= pr_metrics["recall"].keys()
        
        ## transform to array
        #class_labels = np.array(list(class_labels))
        ## make colors for each class
        #colors = sns.color_palette("husl", len(class_labels))
        #for class_label, color in zip(class_labels, colors):

            #display = PrecisionRecallDisplay(recall=pr_metrics["recall"][class_label], 
                                            #precision=pr_metrics["precision"][class_label], 
                                            #average_precision=pr_metrics["average_precision"][class_label])
            #display.plot(ax=ax, name=f"Precision-recall for class {class_label}", color=color)
            #avg_prec = pr_metrics["average_precision"][class_label]
            #opt_thresh = pr_metrics["optimal_threshold"][class_label]
            ## get position of class label
            #class_idx = np.where(class_labels == class_label)[0][0]
        ## plt.annotate(f"Optimal threshold for {class_label}={opt_thresh:.2f}", 
            ##           xy=(0.6, 0.2 + class_idx * 0.05), fontsize=12, color=color)

        ## Add the legend
        #handles, labels = display.ax_.get_legend_handles_labels()
        #handles.extend([l])
        #labels.extend(["iso-f1 curves"])
        #ax.legend(handles=handles, labels=labels, loc="best", bbox_to_anchor=(1.3, 1), borderaxespad=0., fontsize =12)

        #ax.set_xlabel('Recall', fontsize=16)
        #ax.set_ylabel('Precision', fontsize=16)
        #ax.set_title(f"Precision-Recall Curve for {query_name} vs {ref_name}", fontsize=18)

        #plt.savefig(os.path.join("pr_curves",f"{query_name}_{ref_name}_pr_curve.png"), bbox_inches='tight')
        #plt.close()

         
    
    # Classify cells and evaluate
    query = classify_cells(query, ref_keys, cutoff=cutoff, probabilities=prob_df, mapping_df=mapping_df)

    outdir = os.path.join("predicted_meta")
    os.makedirs(outdir, exist_ok=True)

    # map valid labels for given query granularity and evaluate
    query = map_valid_labels(query, ref_keys, mapping_df)  
    class_metrics = eval(query, ref_keys, mapping_df)
    
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
                    'key': key,
                    'cutoff': cutoff,
                    'ref_region': ref_region,
                }
        )
                
    #celltype_precisions = class_metrics["precision"]
    #celltype_recalls = class_metrics["recall"]

    # Save F1 scores to a file
    df = pd.DataFrame(f1_data)
    
    fields_dict = {'disease': disease, 'sex': sex, 'dev_stage': dev_stage, 
                   'query_region': query_region, 
                   'treatment': treatment, 'genotype': genotype, 'strain': strain, 'age': age}   
    for field, value in fields_dict.items():
        df[field] = value if value is not None else np.nan


    outdir = "f1_results"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.f1.scores.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
    

