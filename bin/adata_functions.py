import pandas as pd
import numpy as np
import scanpy as sc
import random
import cellxgene_census
import cellxgene_census.experimental
import os
import anndata as ad
import scvi
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import *
from sklearn.preprocessing import label_binarize
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
#current_directory = Path.cwd()
projPath = "/space/grp/rschwartz/rschwartz/cpsc545_proj/"

import subprocess


def setup(organism="homo_sapiens", version="2024-07-01"):
    organism=organism.replace(" ", "_") 
    #census = cellxgene_census.open_soma(census_version=version)
    outdir = f"scvi-{organism}-{version}"  # Concatenate strings using f-string
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Check if the model file exists
    model_file_path = os.path.join(outdir, "model.pt")
    #if not os.path.exists(model_file_path):
        # Get scVI model info
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
            embedding_name="scvi",
            organism=organism,
            census_version=version,
        )

        # Extract the model link
    model_link = scvi_info["model_link"]
    date = model_link.split("/")[5]
    url = os.path.join("https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/", date, organism, "model.pt")

    # Download the model using wget if it doesn't already exist
    subprocess.run(["wget", "--no-check-certificate", "-q", "-O", model_file_path, url])
# else:
     #   print(f"File already exists at {model_file_path}, skipping download.")

    return(outdir)

def subsample_and_save(dataset_path, n_cells=1000):
    dataset = ad.read_h5ad(dataset_path)
    subsampled_dataset = dataset[np.random.choice(dataset.n_obs, size=n_cells, replace=False), :] if dataset.n_obs > n_cells else dataset
    dir_name, base_name = os.path.split(dataset_path)
    file_name, ext = os.path.splitext(base_name)
    output_path = os.path.join(dir_name, f"{file_name}_subsampled_{n_cells}{ext}")
    subsampled_dataset.write_h5ad(output_path)
    print(f"Saved to {output_path}")
    
# Subsample x cells from each cell type if there are n>x cells present
#ensures equal representation of cell types in reference
def subsample_cells(data, filtered_ids, subsample=500, relabel_path="/biof501_proj/meta/relabel/census_map_human.tsv", 
                    ref_keys=["rachel_subclass","rachel_class","rachel_family"], seed=42):
    random.seed(seed)         # For `random`
    np.random.seed(seed)      # For `numpy`
    scvi.settings.seed = seed # For `scvi`
    
    # Filter data based on filtered_ids
    obs = data[data['soma_joinid'].isin(filtered_ids)]
    relabel_df = pd.read_csv(relabel_path, sep='\t')  # Adjust the separator as needed
    # Take the first column as the join key
    join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in obs.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Perform the left join to update the metadata
    obs = obs.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
    celltypes = obs[ref_keys[0]].unique()
    final_idx = []
    for celltype in celltypes:
        celltype_ids = obs[obs[ref_keys[0]] == celltype]['soma_joinid'].values
        # Sample if there are enough observations, otherwise take all
        if len(celltype_ids) > subsample:
            subsampled_cell_idx = random.sample(list(celltype_ids), subsample)
        else:
            subsampled_cell_idx = celltype_ids.tolist()
        # Append subsampled indices to final list
        final_idx.extend(subsampled_cell_idx)

    # Return final indices
    return final_idx

def relabel(adata, relabel_path, join_key="", sep="\t"):
    # Read the relabel table from the file
    relabel_df = pd.read_csv(relabel_path, sep=sep)  # Adjust the separator as needed
    # Take the first column as the join key
    if join_key=="":
        join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in adata.obs.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Perform the left join to update the metadata
    adata.obs = adata.obs.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
    columns_to_drop = [col for col in adata.obs.columns if col.endswith('_y')]
    adata.obs.drop(columns=columns_to_drop, inplace=True)
    return adata


def extract_data(data, filtered_ids, subsample=10, organism=None, census=None, 
    obs_filter=None, cell_columns=None, dataset_info=None, dims=20, relabel_path="biof501_proj/meta/relabel/census_map_human.tsv'", 
    ref_keys=["rachel_subclass","rachel_class","rachel_family"], seed=42):
    
    brain_cell_subsampled_ids = subsample_cells(data, filtered_ids, subsample, relabel_path=relabel_path, ref_keys=ref_keys, seed=seed)
    # Assuming get_seurat is defined to return an AnnData object
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_value_filter=obs_filter,  # Ensure this is constructed correctly
        obs_column_names=cell_columns,
        obs_coords=brain_cell_subsampled_ids,
        var_value_filter = "nnz > 10",
        obs_embeddings=["scvi"])
    
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_genes(adata, min_counts=200)
    
    print("Subsampling successful.")
    newmeta = adata.obs.merge(dataset_info, on="dataset_id", suffixes=(None,"y"))
    adata.obs = newmeta
    # Assuming relabel_wrapper is defined
    adata = relabel(adata, relabel_path=relabel_path, join_key="cell_type", sep='\t')
    # Convert all columns in adata.obs to factors (categorical type in pandas)
    return adata

def split_and_extract_data(data, split_column, subsample=500, organism=None, census=None, 
                           cell_columns=None, dataset_info=None, dims=20, relabel_path="/biof501_proj/meta/relabel/census_map_human.tsv",
                           ref_keys=["rachel_subclass","rachel_class","rachel_family"], seed=42):
    # Get unique split values from the specified column
    unique_values = data[split_column].unique()
    refs = {}

    for split_value in unique_values:
        # Filter the brain observations based on the split value
        filtered_ids = data[data[split_column] == split_value]['soma_joinid'].values
        obs_filter = f"{split_column} == '{split_value}'"
        
        adata = extract_data(data, filtered_ids, subsample, organism, census, obs_filter, 
                             cell_columns, dataset_info, dims=dims, relabel_path=relabel_path, ref_keys=ref_keys, seed=seed)
        dataset_titles = adata.obs['dataset_title'].unique()

        if split_column == "tissue": 
            name_to_use = split_value
        elif split_column == "dataset_id":
            name_to_use = dataset_titles[0]
        else:
            name_to_use = split_value

        refs[name_to_use] = adata

    return refs

def get_filtered_obs(census, organism, organ="brain", is_primary=True, disease="normal"):
    value_filter = (
        f"tissue_general == '{organ}' and "
        f"is_primary_data == {str(is_primary)} and "
        f"disease == '{disease}'"
    )
    return cellxgene_census.get_obs(census, organism, value_filter=value_filter)

def get_census(census_version="2024-07-01", organism="homo_sapiens", subsample=5, split_column="dataset_id", dims=50, organ="brain",
               ref_collections=["Transcriptomic cytoarchitecture reveals principles of human neocortex organization"," SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas"],
               relabel_path=f"{projPath}meta/census_map_human.tsv", seed=42, ref_keys=["rachel_subclass","rachel_class","rachel_family"]):

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    brain_obs = get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
    
    brain_obs = brain_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    brain_obs.drop(columns=['soma_joinid_y'], inplace=True)
    brain_obs_filtered = brain_obs[brain_obs['collection_name'].isin(ref_collections)] 
    # Filter based on organism
    #if organism == "homo_sapiens":
        #brain_obs_filtered = brain_obs_filtered[~brain_obs_filtered['cell_type'].isin(["unknown", "glutamatergic neuron"])] # remove non specific cells
    #elif organism == "mus_musculus":
        #brain_obs_filtered = brain_obs_filtered[~brain_obs_filtered['cell_type'].isin(["unknown",
                                                                                        #"hippocampal neuron", 
                                                                                        #"cortical interneuron", 
                                                                                        #"glutamatergic neuron",
                                                                                        #"GABAergic neuron"])]
    #else:
       #raise ValueError("Unsupported organism")

    # Adjust organism naming for compatibility
    organism_name_mapping = {
        "homo_sapiens": "Homo sapiens",
        "mus_musculus": "Mus musculus"
    }
    organism = organism_name_mapping.get(organism, organism)

    cell_columns = [
        "assay", "cell_type", "tissue",
        "tissue_general", "suspension_type",
        "disease", "dataset_id", "development_stage",
        "soma_joinid","observation_joinid"
    ]
    # Get individual datasets and embeddings
    refs = split_and_extract_data(
        brain_obs_filtered, split_column=split_column,
        subsample=subsample, organism=organism,
        census=census, cell_columns=cell_columns,
        dataset_info=dataset_info, dims=dims,
        relabel_path=relabel_path,
        ref_keys=ref_keys, seed=seed
    )
    # Get embeddings for all data together
    filtered_ids = brain_obs_filtered['soma_joinid'].values
    adata = extract_data(
        brain_obs_filtered, filtered_ids,
        subsample=subsample, organism=organism,
        census=census, obs_filter=None,
        cell_columns=cell_columns, dataset_info=dataset_info, dims=dims,
        relabel_path=relabel_path, ref_keys=ref_keys, seed = seed
    )
    refs["whole cortex"] = adata
    for name, ref in refs.items():
        dataset_title = name.replace(" ", "_")
        for col in ref.obs.columns:
            if ref.obs[col].dtype.name =='category':
    # Convert to Categorical and remove unused categories
                ref.obs[col] = pd.Categorical(ref.obs[col].cat.remove_unused_categories())
    
    return refs



def process_query(query, model_file_path, batch_key="sample"):
    # Ensure the input AnnData object is valid
    if not isinstance(query, ad.AnnData):
        raise ValueError("Input must be an AnnData object.")

    # Assign ensembl_id to var
    #query.var["ensembl_id"] = query.var["feature_id"]
    if "feature_id" in query.var.columns:
        query.var.set_index("feature_id", inplace=True)

    query.obs["n_counts"] = query.X.sum(axis=1)
    query.obs["joinid"] = list(range(query.n_obs))
    query.obs["batch"] = query.obs[batch_key]

    # Filter out missing HGNC features
    #query = query[:, query.var["feature_name"].notnull().values].copy()

    # Prepare the query AnnData for scVI
    scvi.model.SCVI.prepare_query_anndata(query, model_file_path)
    vae_q = scvi.model.SCVI.load_query_data(query, model_file_path)

    # Set the model to trained and get latent representation
    vae_q.is_trained = True
    latent = vae_q.get_latent_representation()
    query.obsm["scvi"] = latent

    return query


def map_valid_labels(query, ref_keys, mapping_df):
    # deal with differing levels of granularity
    for key in ref_keys:
        original=query[key].unique()
        for og in original:
            if og not in mapping_df[key].unique(): # this will cause a problem in mouse
                # need to create dummy labels for mouse
                level = mapping_df.columns[mapping_df.apply(lambda col: og in col.values, axis=0)]
                # get the highest level in the hierarchy 
                if level.empty: # handle cases where original label is not in mapping file- should only be "unknown"
                    # this doesn't handle individual references-
                    # if a reference is missing a label that is present in other references, it will have 0 predictions
                    # this will therefore penalize references with missing labels
                    continue
                og_index = query.index[query[key] == og]
                # Replace the value in "predicted_" column with corresponding predicted value at `level`
                for idx in og_index:
                    # Find the replacement value from `mapping_df` for this level
                    replacement = query.loc[idx, "predicted_" + level]
                    # replace predicted id with appropriate level
                    query["predicted_" + key] = query["predicted_" + key].astype("object")
                    query.loc[idx, "predicted_" + key] = replacement.iloc[0]
                    query["predicted_" + key] = query["predicted_" + key].astype("category")

    return query            






def rfc_pred(ref, query, ref_keys, seed):
    """
    Fit a RandomForestClassifier at the most granular level and aggregate probabilities for higher levels.
    
    Parameters:
    - ref: Reference data with labels.
    - query: Query data for prediction.
    - ref_keys: List of ordered keys from most granular to highest level (e.g., ["rachel_subclass", "rachel_class", "rachel_family"]).
    - tree: Dictionary representing the hierarchy of classes.
    
    Returns:
    - probabilities: Dictionary with probabilities for each level of the hierarchy.
    """
    probabilities = {}
    
    # The most granular key is the first in the ordered list
    granular_key = ref_keys[0]
    
    # Initialize and fit the RandomForestClassifier at the most granular level
    rfc = RandomForestClassifier(class_weight='balanced', random_state=seed)
    rfc.fit(ref.obsm["scvi"], ref.obs[granular_key].values)
    # Predict probabilities at e most granular level
    probs_granular = rfc.predict_proba(query.obsm["scvi"])
    class_labels_granular = rfc.classes_
    query.obs[granular_key] = query.obs[granular_key].astype(str)
    base_score = rfc.score(query.obsm["scvi"], query.obs[granular_key].values)

    # Store granular level probabilities
    probabilities[granular_key] = {
        "probabilities": probs_granular,
        "class_labels": class_labels_granular,
        "accuracy": base_score
    }
    
    return probabilities


def roc_analysis(probabilities, query, key):
    optimal_thresholds = {}
    metrics={}
  #  for key in ref_keys:
       # print(key) 
    #probs = probabilities[key]["probabilities"]
    probs = np.array(probabilities)
    class_labels = probabilities.columns.values
    optimal_thresholds[key] = {}
        
    # Binarize the class labels for multiclass ROC computation
    true_labels = label_binarize(query[key].values, classes=class_labels)
        
    # Find the optimal threshold for each class
    metrics[key] = {}
    for i, class_label in enumerate(class_labels):
        optimal_thresholds[key][class_label] = {}
        # check for positive samples
        # usually positive samples are 0 when a ref label is
        # replaced with a parent label
        # since it is not in the original query labels
        # but it is being annotated during the label transfer
        # these should not be evaluated ?
        # incorrect classifications will be reflected in the AUC and F1 of the og label
        # eg. ET is not in query so it is changed to "deep layer non-IT"
        # but these cells are CT or NP in the ref, so they're still incorrect
        # not entirely sure how to deal with this
        positive_samples = np.sum(true_labels[:, i] == 1)
        if positive_samples == 0:
            print(f"Warning: No positive samples for class {class_label}, skipping eval and setting threshold to 0.5")
            optimal_thresholds[key][class_label] = 0.5
        elif positive_samples > 0:
            metrics[key][class_label]={}
            # True label one hot encoding at class label index = 
            # vector of all cells which are either 1 = label or 0 = not label
            # probs = probability vector for all cells given class label
            fpr, tpr, thresholds = roc_curve(true_labels[:, i], probs[:, i])
            roc_auc = auc(fpr, tpr) # change to roc_auc_score, ovo, average= macro, labels               
            optimal_idx = np.argmax(tpr - fpr)
            optimal_threshold = thresholds[optimal_idx]
            if optimal_threshold == float('inf'):
                optimal_threshold = 0 
            optimal_thresholds[key][class_label]=optimal_threshold
            metrics[key][class_label]["tpr"] = tpr
            metrics[key][class_label]["fpr"] = fpr
            metrics[key][class_label]["auc"] = roc_auc
            metrics[key][class_label]["optimal_threshold"] = optimal_threshold
    return metrics


def pr_analysis(prob_df, query, key, mapping_df):

    metrics = defaultdict(dict)

    probs = np.array(prob_df)
    class_labels = prob_df.columns.values         
    original_labels = query[key].unique()
    
    # map back to correct granularity
    #for label in class_labels:
        #if label not in original_labels:
            
            #parent_label = mapping_df.loc[mapping_df[key] == label, key].values[0]
            #if parent_label:
                #probs[:, class_labels == label] = probs[:, class_labels == parent
                #_label]
    
    # Binarize the class labels for multiclass ROC computation 
    
    true_labels = label_binarize(query[key].values, classes=class_labels)
    
    # Find the optimal threshold for each class

    for i, class_label in enumerate(class_labels):
        positive_samples = np.sum(true_labels[:, i] == 1)
        if positive_samples == 0:
            print(f"Warning: No positive samples for class {class_label}, skipping eval")
        elif positive_samples > 0:
            precision, recall, thresholds = precision_recall_curve(true_labels[:, i], probs[:, i])
            avg_precision = average_precision_score(true_labels[:, i], probs[:, i])

            # Find optimal threshold as the one maximizing F1 score
            f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
            optimal_idx = np.argmax(f1_scores)

            # Store metrics
            metrics["precision"][class_label] = precision
            metrics["recall"][class_label] = recall
            metrics["average_precision"][class_label] = avg_precision
            metrics["optimal_threshold"][class_label] = thresholds[optimal_idx]
    return metrics

def check_column_ties(probabilities, class_labels):
    """
    Checks for column ties (multiple columns with the same maximum value) in each row.

    Parameters:
    - probabilities (numpy.ndarray): 2D array where rows represent samples and columns represent classes.
    - class_labels (list): List of class labels corresponding to the columns.

    Returns:
    - tie_rows (list): List of row indices where ties occur.
    - tie_columns (dict): Dictionary where the key is the row index and the value is a list of tied column indices and their class labels.
    """
    # Find the maximum probability for each row
    max_probs = np.max(probabilities, axis=1)

    # Check if there are ties (multiple columns with the maximum value in the row)
    ties = np.sum(probabilities == max_probs[:, None], axis=1) > 1

    # Get the indices of rows with ties
    tie_rows = np.where(ties)[0]

    # Find the columns where the tie occurs and associated class labels
    tie_columns = {}
    for row in tie_rows:
        tied_columns = np.where(probabilities[row] == max_probs[row])[0]
        tie_columns[row] = [(col, class_labels[col]) for col in tied_columns]
    
    return tie_rows, tie_columns

def classify_cells(query, ref_keys, cutoff, probabilities, mapping_df):
    class_metrics = {}
    
    # Only use the first ref_key
    key = ref_keys[0]
    class_metrics[key] = {}

    # Extract the class labels and probabilities (DataFrame structure)
    class_labels = probabilities.columns.values  # Class labels are the column names
    class_probs = probabilities.values  # Probabilities as a numpy array
    
    predicted_classes = []
    
    if cutoff > 0:
        # Find the class with the maximum probability for each cell
        max_class_indices = np.argmax(class_probs, axis=1)  # Get the index of the max probability
        max_class_probs = np.max(class_probs, axis=1)  # Get the max probability
        
        # Set predicted classes to "unknown" if the max probability does not meet the threshold
        predicted_classes = [
            class_labels[i] if prob > cutoff else "unknown"
            for i, prob in zip(max_class_indices, max_class_probs)
        ]
    else:
        # Direct prediction without threshold filtering
        predicted_classes = class_labels[np.argmax(class_probs, axis=1)]
    
    # Store predictions and confidence in `query`
    query["predicted_" + key] = predicted_classes
    query["confidence"] = np.max(class_probs, axis=1)  # Store max probability as confidence
    
    # Aggregate predictions (you can keep this logic as needed)
    query = aggregate_preds(query, ref_keys, mapping_df)
    
    return query


def aggregate_preds(query, ref_keys, mapping_df):
    preds = np.array(query["predicted_" + ref_keys[0]])
    query.index = query.index.astype(int)

    # Reorder ref_keys based on the order in tree_df columns (ignoring "cell_type")
    ref_keys = [col for col in mapping_df.columns if col in ref_keys]

    for higher_level_key in ref_keys[1:]:  # Skip the first (most granular) level
        # Get mapping from subclass to higher-level class
        mapping = mapping_df.set_index(ref_keys[0])[higher_level_key].to_dict()

        # Assign higher-level labels based on mapping
        query["predicted_" + higher_level_key] = query["predicted_" + ref_keys[0]].map(mapping)

        # Fill NA values with original subclass labels
        query["predicted_" + higher_level_key] = query["predicted_" + higher_level_key].fillna(query["predicted_" + ref_keys[0]])

    return query

def eval(query, ref_keys, mapping_df):
    class_metrics = defaultdict(lambda: defaultdict(dict))
    for key in ref_keys:
        
       #threshold = kwargs.get('threshold', True)  # Or some other default value    
        query = map_valid_labels(query, ref_keys, mapping_df)  
        class_labels = query[key].unique()
        pred_classes = query[f"predicted_{key}"].unique()
        true_labels= query[key].astype(str)
        predicted_labels = query["predicted_" + key].astype(str)
        labels = list(set(class_labels).union(set(pred_classes)))

    # Calculate accuracy and confusion matrix after removing "unknown" labels
        accuracy = accuracy_score(true_labels, predicted_labels)
        # Add accuracy to classification report
        class_metrics[key]["accuracy"] = accuracy
        
        ## Calculate accuracy for each label
              
        conf_matrix = confusion_matrix(
            true_labels, predicted_labels, 
            labels=labels
        )
        class_metrics[key]["confusion"] = {
            "matrix": conf_matrix,
            "labels": labels
            #"accuracy": accuracy
        }
            # Compute per-label precision, recall, F1-score, and support
        precision, recall, f1, support = precision_recall_fscore_support(
            true_labels, predicted_labels, labels=labels, zero_division=np.nan
        )
        # Compute weighted average metrics
        avg_precision, avg_recall, avg_f1, _ = precision_recall_fscore_support(
            true_labels, predicted_labels, average="weighted", zero_division=np.nan
        )

        # Compute per-label accuracy and store all metrics
        label_metrics = {}
        for i, label in enumerate(labels):
            label_mask = true_labels == label

            # Handle case where there are no true instances
            # if label mask is all False
            if not label_mask.any():
                label_accuracy = "nan"
                precision[i] = "nan"
                recall[i] = "nan"
                f1[i] = "nan"
                label_accuracy = "nan"
            else:
                label_accuracy = accuracy_score(true_labels[label_mask], predicted_labels[label_mask])

            label_metrics[label] = {
                "accuracy": label_accuracy,
                "precision": precision[i],
                "recall": recall[i],
                "f1_score": f1[i],
                "support": support[i]
            } 

        # Store per-label metrics
        class_metrics[key]["label_metrics"] = label_metrics

        # Store weighted averages
        class_metrics[key]["weighted_metrics"] = {
            "precision": avg_precision,
            "recall": avg_recall,
            "f1_score": avg_f1
        } 

    return class_metrics

def update_classification_report(class_metrics, ref_keys):
    #for query_name, query in class_metrics.items():
       # for ref_name, ref in query.items():
    for key in ref_keys:
        for key, val in class_metrics[key]["classification_report"].items():
            if isinstance(val, dict):
                if val['support'] == 0.0: 
                    val['f1-score'] = "nan" 
    return class_metrics

def plot_confusion_matrix(query_name, ref_name, key, confusion_data, output_dir):

    #new_query_name = query_name.replace(" ", "_").replace("/", "_").replace("(","").replace(")","")
   # new_ref_name = ref_name.replace(" ", "_").replace("/", "_").replace("(","").replace(")","")
 #   output_dir = os.path.join(new_query_name, new_ref_name, output_dir)
    os.makedirs(output_dir, exist_ok=True) 
    # Extract confusion matrix and labels from the confusion data
    conf_matrix = confusion_data["matrix"]
    labels = confusion_data["labels"]

    # Plot the confusion matrix
    plt.figure(figsize=(28, 12))
    sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Reds', xticklabels=labels, yticklabels=labels, annot_kws={"size": 20})
    plt.xlabel('Predicted', fontsize =20)
    plt.ylabel('True', fontsize= 20)
    plt.title(f'Confusion Matrix: {query_name} vs {ref_name} - {key}', fontsize=17)
    # Adjust tick parameters for both axes
    #plt.tick_params(axis='both', which='major', labelsize=15, width=1)  # Increase tick label size and make ticks thicker

    # Rotate both x and y tick labels by 90 degrees
    plt.xticks(rotation=45, fontsize=15)  # Rotate x-axis labels by 90 degrees
    plt.yticks(rotation=45, fontsize=15)  # Rotate y-axis labels by 90 degrees

    # Save the plot
   # output_dir = os.path.join(projPath, 'results', 'confusion')
    
    #os.makedirs(os.path.join(output_dir, new_query_name, new_ref_name), exist_ok=True)  # Create the directory if it doesn't exist
    plt.savefig(os.path.join(output_dir,f"{key}_confusion.png"))
    plt.close() 

def plot_roc_curves(metrics, title="ROC Curves for All Keys and Classes", save_path=None):
    """
    Plots ROC curves for each class at each key level from the metrics dictionary on the same figure.
    
    Parameters:
    metrics (dict): A dictionary with structure metrics[key][class_label] = {tpr, fpr, auc, optimal_threshold}.
    title (str): The title of the plot.
    save_path (str, optional): The file path to save the plot. If None, the plot is not saved.
    """
    fig, ax = plt.subplots(figsize=(10, 8))  # Create a figure and axis

    # Create a subplot for each key
    for key in metrics:
        for class_label in metrics[key]:

            if isinstance(metrics[key][class_label], dict):
                if all(k in metrics[key][class_label] for k in ["tpr", "fpr", "auc"]):
                    tpr = metrics[key][class_label]["tpr"]
                    fpr = metrics[key][class_label]["fpr"]
                    roc_auc = metrics[key][class_label]["auc"]

                    # Find the index of the optimal threshold
                    optimal_idx = np.argmax(tpr - fpr)
                    optimal_fpr = fpr[optimal_idx]
                    optimal_tpr = tpr[optimal_idx]

                    # Plot the ROC curve for the current class
                    #plt.plot(fpr, tpr, lw=2, label=f"Class {class_label} (AUC = {roc_auc:.3f})")
                 #   curve = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc, estimator_name=class_label) #drop_intermediate=False)
                   # curve.plot(ax=ax)  # Add to the shared axis
                    
                    ax.step(fpr, tpr, where='post', lw=2, label=f"{key}: {class_label} (AUC = {roc_auc:.3f})")

                    # Plot the optimal threshold as a point
                    ax.scatter(optimal_fpr, optimal_tpr, color='red', marker='o') 
                          #  label=f"Optimal Threshold (Class {class_label})")
                    
    # Plot the reference line (random classifier)
    ax.plot([0, 1], [0, 1], 'k--', lw=2, label="Random Classifier")

    # Add title, labels, legend, and grid
    ax.set_title(title, fontsize=16)
    ax.set_xlabel('False Positive Rate', fontsize = 15)
    ax.set_ylabel('True Positive Rate', fontsize = 15)
    ax.legend(loc='lower right', bbox_to_anchor=(1.05, 0), fontsize='medium', borderaxespad=0)
    ax.grid(True)

    # Adjust layout and save the plot if a path is provided
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
        plt.close()
    
    
def combine_f1_scores(class_metrics, ref_keys):
   # metrics = class_metrics
    # Dictionary to store DataFrames for each key
    all_f1_scores = {}
    #cutoff = class_metrics["cutoff"]
    # Iterate over each key in ref_keys
    for key in ref_keys:
        # Create a list to store F1 scores for each query-ref combo
        f1_data = [] 
        # Iterate over all query-ref combinations
        for query_name in class_metrics:
            for ref_name in class_metrics[query_name]:
                # Extract the classification report for the current query-ref-key combination
                classification_report = class_metrics[query_name][ref_name][key]["classification_report"]
                # Extract F1 scores for each label
                if classification_report:
                    for label, metrics in classification_report.items():
                        if label not in ["macro avg","micro avg","weighted avg","accuracy"]:
                         #   if isinstance(metrics, dict) and 'f1-score' in metrics:
                                f1_data.append({
                                    'query': query_name,
                                    'reference': ref_name,
                                    'label': label,
                                    'f1_score': metrics['f1-score'],                         
                                    'macro_f1': classification_report.get('macro avg', {}).get('f1-score', None),
                                    'micro_f1': classification_report.get('micro avg', {}).get('f1-score', None),
                                    'weighted_f1': classification_report.get('weighted avg', {}).get('f1-score', None), #,
                                    'precision': metrics['precision'],
                                    'recall': metrics['recall']        
                                })

        # Create DataFrame for the current key
        df = pd.DataFrame(f1_data)

        # Store the DataFrame in the dictionary for the current key
        all_f1_scores[key] = df

    return all_f1_scores



def plot_label_f1_heatmaps(all_f1_scores, threshold, outpath, widths=[1,0.8,0.5], acronym_mapping=None):
    """
    Plot horizontally stacked heatmaps for label-level F1 scores for each query across multiple keys with variable widths,
    shared y-axis labels, a shared color bar, and a title for the query.

    Parameters:
        all_f1_scores (dict): Dictionary with keys as reference names and values as DataFrames containing F1 scores.
        threshold (float): Threshold value to display in plot titles.
        outpath (str): Directory to save the generated heatmaps.
        widths (list of float): Proportional widths for subplots. If None, defaults to equal widths.
    """
    sns.set(style="whitegrid")
    os.makedirs(outpath, exist_ok=True)

    # Determine global F1 score range
    all_scores = pd.concat([df['f1_score'] for df in all_f1_scores.values()])
    vmin, vmax = all_scores.min(), all_scores.max()
    
    # Initialize an empty set to store unique queries
    queries = set() 
    keys = list(all_f1_scores.keys()) # get levels of hierarchy
    for key, df in all_f1_scores.items():
        queries.update(df['query'].unique())  # add unique queries to the set
    queries = sorted(queries) # sort queries

    if widths is None:
        widths = [1] * len(keys)  # Equal widths by default for each level

    for query in queries:
        # Create a figure with variable subplot widths
        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(keys),
            figsize=(sum(widths) * 15, 8),
            gridspec_kw={'width_ratios': widths},
            constrained_layout=True
        )

        # Add a figure title for the query
        fig.suptitle(f'Class-level F1 for Query: {query}\nThreshold = {threshold:.2f}', fontsize=25, y=1.1)

        if len(keys) == 1:
            axes = [axes]  # Ensure axes is always iterable

        for i, key in enumerate(keys):
            if key not in all_f1_scores:
                continue
            
            df = all_f1_scores[key]
            query_df = df[df['query'] == query]

            # Pivot DataFrame to create the heatmap
            pivot_df = query_df.pivot_table(index='reference_acronym', columns='label', values='f1_score')
            mask = pivot_df.isnull() | (pivot_df == "nan")

            sns.heatmap(
                pivot_df,
                annot=True,
                cmap='YlOrRd',
                cbar=i == len(keys) - 1,  # Add cbar only for the last subplot
                cbar_kws={'label': 'F1 Score'} if i == len(keys) - 1 else None,
                mask=mask,
                ax=axes[i],
                linewidths=0.5,
                annot_kws={"size": 8},
                vmin=vmin,  # Use global vmin
                vmax=vmax   # Use global vmax
            )

            for i, key in enumerate(keys):
                axes[i].set_title(f'{key}', fontsize=25)
                
                # Set x-axis tick labels and rotation
                axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=90, fontsize=20)

                # Only add y-axis labels to the leftmost subplot
                if i == 0:
                    axes[i].set_ylabel('Reference', fontsize=25)
                    axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=20)  # Set font size for y-axis labels
                    if acronym_mapping:
                        legend_text = "\n".join([f"{k}: {v}" for k, v in acronym_mapping.items()])
                        axes[i].text(
                            0.1, 0.5, legend_text,  # Position text to the left of the subplot
                            fontsize=14,
                            verticalalignment='center',
                            horizontalalignment='right',
                            transform=axes[i].transAxes,  # Use subplot's coordinate system
                            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.1')
                        )
                else:
                    axes[i].set_ylabel("", fontsize=11)
                    axes[i].set_yticks([])  # Remove y-axis labels for other subplots
                    # Save the figure
        plt.savefig(os.path.join(outpath, f'f1_heatmaps_{query}_threshold_{threshold:.2f}.png'), bbox_inches='tight')
        plt.close()


def plot_f1_heatmaps_by_level(weighted_f1_data, threshold, outpath, ref_keys, acronym_mapping=None):
    """
    Plot one heatmap for each level of the hierarchy with references as rows and queries as columns.

    Parameters:
        weighted_f1_data (pd.DataFrame): DataFrame containing F1 scores with columns:
                                         'level', 'reference', 'query', 'weighted_f1'.
        threshold (float): Threshold value to display in plot titles.
        outpath (str): Directory to save the heatmaps.
        levels (list): List of levels to plot (e.g., ['Rachel_class', 'Rachel_subclass', 'Rachel_family']).
        ref_keys (list): List of reference keys to ensure consistent ordering of rows.
    """
    sns.set(style="whitegrid")
    os.makedirs(outpath, exist_ok=True)

    for level in ref_keys:
        # Filter data for the current level
        level_data = weighted_f1_data[weighted_f1_data['key'] == level]

        # Pivot the data for heatmap
        pivot_f1 = level_data.pivot_table(
            index='reference_acronym',
            columns='query',
            values='weighted_f1'
        )

        # Define color map limits
        vmin = weighted_f1_data['weighted_f1'].min()
        vmax = weighted_f1_data['weighted_f1'].max()

        # Create the heatmap
        plt.figure(figsize=(20, 15))
        sns.heatmap(
            pivot_f1,
            annot=True,
            cmap='YlOrRd',
            cbar_kws={'label': 'Weighted F1 Score'},
            fmt='.3f',
            annot_kws={"size": 15},
            vmin=vmin,
            vmax=vmax
        )

        # Add title and axis labels
        plt.title(f'Weighted F1 Scores for {level}\nThreshold = {threshold:.2f}', fontsize=25)
        plt.ylabel('Reference', fontsize=25)
        plt.xlabel('Query', fontsize=25)
        plt.xticks(rotation=90, ha='right', fontsize=25)
        plt.yticks(fontsize=25, rotation=90)
        
        if acronym_mapping:
            # Add an annotation box with the acronym legend
            legend_text = "\n".join([f"{k}: {v}" for k, v in acronym_mapping.items()])
            plt.gcf().text(
                0.85, 0.5, legend_text, 
                fontsize=14, 
                verticalalignment='center', 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.1')
            )
        
        plt.tight_layout()
        # Save the heatmap
        plt.savefig(os.path.join(outpath, f'{level}_weighted_f1_heatmap_threshold_{threshold:.2f}.png'))
        plt.close()

    
 
def split_anndata_by_obs(adata, obs_key="dataset_title"):
    """
    Split an AnnData object into multiple AnnData objects based on unique values in an obs key.

    Parameters:
    - adata: AnnData object to split.
    - obs_key: Key in `adata.obs` on which to split the data.

    Returns:
    - A dictionary where keys are unique values in `obs_key` and values are corresponding AnnData subsets.
    """
    # Dictionary comprehension to create a separate AnnData for each unique value in obs_key
    split_data = {
        value: adata[adata.obs[obs_key] == value].copy() 
        for value in adata.obs[obs_key].unique()
    }
    
    return split_data


