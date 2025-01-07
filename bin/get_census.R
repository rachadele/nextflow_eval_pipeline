subsample_cells <- function(data, filtered_ids, census, subsample = 500) {
    # Filter data based on filtered_ids
    obs <- data[data$soma_joinid %in% filtered_ids, ]
    celltypes <- unique(obs$cell_type)

    final_idx <- c()

    for (celltype in celltypes) {
        celltype_ids <- obs[obs$cell_type == celltype, ]$soma_joinid
        # Sample if there are enough observations, otherwise take all
        if (length(celltype_ids) >= subsample) {
            subsampled_cell_idx <- sample(celltype_ids, subsample)
        } else {
            subsampled_cell_idx <- celltype_ids
        }
        # Append subsampled indices to final list
        final_idx <- c(final_idx, subsampled_cell_idx)
    }
    # Return final indices
    return(final_idx)
}


extract_data <- function(data, filtered_ids, subsample=500, organism, census, obs_filter, cell_columns, dataset_info, relabel_path,
      dims=20) {

  brain_cell_subsampled_ids = subsample_cells(data, filtered_ids, census, subsample)

          seurat_obj <- get_seurat(
            census = census,
            organism = organism,
            # var_value_filter = gene_filter, # Uncomment if needed
            obs_value_filter = obs_filter,  # Ensure this is constructed correctly
            obs_column_names = cell_columns,
            obs_coords = brain_cell_subsampled_ids,
            var_index="feature_id"
          )
          print("Subsampling successful.")
      # Subset the Seurat object to keep only genes with HGNC"
          # Filter out genes that are not expressed in at least 3 cells
          seurat_obj <- seurat_obj[rowSums(seurat_obj[["RNA"]]@counts > 0) >= 3, ]
          #seurat_obj <- seurat_obj %>% NormalizeData() %>% 
                    #ScaleData() %>% 
                    #FindVariableFeatures() %>% 
                    #RunPCA(npcs=dims) %>% RunUMAP(dims=1:dims)
          #batch_key="dataset_id"
          #seurat_obj <- harmony_integration(seurat_obj, dims=20, integration_key=batch_key)

          newmeta <- left_join(seurat_obj@meta.data , dataset_info, by="dataset_id")  %>%
              select(-soma_joinid.y, everything()) %>%  # Remove soma_joinid.y and keep all other columns
              rename(soma_joinid = soma_joinid.x) 
          rownames(newmeta) <- colnames(seurat_obj)
          seurat_obj@meta.data  <- newmeta
          seurat_obj <- relabel_wrapper(seurat_obj, relabel_path=relabel_path)
          return(seurat_obj)
 }

# Function to extract data based on a specified split value
split_and_extract_data <- function(data, split_column, subsample=500, organism, census, cell_columns, dataset_info, relabel_path, dims=20) {
  # Get unique split values from the specified column
  unique_values <- unique(data[[split_column]])
  # Initialize an empty list to store results
  refs <- list()
  # Loop through each unique split value
  for (split_value in unique_values) {
    # Filter the brain observations based on the split value
    filtered_ids <- data[data[[split_column]] == split_value, ]$soma_joinid
    # Create a filter expression for the Seurat object
    obs_filter <- paste0(split_column, " == '", split_value, "'")
    # Check if the number of filtered ids is greater than or equal to the subsample size
    seurat_obj <- extract_data(data, filtered_ids, subsample, organism, census, obs_filter, cell_columns, dataset_info, relabel_path, dims=dims)
    # Get the unique dataset title
    dataset_titles <- unique(seurat_obj@meta.data[["dataset_title"]])
    # Determine the name to use based on the number of dataset titles
    if (length(dataset_titles) > 1) {
      # Use the split column value as the name
      name_to_use <- split_value
    } else {
      # Use the unique dataset title
      name_to_use <- dataset_titles
    }
    # Store the Seurat object in the refs list using the determined name
    refs[[name_to_use]] <- seurat_obj
  } 
  return(refs)
}

# Call the function with your parameters
get_census <- function(organism="homo_sapiens", census_version="2024-07-01", subsample=500, split_column="dataset_id", dims=20, relabel_path,
                                      ref_collections=c("Transcriptomic cytoarchitecture reveals principles of human neocortex organization",
                                                                                   "SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas")) {
      census <- open_soma(census_version)
      # Open obs SOMADataFrame
      dataset_info <- census$get("census_info")$get("datasets")$read()$concat() %>% as.data.frame()
      brain_obs <- census$get("census_data")$get("homo_sapiens")$get("obs")$read(
        value_filter = "tissue_general == 'brain' & is_primary_data == True & disease == 'normal'",
        column_names = c("assay", "cell_type", "sex", "tissue", "tissue_general", "suspension_type", "disease","dataset_id",
        "development_stage","soma_joinid")
      )
      # Concatenates results to an Arrow Table
      brain_obs <-  brain_obs$concat() %>% as.data.frame()
      brain_obs <- left_join(brain_obs, dataset_info, by="dataset_id")
      brain_obs <- brain_obs %>%
        select(-soma_joinid.y, everything()) %>%  # Remove soma_joinid.y and keep all other columns
        rename(soma_joinid = soma_joinid.x) 
      # Allen datasets ?
     # if (organism == "homo_sapiens") {

      brain_obs_filtered <- brain_obs[
          brain_obs$collection_name %in% ref_collections, ]

      brain_obs_filtered <- brain_obs_filtered[!brain_obs_filtered$cell_type %in% c("unknown","glutamatergic neuron"), ]

      #} else if (organism == "mus_musculus") {
      #  brain_obs_filtered <- brain_obs[
      #    brain_obs$collection_name %in% c("x"), ]  # Replace "x" with actual collection names for mouse data
     # } else {
     #   stop("Unsupported organism")
     # }

      # Subsampling cells and extracting the soma_joinid
      set.seed(1)
      # if organism="homo_sapiens set to Homo sapiens"
      #if mus_musculus set to Mus musculus
      if (organism=="homo_sapiens"){
          organism="Homo sapiens"
      } else if (organism=="mus_musculus") {
          organism=="Mus musculus"
      }
      cell_columns <-  c("assay", "cell_type", "tissue", 
      "tissue_general", "suspension_type", 
      "disease","dataset_id","development_stage",
      "soma_joinid")
      refs <- list()
      #everything except whole cortex
      refs <- split_and_extract_data(brain_obs_filtered, split_column=split_column, subsample, organism, census, cell_columns, dataset_info, 
        relabel_path=relabel_path, dims=dims)
      filtered_ids = brain_obs_filtered$soma_joinid
      #whole cortex
      seurat_obj <- extract_data(brain_obs_filtered, filtered_ids, subsample, organism, census, obs_filter=NULL, 
        cell_columns, dataset_info, relabel_path,
        dims=dims) 
      refs[["whole cortex"]] <- seurat_obj

     # names(refs) <- c("whole cortex", unique(brain_obs_filtered$tissue))
    lapply(names(refs), function(x){
      dataset_title <- gsub(" ","_",x)
      #p <- DimPlot(refs[[x]], group.by=c("rachel_subclass","assay","tissue","dataset_title"),ncol=2, pt.size=2) +
        #theme(
          #text = element_text(size = 20),  # Adjust the base text size
          #axis.title = element_text(size = 25),  # Adjust axis title size
          #axis.text = element_text(size = 20),  # Adjust axis label size
          #legend.text = element_text(size = 20),  # Adjust legend text size
          #legend.title = element_text(size = 20)  # Adjust legend title size
      #)
      #switch back
      if (organism=="Homo sapiens"){
        organism="homo_sapiens"
      } else if (organism=="Mus musculus") {
        organism=="mus_musculus"
      }
      #ggsave(plot=p, file=file.path(projPath, paste0("/refs/census/",dataset_title,"_",subsample,"_umap.png")), width=35, height=25)
      dir.create("refs", showWarnings = FALSE)
      saveRDS(refs[[x]], file=file.path("refs",paste0(dataset_title,".rds")))
      #meta <- refs[[x]]@meta.data %>% select("cell_type","rachel_class","rachel_subclass","rachel_family") %>% unique()
      #write.table(meta, file=file.path(projPath, paste0("/meta/relabel/",dataset_title,"_relabel.tsv")), sep="\t",row.names=FALSE)
  })
      
#return(refs)
}
#refs <- get_census()



