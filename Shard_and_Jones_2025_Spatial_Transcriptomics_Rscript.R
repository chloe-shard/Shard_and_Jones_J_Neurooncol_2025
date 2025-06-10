################################# Spatial Transcriptomics Analysis ########################################################################################

#Load libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

#Load samples
YY199_annot <- readRDS("path to file /YY199_annot.rds") #GBM
YY200_annot <- readRDS("path to file /YY200_annot.rds") #GBM
YY201_annot <- readRDS("path to file /YY201_annot.rds") #GBM
YY202_annot <- readRDS("path to file /YY202_annot.rds") #GBM
DMG1_annot <- readRDS("path to file /DMG1_annot.rds") #DIPG
DMG2_annot <- readRDS("path to file /DMG2_annot.rds") #DIPG
DMG3_annot <- readRDS("path to file /DMG3_annot.rds") #DIPG
DMG4_annot <- readRDS("path to file /DMG4_annot.rds") #DIPG
DMG5_annot <- readRDS("path to file /DMG5_annot.rds") #DIPG

# Log normalize each sample
DefaultAssay(DMG1_annot) <- "Spatial"
DMG1_annot <- NormalizeData(DMG1_annot, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE, assay = "Spatial")

###### annotate spatial regions using signatures from Puchalski et al., 2018 (Science) - Fig.1E, 3A, and 4A ########

# Load the signature gene lists
Puchalski_LE_sig <- read_csv("path to file /Puchalski_LE_sig.csv")
Puchalski_CTpan_sig <- read_csv("path to file/Puchalski_CTpan_sig.csv")
Puchalski_CTmvp_sig <- read_csv("path to file/Puchalski_CTmvp_sig.csv")
Puchalski_CT_sig <- read_csv("path to file/Puchalski_CT_sig.csv")

# Extract the gene lists
LE_genes <- Puchalski_LE_sig[["LE"]]
CTmvp_genes <- Puchalski_CTmvp_sig[["CTmvp"]]
CTpan_genes <- Puchalski_CTpan_sig[["CTpan"]]
CT_genes <- Puchalski_CT_sig[["CT"]]

# Function to add region metadata to a Seurat object
add_tumour_region_metadata <- function(seurat_object, LE_genes, CTmvp_genes, CTpan_genes, CT_genes) {
  all_genes <- rownames(seurat_object)
  
  LE_genes <- LE_genes[LE_genes %in% all_genes]
  CTmvp_genes <- CTmvp_genes[CTmvp_genes %in% all_genes]
  CTpan_genes <- CTpan_genes[CTpan_genes %in% all_genes]
  CT_genes <- CT_genes[CT_genes %in% all_genes]
  
  if(length(LE_genes) == 0) print("No LE genes found")
  if(length(CTmvp_genes) == 0) print("No CTmvp genes found")
  if(length(CTpan_genes) == 0) print("No CTpan genes found")
  if(length(CT_genes) == 0) print("No CT genes found")
  
  LE_genes_features <- FetchData(object = seurat_object, vars = LE_genes, slot = "data")
  LE_genes_features <- t(LE_genes_features)
  
  CTmvp_genes_features <- FetchData(object = seurat_object, vars = CTmvp_genes, slot = "data")
  CTmvp_genes_features <- t(CTmvp_genes_features)
  
  CTpan_genes_features <- FetchData(object = seurat_object, vars = CTpan_genes, slot = "data")
  CTpan_genes_features <- t(CTpan_genes_features)
  
  CT_genes_features <- FetchData(object = seurat_object, vars = CT_genes, slot = "data")
  CT_genes_features <- t(CT_genes_features)
  
  LE_avg_Score <- t(colMeans(LE_genes_features, na.rm = TRUE))
  LE_avg_Score <- t(LE_avg_Score)
  colnames(LE_avg_Score)[1] = "LE_region"
  
  CTmvp_avg_Score <- t(colMeans(CTmvp_genes_features, na.rm = TRUE))
  CTmvp_avg_Score <- t(CTmvp_avg_Score)
  colnames(CTmvp_avg_Score)[1] = "CTmvp_region"
  
  CTpan_avg_Score <- t(colMeans(CTpan_genes_features, na.rm = TRUE))
  CTpan_avg_Score <- t(CTpan_avg_Score)
  colnames(CTpan_avg_Score)[1] = "CTpan_region"
  
  CT_avg_Score <- t(colMeans(CT_genes_features, na.rm = TRUE))
  CT_avg_Score <- t(CT_avg_Score)
  colnames(CT_avg_Score)[1] = "CT_region"
  
  Tumour_region <- cbind(LE_avg_Score, CTmvp_avg_Score, CTpan_avg_Score, CT_avg_Score)
  Tumour_region <- as.data.frame(Tumour_region)
  
  seurat_object <- AddMetaData(seurat_object, Tumour_region)
  return(seurat_object)
}

# List of Seurat objects
seurat_objects <- list(YY199_annot, YY200_annot, YY201_annot, YY202_annot, DMG1_annot, DMG2_annot, DMG3_annot, DMG4_annot, DMG5_annot)

# Apply the function to each Seurat object in the list
seurat_objects <- lapply(seurat_objects, add_tumour_region_metadata, LE_genes = LE_genes, CTmvp_genes = CTmvp_genes, CTpan_genes = CTpan_genes, CT_genes = CT_genes)

# Assign the updated objects back to their original names
YY199_annot <- seurat_objects[[1]]
YY200_annot <- seurat_objects[[2]]
YY201_annot <- seurat_objects[[3]]
YY202_annot <- seurat_objects[[4]]
DMG1_annot <- seurat_objects[[5]]
DMG2_annot <- seurat_objects[[6]]
DMG3_annot <- seurat_objects[[7]]
DMG4_annot <- seurat_objects[[8]]
DMG5_annot <- seurat_objects[[9]]

# Print metadata to check
print(seurat_objects[[1]]@meta.data)

# Spatial region plots
SpatialFeaturePlot(object = DMG1_annot, features = c("LE_region", "CTmvp_region", "CTpan_region", "CT_region"))

############# Manders' coefficent analysis (Fig.3B,3D,4B, 4D, 5C) ##############

# Load packages
library(openxlsx)

# Define genes of interest [e.g. GABA genes or synapse genes]
genes_of_interest <- c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", 
                       "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3", 
                       "GABRR1", "GABRR2", "GABRR3", "GABRD", "GABRQ", "GABRP", 
                       "GABRE")

# List of Seurat objects
#Option 1 DIPG
seurat_objects <- list(DMG1 = DMG1_annot, DMG2 = DMG2_annot, DMG3 = DMG3_annot, DMG4 = DMG4_annot, DMG5 = DMG5_annot)

#Option 2 GBM
#seurat_objects <- list(YY199 = YY199_annot, YY200 = YY200_annot, YY201 = YY201_annot, YY202 = YY202_annot)

# Create workbook
wb <- createWorkbook()

# Initialize summary info table
slides_info <- data.frame(slide_name = character(), Data_range = character(), stringsAsFactors = FALSE)

# First add 'slides_info' sheet (empty for now)
addWorksheet(wb, "slides_info")

# Loop over each Seurat object
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  
  # Coordinates
  coords <- GetTissueCoordinates(seurat_obj)[, c("imagecol", "imagerow")]
  colnames(coords) <- c("x", "y")
  
  # Metadata
  meta_data <- seurat_obj@meta.data[, c("LE_region", "CTmvp_region", "CTpan_region", "CT_region")]
  
  # Gene expression
  gene_expr <- FetchData(seurat_obj, vars = genes_of_interest)
  
  # Combine
  output_df <- cbind(coords, meta_data, gene_expr)
  output_df$spot <- rownames(output_df)
  output_df <- output_df[, c("spot", setdiff(colnames(output_df), "spot"))]
  
  # Write data to sheet
  addWorksheet(wb, sample_name)
  writeData(wb, sheet = sample_name, x = output_df)
  
  # Calculate data range (excluding header row)
  n_rows <- nrow(output_df) + 1
  n_cols <- ncol(output_df)
  end_col <- int2col(n_cols)
  data_range <- paste0("A2:", end_col, n_rows)
  
  # Record in slides_info
  slides_info <- rbind(slides_info, data.frame(slide_name = sample_name, Data_range = data_range))
}

# Now write 'slides_info' data
writeData(wb, sheet = "slides_info", x = slides_info)

# Save the workbook - 
saveWorkbook(wb, file = "path to output folder/GABA_Receptor_SpatialData.xlsx", overwrite = TRUE)

#use the above workbook as input for Manders' coefficient analysis (generate a separate workbook for GBM and DMG samples)
# 1) first run through scale correction factor pipeline! 
# 2) once the scale is corrected then run Manders' pipeline

############# Heatmap - neuronal module expression in tumour regions Fig.1I, 1J ##############

###Threshold each region signature

# Define the regions
Regions <- c("LE_region", "CTmvp_region", "CTpan_region", "CT_region")

process_seurat_object <- function(spatial_obj) {
  for (Region in Regions) {
    # Convert the column to numeric (if not already)
    numeric_values <- as.numeric(as.character(spatial_obj@meta.data[[Region]]))
    
    # Check if conversion resulted in any non-NA numeric values
    if (all(is.na(numeric_values))) {
      warning(paste("Column", Region, "has no valid numeric values. Skipping."))
      next
    }
    
    # Calculate the mean and standard deviation of the region
    mean_value <- mean(numeric_values, na.rm = TRUE)
    sd_value <- sd(numeric_values, na.rm = TRUE)
    
    # Calculate the threshold (mean + 1 standard deviation)
    threshold <- mean_value + sd_value
    
    # Create new annotation column
    new_annotation <- paste0(Region, "_thresholded")
    spatial_obj@meta.data[[new_annotation]] <- ifelse(numeric_values > threshold, "positive", "negative")
  }
  return(spatial_obj)
}

# List of Seurat objects and their names
seurat_names <- c("YY199_annot", "YY200_annot", "YY201_annot", "YY202_annot", "DMG1_annot", "DMG2_annot", "DMG3_annot", "DMG4_annot", "DMG5_annot")

# Apply the function to each Seurat object and assign back to their original names
for (seurat_name in seurat_names) {
  # Get the Seurat object by name
  spatial_obj <- get(seurat_name)
  
  # Process the Seurat object
  processed_object <- process_seurat_object(spatial_obj)
  
  # Assign the processed object back to its original name
  assign(seurat_name, processed_object)
}

# Print metadata to check
print(DMG1_annot@meta.data)

#spatial plot
SpatialDimPlot(DMG1_annot, group.by = "CT_region_thresholded", pt.size.factor = 3)

##### Average expression in tumour regions

# Option 1 DIPG datasets
seurat_list <- list(DMG1 = DMG1_annot, DMG2 = DMG2_annot, DMG3 = DMG3_annot, DMG4 = DMG4_annot, DMG5 = DMG5_annot)

# Option 2 GBM datasets
#seurat_list <- list(YY199 = YY199_annot, YY200 = YY200_annot, YY201 = YY201_annot, YY202 = YY202_annot)

# Metadata region columns
region_columns <- c("LE_region_thresholded", "CTmvp_region_thresholded", "CTpan_region_thresholded", "CT_region_thresholded")

Genesofinterest <- Genesofinterest_GBM #Neuronal regulation module gene list
# Genesofinterest <- c("GABRA1", "GABRB2", ...)  # example

# Store results for all samples
all_results <- list()

# Loop through each Seurat object
for (sample_name in names(seurat_list)) {
  seurat_obj <- seurat_list[[sample_name]]
  
  # Get expression matrix (assumed to be log-normalised)
  expr <- as.data.frame(t(FetchData(seurat_obj, vars = Genesofinterest)))
  
  # Add region annotations
  meta <- seurat_obj@meta.data[, region_columns, drop = FALSE]
  expr <- cbind(meta, t(expr))
  
  # Prepare to collect results per region
  sample_results <- list()
  
  for (region in region_columns) {
    # Subset to "positive" spots for the region
    region_df <- expr %>% filter(.data[[region]] == "positive")
    
    # Only keep gene columns (remove region annotation columns)
    region_gene_df <- region_df[, !(names(region_df) %in% region_columns)]
    
    # Calculate mean expression per gene
    avg_expr <- region_gene_df %>%
      summarise(across(everything(), mean, na.rm = TRUE)) %>%
      t() %>%
      as.data.frame()
    
    colnames(avg_expr) <- paste0(sample_name, "_", region)
    avg_expr$Gene <- rownames(avg_expr)
    rownames(avg_expr) <- NULL
    
    sample_results[[region]] <- avg_expr
  }
  
  # Merge all regions for this sample
  sample_merged <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), sample_results)
  
  # Min-max normalisation for each gene (row) across the regions (columns 2:n)
  numeric_cols <- colnames(sample_merged)[-1]
  sample_merged[numeric_cols] <- t(apply(sample_merged[numeric_cols], 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) {
      return(rep(0, length(x))) # No variation, return 0s
    } else {
      return((x - rng[1]) / (rng[2] - rng[1]))
    }
  }))
  
  # Store in all_results
  all_results[[sample_name]] <- sample_merged
}

# Merge across all samples
final_output <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), all_results)

# View result
head(final_output)

write.csv(final_output, "path to output folder/final_output.csv", row.names=TRUE)

