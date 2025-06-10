################################# Single Cell Analysis ########################################################################################

###### Heatmaps for Fig 1C and 1D [Neuronal regulation module - cell type expression]

#Load necessary libraries
library(Seurat)
library(readr)
library(Seurat)
library(dplyr)

#load neuronal module gene lists (DIPG = black, GBM = blue)
GABA_black_DIPG_module_genes <- read_csv("path to file/GABA_black_DIPG_module_genes.csv")
GABA_blue_GBM_module_genes <- read_csv("path to file /GABA_blue_GBM_module_genes.csv")

Genesofinterest_DIPG <- GABA_black_DIPG_module_genes[["genes"]]
Genesofinterest_GBM <- GABA_blue_GBM_module_genes[["genes"]]

# MIN/MAX Norm [use this for heatmaps - Fig.1C, 1D]
normalize_function <- function(x) {
  # Remove NA values for min/max calculation
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  # Normalize only if max_val is greater than min_val
  if (max_val > min_val) {
    (x - min_val) / (max_val - min_val)
  } else {
    # If max_val equals min_val, return NA for all entries
    rep(NA, length(x))
  }
}

# Z Score Norm [use this for overlap analysis - Sup.Fig.1B]
normalize_function <- function(x) {
  # Remove NA values for min/max calculation
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  Mean <- mean(x, na.rm = TRUE)
  SD <- sd(x, na.rm = TRUE)
  # Normalize only if max_val is greater than min_val
  if (max_val > min_val) {
    (x - Mean) / (SD)
  } else {
    # If max_val equals min_val, return NA for all entries
    rep(NA, length(x))
  }
}

# DIPG datasets ###########################

#option 1 neuronal regulation module genes [Fig.1C]
Genesofinterest <- Genesofinterest_DIPG

#option 2 GABA genes [Fig.2C]
#Genesofinterest <- c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3", "GABRR1", "GABRR2", "GABRR3", "GABRD", "GABRQ", "GABRP", "GABRE")

##### Filbin ############

obj_Filbin <- readRDS("path to file/Filbin.rds")

obj_Filbin.df <- as.matrix(GetAssayData(object = obj_Filbin[["RNA"]], slot = "data"))

celltype_Filbin <- obj_Filbin@meta.data$Type
celltype_Filbin <- as.data.frame(celltype_Filbin)

Genesofinterest_Filbin <- subset(obj_Filbin.df, subset = rownames(obj_Filbin.df) %in% Genesofinterest)
Genesofinterest_Filbin <- t(Genesofinterest_Filbin)

Genesofinterest_Filbin_celltype_merge <- data.frame(celltype_Filbin, Genesofinterest_Filbin)

library(dplyr)

Average_Geneofinterest_per_Filbin_celltype <- Genesofinterest_Filbin_celltype_merge %>%
  group_by(celltype_Filbin) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column (make sure the length matches the number of rows in 'data')
new_values <- c("Other", "Macrophages", "Malignant_cells", "Oligodendrocytes") # Update this list as needed
new_values_V2 <- c("Macrophages", "Malignant_cells", "Oligodendrocytes") # Update this list as needed

# Assign new values to the first column
Average_Geneofinterest_per_Filbin_celltype <- Average_Geneofinterest_per_Filbin_celltype[,-1]
rownames(Average_Geneofinterest_per_Filbin_celltype) <- new_values

Average_Geneofinterest_per_Filbin_celltype <- Average_Geneofinterest_per_Filbin_celltype[rownames(Average_Geneofinterest_per_Filbin_celltype) != "Other", ]

# Apply the normalization across columns (genes)
Filbin_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_Filbin_celltype, normalize_function))
rownames(Filbin_normalized_data_celltype) <- new_values_V2

##### Liu ############

obj_Liu <- readRDS("path to file/Liu.rds")

Idents(obj_Liu) <- "age"
obj_Liu <- subset(obj_Liu, idents = c("pediatric"))

obj_Liu.df <- as.matrix(GetAssayData(object = obj_Liu[["RNA"]], slot = "data"))

celltype_Liu <- obj_Liu@meta.data$annotation
celltype_Liu <- as.data.frame(celltype_Liu)
celltype_Liu[celltype_Liu == "AC-like"] <- "malignant"
celltype_Liu[celltype_Liu == "S"] <- "malignant"
celltype_Liu[celltype_Liu == "OPC-like-1"] <- "malignant"
celltype_Liu[celltype_Liu == "OPC-like-2"] <- "malignant"
celltype_Liu[celltype_Liu == "OPC-like-3"] <- "malignant"
celltype_Liu[celltype_Liu == "OC-like"] <- "malignant"
celltype_Liu[celltype_Liu == "G2M"] <- "malignant"
celltype_Liu[celltype_Liu == "MES-like"] <- "malignant"

Genesofinterest_Liu <- subset(obj_Liu.df, subset = rownames(obj_Liu.df) %in% Genesofinterest)
Genesofinterest_Liu <- t(Genesofinterest_Liu)

Genesofinterest_Liu_celltype_merge <- data.frame(celltype_Liu, Genesofinterest_Liu)

library(dplyr)

Average_Geneofinterest_per_Liu_celltype <- Genesofinterest_Liu_celltype_merge %>%
  group_by(celltype_Liu) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column (make sure the length matches the number of rows in 'data')
new_values <- c("Endothelial_cells", "Microglia", "Macrophages", "Oligodendrocytes", "T_cells","Other1", "Malignant_cells", "Other2") # Update this list as needed
new_values_V2 <- c("Endothelial_cells", "Microglia", "Macrophages", "Oligodendrocytes", "T_cells", "Malignant_cells", "Other2") # Update this list as needed
new_values_V3 <- c("Endothelial_cells", "Microglia", "Macrophages", "Oligodendrocytes", "T_cells", "Malignant_cells") # Update this list as needed

# Assign new values to the first column
Average_Geneofinterest_per_Liu_celltype <- Average_Geneofinterest_per_Liu_celltype[,-1]
rownames(Average_Geneofinterest_per_Liu_celltype) <- new_values

Average_Geneofinterest_per_Liu_celltype <- Average_Geneofinterest_per_Liu_celltype[rownames(Average_Geneofinterest_per_Liu_celltype) != "Other1", ]
rownames(Average_Geneofinterest_per_Liu_celltype) <- new_values_V2
Average_Geneofinterest_per_Liu_celltype <- Average_Geneofinterest_per_Liu_celltype[rownames(Average_Geneofinterest_per_Liu_celltype) != "Other2", ]

# Apply the normalization across columns (genes)
Liu_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_Liu_celltype, normalize_function))
rownames(Liu_normalized_data_celltype) <- new_values_V3

# Merge normalized dataframes from all DIPG datasets ###########################

# Dataset names you define
Filbin <- Filbin_normalized_data_celltype
Liu <- Liu_normalized_data_celltype

Filbin$cell_type <- c("Macrophages_Filbin", "Malignant_cells_Filbin", "Oligodendrocytes_Filbin")
Filbin <- Filbin %>% select(cell_type, everything())

Liu$cell_type <- c("Endothelial_cells_Liu", "Microglia_Liu", "Macrophages_Liu", "Oligodendrocytes_Liu", "T_cells_Liu", "Malignant_cells_Liu")
Liu <- Liu %>% select(cell_type, everything())

# Load necessary libraries
library(dplyr)
library(tidyr)

# List of dataframes
data_frames <- list(Filbin, Liu)

# Function to transform a dataset
transform_data <- function(data) {
  # Transform data: first column should be the cell type, remaining columns are gene expressions
  # Gather all columns except the first one to create a long format dataframe
  data_long <- pivot_longer(data, cols = -1, names_to = "Gene", values_to = "Expression")
  
  # Spread the gene expression values to have cell types as columns and genes as rows
  data_wide <- pivot_wider(data_long, names_from = 1, values_from = Expression)
  
  # Ensure the gene column is the first column for merging
  data_wide <- data_wide %>% select(Gene, everything())
  
  return(data_wide)
}

# Apply the function to each dataset and merge them
DIPG_combined_data <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)
}, lapply(data_frames, transform_data))

# Replace any leftover column name duplications (from merging)
names(DIPG_combined_data) <- make.unique(names(DIPG_combined_data))

# Keep the "Gene" column first, sort the rest alphabetically
sorted_columns <- c("Gene", sort(setdiff(names(DIPG_combined_data), "Gene")))

# Reorder the columns based on the sorted list
DIPG_combined_data <- DIPG_combined_data[, sorted_columns]

write.csv(DIPG_combined_data, file = "Path to output folder/DIPG_GeneExpression_Normalized.csv", row.names = FALSE)

#Use DIPG_combined_data for heatmap using Morpheus available at: https://software.broadinstitute.org/morpheus/

######## GBM datasets

#option 1 neuronal regulation module genes [Fig.1D]
Genesofinterest <- Genesofinterest_GBM

#option 2 GABA genes [Fig.2D]
#Genesofinterest <- c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3", "GABRR1", "GABRR2", "GABRR3", "GABRD", "GABRQ", "GABRP", "GABRE")

#Ebert####

obj_Ebert <- readRDS("path to file/Ebert.rds")

obj_Ebert.df <- as.matrix(GetAssayData(object = obj_Ebert[["RNA"]], slot = "data"))

celltype_Ebert <- obj_Ebert@meta.data$combined_clusters_annot
celltype_Ebert <- as.data.frame(celltype_Ebert)

Genesofinterest_Ebert <- subset(obj_Ebert.df, subset = rownames(obj_Ebert.df) %in% Genesofinterest)
Genesofinterest_Ebert <- t(Genesofinterest_Ebert)

Genesofinterest_Ebert_celltype_merge <- data.frame(celltype_Ebert, Genesofinterest_Ebert)

library(dplyr)

Average_Geneofinterest_per_Ebert_celltype <- Genesofinterest_Ebert_celltype_merge %>%
  group_by(celltype_Ebert) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column
new_values <- c("Endothelial_cells", "Macrophages", "Pericytes", "T_cells", "Malignant_cells", "Prolif_immune_cells")

# Assign new values to the first column
Average_Geneofinterest_per_Ebert_celltype <- Average_Geneofinterest_per_Ebert_celltype[,-1]
rownames(Average_Geneofinterest_per_Ebert_celltype) <- new_values


# Apply the normalization across columns (genes)
Ebert_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_Ebert_celltype, normalize_function))
rownames(Ebert_normalized_data_celltype) <- new_values

# LeBlanc ####

obj_LeBlanc <- readRDS("path to file/LeBlanc_tissue.rds")

Idents(obj_LeBlanc) <- "sample"
table(obj_LeBlanc$sample)
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1"))

obj_LeBlanc.df <- as.matrix(GetAssayData(object = obj_LeBlanc[["RNA"]], slot = "data"))

celltype_LeBlanc <- obj_LeBlanc@meta.data$cell_type
celltype_LeBlanc <- as.data.frame(celltype_LeBlanc)

Genesofinterest_LeBlanc <- subset(obj_LeBlanc.df, subset = rownames(obj_LeBlanc.df) %in% Genesofinterest)
Genesofinterest_LeBlanc <- t(Genesofinterest_LeBlanc)

Genesofinterest_LeBlanc_celltype_merge <- data.frame(celltype_LeBlanc, Genesofinterest_LeBlanc)

library(dplyr)

Average_Geneofinterest_per_LeBlanc_celltype <- Genesofinterest_LeBlanc_celltype_merge %>%
  group_by(celltype_LeBlanc) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column
new_values <- c("Endothelial_cells", "Pericytes", "Macrophages", "Malignant_cells", "Neurons", "Oligodendrocytes", "Other")
new_values_V2 <- c("Endothelial_cells", "Pericytes", "Macrophages", "Malignant_cells", "Neurons", "Oligodendrocytes")

# Assign new values to the first column
Average_Geneofinterest_per_LeBlanc_celltype <- Average_Geneofinterest_per_LeBlanc_celltype[,-1]
rownames(Average_Geneofinterest_per_LeBlanc_celltype) <- new_values

Average_Geneofinterest_per_LeBlanc_celltype <- Average_Geneofinterest_per_LeBlanc_celltype[rownames(Average_Geneofinterest_per_LeBlanc_celltype) != "Other", ]


# Apply the normalization across columns (genes)
LeBlanc_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_LeBlanc_celltype, normalize_function))
rownames(LeBlanc_normalized_data_celltype) <- new_values_V2


# Abdelfattah ####
obj_Abdelfattah <- readRDS("path to file/Abdelfattah.rds")

Idents(obj_Abdelfattah) <- "Type"
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("GBM"))

obj_Abdelfattah.df <- as.matrix(GetAssayData(object = obj_Abdelfattah[["RNA"]], slot = "data"))

celltype_Abdelfattah <- obj_Abdelfattah@meta.data$Assignment
celltype_Abdelfattah <- as.data.frame(celltype_Abdelfattah)

Genesofinterest_Abdelfattah <- subset(obj_Abdelfattah.df, subset = rownames(obj_Abdelfattah.df) %in% Genesofinterest)
Genesofinterest_Abdelfattah <- t(Genesofinterest_Abdelfattah)

Genesofinterest_Abdelfattah_celltype_merge <- data.frame(celltype_Abdelfattah, Genesofinterest_Abdelfattah)

library(dplyr)

Average_Geneofinterest_per_Abdelfattah_celltype <- Genesofinterest_Abdelfattah_celltype_merge %>%
  group_by(celltype_Abdelfattah) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column (make sure the length matches the number of rows in 'data')
new_values <- c("B_cells", "Endothelial_cells",  "Malignant_cells", "Macrophages", "Oligodendrocytes", "Other", "Pericytes", "T_cells") # Update this list as needed
new_values_V2 <- c("B_cells", "Endothelial_cells",  "Malignant_cells", "Macrophages", "Oligodendrocytes", "Pericytes", "T_cells") # Update this list as needed

# Assign new values to the first column
Average_Geneofinterest_per_Abdelfattah_celltype <- Average_Geneofinterest_per_Abdelfattah_celltype[,-1]
rownames(Average_Geneofinterest_per_Abdelfattah_celltype) <- new_values

Average_Geneofinterest_per_Abdelfattah_celltype <- Average_Geneofinterest_per_Abdelfattah_celltype[rownames(Average_Geneofinterest_per_Abdelfattah_celltype) != "Other", ]

# Apply the normalization across columns (genes)
Abdelfattah_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_Abdelfattah_celltype, normalize_function))
rownames(Abdelfattah_normalized_data_celltype) <- new_values_V2

# Chen ####

obj_Chen <- readRDS("path to file/Chen.rds")

#remove recurrent patient
Idents(obj_Chen) <- "orig.ident"
obj_Chen <- subset(obj_Chen, idents = c("PDC001", "PJ052", "PJ053", "PW016-703", "PW017-703", "PW032-710", "PW035-710", "PW039-705"))

obj_Chen.df <- as.matrix(GetAssayData(object = obj_Chen[["RNA"]], slot = "data"))

celltype_Chen <- obj_Chen@meta.data$cellkb
celltype_Chen <- as.data.frame(celltype_Chen)

Genesofinterest_Chen <- subset(obj_Chen.df, subset = rownames(obj_Chen.df) %in% Genesofinterest)
Genesofinterest_Chen <- t(Genesofinterest_Chen)

Genesofinterest_Chen_celltype_merge <- data.frame(celltype_Chen, Genesofinterest_Chen)

library(dplyr)

Average_Geneofinterest_per_Chen_celltype <- Genesofinterest_Chen_celltype_merge %>%
  group_by(celltype_Chen) %>%
  summarise_at(vars(everything()), list(name = mean))

# Define new values for the first column (make sure the length matches the number of rows in 'data')
new_values <- c("Malignant_cells", "Macrophages", "Oligodendrocytes", "Endothelial_cells", "Pericytes") # Update this list as needed

# Assign new values to the first column
Average_Geneofinterest_per_Chen_celltype <- Average_Geneofinterest_per_Chen_celltype[,-1]
rownames(Average_Geneofinterest_per_Chen_celltype) <- new_values

# Apply the normalization across columns (genes)
Chen_normalized_data_celltype <- as.data.frame(lapply(Average_Geneofinterest_per_Chen_celltype, normalize_function))
rownames(Chen_normalized_data_celltype) <- new_values


# Merge normalized dataframes from all datasets ###########################

# Dataset names you define
Ebert <- Ebert_normalized_data_celltype
LeBlanc <- LeBlanc_normalized_data_celltype
Abdelfattah <- Abdelfattah_normalized_data_celltype
Chen <- Chen_normalized_data_celltype

Ebert$cellstate <- c("Endothelial_cells_Ebert", "Macrophages_Ebert", "Pericytes_Ebert", "T_cells_Ebert", "Malignant_cells_Ebert", "Prolif_immune_cells_Ebert")
Ebert <- Ebert %>% select(cellstate, everything())

LeBlanc$cellstate <- c("Endothelial_cells_LeBlanc", "Pericytes_LeBlanc", "Macrophages_LeBlanc", "Malignant_cells_LeBlanc", "Neurons_LeBlanc", "Oligodendrocytes_LeBlanc")
LeBlanc <- LeBlanc %>% select(cellstate, everything())

Abdelfattah$cellstate <- c("B_cells_Abdelfattah", "Endothelial_cells_Abdelfattah",  "Malignant_cells_Abdelfattah", "Macrophages_Abdelfattah", "Oligodendrocytes_Abdelfattah", "Pericytes_Abdelfattah", "T_cells_Abdelfattah")
Abdelfattah <- Abdelfattah %>% select(cellstate, everything())

Chen$cellstate <- c("Malignant_cells_Chen", "Macrophages_Chen", "Oligodendrocytes_Chen", "Endothelial_cells_Chen", "Pericytes_Chen")
Chen <- Chen %>% select(cellstate, everything())

# Load necessary libraries
library(dplyr)
library(tidyr)

# List of dataframes
data_frames <- list(Ebert, LeBlanc, Abdelfattah, Chen)

# Function to transform a dataset
transform_data <- function(data) {
  # Transform data: first column should be the cell type, remaining columns are gene expressions
  # Gather all columns except the first one to create a long format dataframe
  data_long <- pivot_longer(data, cols = -1, names_to = "Gene", values_to = "Expression")
  
  # Spread the gene expression values to have cell types as columns and genes as rows
  data_wide <- pivot_wider(data_long, names_from = 1, values_from = Expression)
  
  # Ensure the gene column is the first column for merging
  data_wide <- data_wide %>% select(Gene, everything())
  
  return(data_wide)
}

# Apply the function to each dataset and merge them
GBM_combined_data <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)
}, lapply(data_frames, transform_data))

# Replace any leftover column name duplications (from merging)
names(GBM_combined_data) <- make.unique(names(GBM_combined_data))

# Keep the "Gene" column first, sort the rest alphabetically
sorted_columns <- c("Gene", sort(setdiff(names(GBM_combined_data), "Gene")))

# Reorder the columns based on the sorted list
GBM_combined_data <- GBM_combined_data[, sorted_columns]

write.csv(GBM_combined_data, file = "Path to output folder/GBM_GeneExpression_Normalized.csv", row.names = FALSE)

#Use GBM_combined_data for heatmap using Morpheus available at: https://software.broadinstitute.org/morpheus/

###### select genes that are max in malignant cells

#DIPG (max in both datasets)
# Assuming your dataframe is named combined_data
filtered_genes <- DIPG_combined_data[DIPG_combined_data$Malignant_cells_Liu == 1 & DIPG_combined_data$Malignant_cells_Filbin == 1, ]

# Display the result
print(filtered_genes)

# gene names:
gene_list_DIPG <- filtered_genes$Gene
print(gene_list_DIPG)

#GBM (max in 3 or more datasets)
filtered_genes <- GBM_combined_data %>%
  rowwise() %>%  # Ensures operations are performed per row
  filter(sum(c_across(c(Malignant_cells_Ebert, 
                        Malignant_cells_Abdelfattah, 
                        Malignant_cells_LeBlanc, 
                        Malignant_cells_Chen)) == 1, na.rm = TRUE) >= 3)

# Extract gene list
gene_list_GBM <- filtered_genes$Gene
print(gene_list_GBM)

# Identify overlapping tumour up reg genes between DIPG and GBM
overlapping_genes <- intersect(gene_list_DIPG, gene_list_GBM)

# Display the overlap
print(overlapping_genes)
print(gene_list_GBM)
print(gene_list_DIPG)

########### Mean of tumour cell Z score in DIPG

# Dataset names you define
Filbin <- Filbin_normalized_data_celltype
Liu <- Liu_normalized_data_celltype

library(tibble)

# Subset Malignant_cells from both datasets
Liu_malignant <- as.data.frame(Liu["Malignant_cells", , drop = FALSE])
Filbin_malignant <- as.data.frame(Filbin["Malignant_cells", , drop = FALSE])

# Add row identifiers
Liu_malignant <- Liu_malignant %>% mutate(Cell_Type = "Malignant_cells_Liu")
Filbin_malignant <- Filbin_malignant %>% mutate(Cell_Type = "Malignant_cells_Filbin")

# Combine rows, allowing different genes to be retained
merged_data <- bind_rows(Liu_malignant, Filbin_malignant, .id = "Dataset")

# View output
head(merged_data)


# Ensure all columns (except 'Dataset' and 'Cell_Type') are numeric
numeric_data <- merged_data %>%
  select(-Dataset, -Cell_Type) %>%  # Exclude non-numeric columns
  mutate(across(everything(), as.numeric))  # Convert to numeric

# Compute the mean expression for each gene (ignoring NAs)
mean_expression <- numeric_data %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Transpose the result so genes are in rows
mean_expression <- as.data.frame(t(mean_expression))

# Add gene names as a column
mean_expression <- mean_expression %>% rownames_to_column(var = "Gene")

# View output
head(mean_expression)

DIPG_mean_exp <- mean_expression

########### Mean of tumour cell Z score in GBM

# Dataset names you define
Ebert <- Ebert_normalized_data_celltype
LeBlanc <- LeBlanc_normalized_data_celltype
Abdelfattah <- Abdelfattah_normalized_data_celltype
Chen <- Chen_normalized_data_celltype

library(dplyr)

# Subset Malignant_cells from both datasets
Ebert_malignant <- as.data.frame(Ebert["Malignant_cells", , drop = FALSE])
LeBlanc_malignant <- as.data.frame(LeBlanc["Malignant_cells", , drop = FALSE])
Abdelfattah_malignant <- as.data.frame(Abdelfattah["Malignant_cells", , drop = FALSE])
Chen_malignant <- as.data.frame(Chen["Malignant_cells", , drop = FALSE])

# Add row identifiers
Ebert_malignant <- Liu_malignant %>% mutate(Cell_Type = "Malignant_cells_Ebert")
LeBlanc_malignant <- LeBlanc_malignant %>% mutate(Cell_Type = "Malignant_cells_LeBlanc")
Abdelfattah_malignant <- Abdelfattah_malignant %>% mutate(Cell_Type = "Malignant_cells_Abdelfattah")
Chen_malignant <- Chen_malignant %>% mutate(Cell_Type = "Malignant_cells_Chen")

# Combine rows, allowing different genes to be retained
merged_data <- bind_rows(Ebert_malignant, LeBlanc_malignant, Abdelfattah_malignant, Chen_malignant, .id = "Dataset")

# View output
head(merged_data)

# Ensure all columns (except 'Dataset' and 'Cell_Type') are numeric
numeric_data <- merged_data %>%
  select(-Dataset, -Cell_Type) %>%  # Exclude non-numeric columns
  mutate(across(everything(), as.numeric))  # Convert to numeric

# Compute the mean expression for each gene (ignoring NAs)
mean_expression <- numeric_data %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Transpose the result so genes are in rows
mean_expression <- as.data.frame(t(mean_expression))

# Add gene names as a column
mean_expression <- mean_expression %>% rownames_to_column(var = "Gene")

# View output
head(mean_expression)

GBM_mean_exp <- mean_expression

# Plot DIPG vs GBM

library(ggplot2)

# Define the genes to highlight (including the '_name' suffix)
highlight_genes <- c("GABRA5_name", "GABRB2_name", "GABRA1_name", "GABRG2_name")

# Create a new column to flag selected genes
scatter_data$highlight <- ifelse(scatter_data$Gene %in% highlight_genes, "Highlighted", "Other")

# Create scatter plot
ggplot(scatter_data, aes(x = V1_DIPG, y = V1_GBM, color = highlight)) +
  geom_point(alpha = 0.7) +  # Scatter points with transparency
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", alpha = 0.3) +  # Linear regression with 95% CI
  geom_text(data = scatter_data %>% filter(Gene %in% highlight_genes), 
            aes(label = Gene), vjust = -1, size = 4, color = "red") +  # Add labels
  geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +  # Vertical line at x = 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +  # Horizontal line at y = 1
  scale_color_manual(values = c("Highlighted" = "red", "Other" = "black")) +  # Color mapping
  labs(x = "DIPG Mean Expression", 
       y = "GBM Mean Expression", 
       title = "DIPG vs GBM Gene Expression") +
  theme_minimal()

lm_model <- lm(V1_GBM ~ V1_DIPG, data = scatter_data)
summary_lm <- summary(lm_model)

# View R-squared and p-value
summary_lm$r.squared       # R-squared
summary_lm$coefficients    # Contains p-value for the slope

write.csv(scatter_data, "path to output folder/scatter_plot_DIPG_black_vs_GBM_blue_neuronal_reg_module_malignant_cells.csv", row.names=TRUE)

###### Heatmaps for Sup.Fig 2B - GABA gene expression in malignant cell states

# First Annotate malignant cell states (use code 'calculate_Neftel_states.Rmd by Ashwini Patil and save objects)

#load datasets and subset malignant cells (annotated with cellstate metadata)
#Ebert
obj_Ebert <- readRDS("path to file/Ebert.rds")
Idents(obj_Ebert) <- "combined_clusters_annot"
obj_Ebert_Mal <- subset(obj_Ebert, idents = c("Malignant_cells"))
Idents(obj_Ebert_Mal) <- "cellstate"
table(obj_Ebert_Mal$cellstate)
obj_Ebert_Mal <- subset(obj_Ebert_Mal, idents = c("AC", "MES1.Hypoxia.independent.", "MES2.Hypoxia.dependent.", "NPC1", "NPC2", "OPC"))

#LeBlanc
obj_LeBlanc <- readRDS("C:/Users/shardcl/OneDrive - University of South Australia/NHMRC Ideas 2021 GNT 2013180/Chloe_R_Scripts/LeBlanc/LeBlanc_tissue.rds")
Idents(obj_LeBlanc) <- "sample"
table(obj_LeBlanc$sample)
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1"))
Idents(obj_LeBlanc) <- "cell_type"
obj_LeBlanc_Mal <- subset(obj_LeBlanc, idents = c("malignant"))
Idents(obj_LeBlanc_Mal) <- "cellstate"
table(obj_LeBlanc_Mal$cellstate)
obj_LeBlanc_Mal <- subset(obj_LeBlanc_Mal, idents = c("AC", "MES1.Hypoxia.independent.", "MES2.Hypoxia.dependent.", "NPC1", "NPC2", "OPC"))

#Abdelfattah
obj_Abdelfattah <- readRDS("C:/Users/shardcl/OneDrive - University of South Australia/NHMRC Ideas 2021 GNT 2013180/Chloe_R_Scripts/Abdelfattah/Abdelfattah.rds")
Idents(obj_Abdelfattah) <- "Type"
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("GBM"))
Idents(obj_Abdelfattah) <- "Assignment"
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah, idents = c("Glioma"))
Idents(obj_Abdelfattah_Mal) <- "cellstate"
table(obj_Abdelfattah_Mal$cellstate)
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah_Mal, idents = c("AC", "MES1_HypoxIndep", "MES2_HypoxDep", "NPC1", "NPC2", "OPC"))

#Chen
obj_Chen <- readRDS("C:/Users/shardcl/OneDrive - University of South Australia/NHMRC Ideas 2021 GNT 2013180/Chloe_R_Scripts/Chen/Chen.rds")
Idents(obj_Chen) <- "orig.ident"
obj_Chen <- subset(obj_Chen, idents = c("PDC001", "PJ052", "PJ053", "PW016-703", "PW017-703", "PW032-710", "PW035-710", "PW039-705"))
Idents(obj_Chen) <- "cellkb"
obj_Chen_Mal <- subset(obj_Chen, idents = c("Tumor"))
Idents(obj_Chen_Mal) <- "cellstate"
table(obj_Chen_Mal$cellstate)
obj_Chen_Mal <- subset(obj_Chen_Mal, idents = c("AC", "MES1_HypoxIndep", "MES2_HypoxDep", "NPC1", "NPC2", "OPC"))

genes_of_interest <- c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", 
                       "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3", 
                       "GABRR1", "GABRR2", "GABRR3", "GABRD", "GABRQ", "GABRP", 
                       "GABRE")

# Map dataset-specific cell state labels to standardized labels
normalize_cellstate <- function(cellstate, dataset_name) {
  mapping <- list(
    "Ebert" = c("MES1.Hypoxia.independent." = "MES1", "MES2.Hypoxia.dependent." = "MES2"),
    "Adelfattah" = c("MES1_HypoxIndep" = "MES1", "MES2_HypoxDep" = "MES2"),
    "Chen" = c("MES1_HypoxIndep" = "MES1", "MES2_HypoxDep" = "MES2"),
    "LeBlanc" = c("MES1.Hypoxia.independent." = "MES1", "MES2.Hypoxia.dependent." = "MES2")
  )
  
  # Ensure the mapping works for all datasets
  standardized <- mapping[[dataset_name]][cellstate]
  return(ifelse(is.na(standardized), cellstate, standardized))
}

# Function to calculate average expression (mean), z-score, and percentage expression per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  # Normalize cellstate labels
  seurat_obj@meta.data$cellstate <- sapply(seurat_obj@meta.data$cellstate, normalize_cellstate, dataset_name)
  
  # Filter for genes present in the dataset
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  if (length(genes_present) == 0) stop(paste("No genes of interest found in", dataset_name))
  
  # Extract normalised expression matrix
  expr_mat <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  expr_df <- as.data.frame(t(expr_mat))  # genes as columns, cells as rows
  expr_df$cell_type <- seurat_obj@meta.data$cellstate
  
  # Calculate average expression per gene per cell type
  avg_expr <- expr_df %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), mean, na.rm = TRUE)) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression")
  
  # Min-max normalisation across cell types for each gene
  avg_expr <- avg_expr %>%
    group_by(gene) %>%
    mutate(minmax_expr = (avg_expression - min(avg_expression)) / (max(avg_expression) - min(avg_expression))) %>%
    ungroup()
  
  # Add dataset name
  avg_expr <- avg_expr %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_")) %>%
    select(gene, cell_type_dataset, minmax_expr)
  
  return(avg_expr)
}

plot_data_Ebert <- process_dataset(obj_Ebert_Mal, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc_Mal, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah_Mal, "Adelfattah")
plot_data_Chen <- process_dataset(obj_Chen_Mal, "Chen")

combined_plot_data <- bind_rows(
  plot_data_Ebert,
  plot_data_LeBlanc,
  plot_data_Adelfattah,
  plot_data_Chen
)

heatmap_matrix <- combined_plot_data %>%
  pivot_wider(names_from = cell_type_dataset, values_from = minmax_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Define your desired cell state order
cellstate_order <- c("NPC2", "NPC1", "OPC", "AC", "MES1", "MES2")

# Extract all unique datasets (e.g., Neftel, Ebert, etc.) from column names
dataset_labels <- unique(gsub(".*_", "", colnames(heatmap_matrix)))

# Build the ordered list of column names
ordered_columns <- unlist(lapply(cellstate_order, function(cs) {
  grep(paste0("^", cs, "_"), colnames(heatmap_matrix), value = TRUE)
}))

# Reorder the heatmap matrix
heatmap_matrix <- heatmap_matrix[, ordered_columns]
