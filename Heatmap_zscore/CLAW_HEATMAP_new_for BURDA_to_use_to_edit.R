library(tidyverse)
library(pheatmap)
library(RColorBrewer)


library(dplyr)
library(RColorBrewer)
# File path

library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(grid)

create_heatmap_Z_score <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  
  
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()

  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_01_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "")
  )
  
  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      df_filtered <- df %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      df_filtered <- df %>% filter(PValue < condition$threshold)
    } else {
      df_filtered <- df
    }
    
    if(nrow(df_filtered) < 1) {
      next
    }
    
    df_filtered <- df_filtered %>% arrange(type)
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-2, 2, length.out = 100)
    annotation_data <- data.frame(type = df_filtered$type)
    first_occurrence <- !duplicated(df_filtered$type)
    df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
    heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                      main = paste0(title, title_suffix),
                                      cluster_rows = FALSE, 
                                      cluster_cols = FALSE, 
                                      show_colnames = TRUE, 
                                      show_rownames = TRUE,
                                      breaks = heatmap_breaks,
                                      border_color = "black", 
                                      fontsize = 10,
                                      fontsize_legend = 25,
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large/", title, title_suffix, "Zscore_2_limit.pdf"), width = 11, height = 10)
  }
}

create_heatmap_by_class_adjusted <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe by lipid
  # df <- df %>%
  #   arrange(lipid)
  df <- df %>%
    arrange(logFC)
  # Read the names to be replaced from the Excel file
  names_to_fix <- readxl::read_xlsx("Names_to_fix.xlsx")
  
  names_to_fix <- names_to_fix[!is.na(names_to_fix$New_name) & names_to_fix$New_name != "",]
  name_changes <- setNames(names_to_fix$New_name, names_to_fix$lipid)
  
  # Replace lipid names in df using the named vector
  df$lipid <- ifelse(df$lipid %in% names(name_changes), name_changes[df$lipid], df$lipid)
  
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Get the blank column (assuming it's the last column in the dataframe)
  blank_col <- tail(names(df), 1)
  
  # Subtract the blank column from the intensity columns and set negative values to 0
  df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
  
  df <- df %>%
    mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
  
  # Remove rows where all values in intensity_columns are zeros
  df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
  
  ## Inserted code to divide each value in intensity_columns by the sum of that column
  df <- df %>% 
    mutate(across(all_of(intensity_columns), 
                  ~ . / sum(. , na.rm = TRUE)))
  
  # Compute Z-scores
  df <- df %>%
    rowwise() %>%
    mutate(across(all_of(intensity_columns),
                  ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
                    sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
    ungroup()
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  conditions <- list(
    list(filter = TRUE, suffix = "_PVALUE_01_by_class"),
    list(filter = FALSE, suffix = "_by_class")
  )
  
  for (lipid_type in unique(df$type)) {
    df_filtered_by_type <- df %>% filter(type == lipid_type)
    
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df_filtered_by_type %>% filter(PValue < 0.01)
      } else {
        df_filtered <- df_filtered_by_type
      }
      if(nrow(df_filtered) < 1) {
        next
      }
      
      title_suffix <- condition$suffix
      heatmap_breaks <- seq(-2, 2, length.out = 100)
      annotation_data <- data.frame(type = df_filtered$type)
      
      # Adjust the row annotations to use lipid:
      annotation_col_df <- data.frame(Labels = annotation_labels)
      rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
      
      heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                        main = paste0(title, "_", lipid_type, title_suffix),
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = TRUE, 
                                        show_rownames = TRUE,
                                        breaks = heatmap_breaks,
                                        border_color = "black", 
                                        fontsize = 10,
                                        fontsize_legend = 25,
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$lipid, # Use lipid for row labels
                                        annotation_col = annotation_col_df)  # Column annotations
      
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_by_class/", title, "_", lipid_type, title_suffix, ".pdf"), width = 11, height = 10)
    }
  }
}


# Check if heatmaps folder exists, if not, create it
if (!dir.exists("heatmaps_large")) {
  dir.create("heatmaps_large")
}
if (!dir.exists("heatmaps_large_by_class")) {
  dir.create("heatmaps_large_by_class")
}
# Get the list of files with "_full.csv"
files <- list.files(path = "results", pattern = "_full.csv$", full.names = TRUE)

# Create heatmap for each file
# lapply(files, create_heatmap_with_blank)

# Create heatmap for each file
lapply(files, create_heatmap_Z_score)



# lapply(files, create_heatmap_by_class)
lapply(files, create_heatmap_by_class_adjusted)


# create_heatmap_by_class <- function(file_path) {
#   
#   # Read and preprocess the data
#   df <- read.csv(file_path, stringsAsFactors = FALSE)
#   
#   # Arrange dataframe by lipid
#   df <- df %>%
#     arrange(lipid)
#   
#   # Identify intensity columns
#   intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
#   
#   # Get the blank column (assuming it's the last column in the dataframe)
#   blank_col <- tail(names(df), 1)
#   
#   # Subtract the blank column from the intensity columns and set negative values to 0
#   df[intensity_columns] <- pmax(df[intensity_columns] - df[[blank_col]], 0)
#   # Replace infinite values with NA for data frames
#   df <- df %>%
#     mutate(across(all_of(intensity_columns), ~replace(., is.infinite(.), NA)))
#   
#   # Replace NAs with 0 in intensity columns
#   df <- df %>%
#     mutate(across(all_of(intensity_columns), ~replace(., is.na(.), 0)))
#   
#   # Remove rows where all values in intensity_columns are zeros
#   df <- df %>% filter(rowSums(.[intensity_columns]) > 0)
#   
#   # Apply transformations
#   # df <- df %>%
#   #   mutate(across(all_of(intensity_columns), log2)) %>%
#   #   rowwise() %>%
#   #   mutate(mean_intensity = mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) %>%
#   #   mutate(across(all_of(intensity_columns), ~ . - mean_intensity)) %>%
#   #   select(-mean_intensity, -tail(names(df), 1))
#   
#   # Compute Z-scores
#   df <- df %>%
#     rowwise() %>%
#     mutate(across(all_of(intensity_columns),
#                   ~ ( . - mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) /
#                     sd(c_across(all_of(intensity_columns)), na.rm = TRUE))) %>%
#     ungroup()
#   
#   # Sort the dataframe by lipid in alphabetical order
#   df <- df %>%
#     arrange(lipid)
#   
#   title <- sub("_full.csv", "", basename(file_path))
#   labels_title_1 <- unique(df$Title_1)
#   labels_title_2 <- unique(df$Title_2)
#   annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
#   
#   # Function to save pheatmap as PDF
#   save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
#     stopifnot(!missing(x))
#     stopifnot(!missing(filename))
#     pdf(filename, width=width, height=height)
#     grid::grid.newpage()
#     grid::grid.draw(x$gtable)
#     dev.off()
#   }
#   
#   # Setup conditions for FDR filtering
#   conditions <- list(
#     list(filter = TRUE, suffix = "_PVALUE_05_by_class"),
#     list(filter = FALSE, suffix = "_by_class")
#   )
#   
#   # Loop through each unique 'type'
#   for (lipid_type in unique(df$type)) {
#     df_filtered_by_type <- df %>% filter(type == lipid_type)
#     
#     # Loop through conditions to generate and save heatmaps
#     for (condition in conditions) {
#       if (condition$filter) {
#         df_filtered <- df_filtered_by_type %>% filter(PValue < 0.05)
#       } else {
#         df_filtered <- df_filtered_by_type
#       }
#       if(nrow(df_filtered) < 1) {
#         next
#       }
#       
#       title_suffix <- condition$suffix
#       heatmap_breaks <- seq(-1, 1, length.out = 100)
#       annotation_data <- data.frame(type = df_filtered$type)
#       
#       # Using lipid as row label
#       rownames(df_filtered) <- df_filtered$lipid
#       
#       annotation_col_df <- data.frame(Labels = annotation_labels)
#       rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
#       
#       heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
#                                         main = paste0(title, "_", lipid_type, title_suffix),
#                                         cluster_rows = FALSE, 
#                                         cluster_cols = FALSE, 
#                                         show_colnames = TRUE, 
#                                         show_rownames = TRUE,
#                                         breaks = heatmap_breaks,
#                                         border_color = "black", 
#                                         fontsize = 10,
#                                         fontsize_legend = 25,
#                                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
#       
#       # Save the heatmap
#       save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_by_class/", title, "_", lipid_type, title_suffix, ".pdf"), width = 11, height = 10)
#     }
#   }
# }


