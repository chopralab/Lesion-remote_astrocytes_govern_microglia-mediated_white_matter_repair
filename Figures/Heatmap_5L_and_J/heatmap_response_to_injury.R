# Load necessary libraries
# library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)


library(dplyr)
library(RColorBrewer)
# File path
# Read the csv file



create_heatmap_logFC_white_non_sig <- function(df, Title1, Title2, title, intensity_columns, limit, limit_lower, min_PValue = 0.01) {
  
  tryCatch({
    
    # Define intensity columns
    intensity_columns <- intensity_columns
    
    # Arrange dataframe by lipid column
    df <- df %>%
      arrange(type)
    
    title <- title
    labels_title_1 <- Title1
    labels_title_2 <- Title2
    Length1 <- 1
    Length2 <- 1
    
    annotation_labels <- c(rep(labels_title_1, each=Length1), rep(labels_title_2, each=Length2))
    
    # Update LogFC values based on PValue conditions
    df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
    df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
    df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
    # Function to save pheatmap as PDF
    save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
      stopifnot(!missing(x))
      stopifnot(!missing(filename))
      pdf(filename, width=width, height=height)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    
    # Setup conditions for FDR filtering
    conditions <- list(
      list(filter = TRUE, suffix = "_PVALUE_01"),
      list(filter = FALSE, suffix = "")
    )
    
    # Loop through conditions to generate and save heatmaps
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df %>% filter(min_PValue < 0.01)
        # df_filtered <- df %>% filter(Min_PValue_all < 0.05)
      } else {
        df_filtered <- df
      }
      df_filtered <- df_filtered %>% arrange(type,LogFC_WT,LogFC_cKO)#%>%(type,PValue_WT)
      # arrange(PValue_cKO)%>%
      # arrange(type)
      # df_filtered <- df_filtered %>% arrange(type, 
      #                                        ifelse(LogFC_cKO == 0, LogFC_WT, LogFC_cKO),
      #                                        ifelse(LogFC_WT == 0, LogFC_cKO, LogFC_WT))
      
      title_suffix <- condition$suffix
      heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
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
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$Row_Label,
                                        annotation_col = annotation_col_df)
      
      # Save the heatmap
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_newer/", title, title_suffix,"Limit_",limit, "_LogFC.pdf"), width = 11, height = 10)
      # save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, ".svg"), width = 7, height = 10)
    }
    
  }, error = function(e) {
    print(paste("Error processing file:Error message:", e))
  })
}
# create_heatmap_logFC_white_non_sig(df, "cKO", "WT", "Unique Lipids scaled sig arragned cKO PValue", intensity_columns,1,-1)


create_heatmap_logFC_Class_non_sig_white <- function(df, Title1, Title2, title, intensity_columns,limit,limit_lower, min_PValue = 0.05) {
  
  tryCatch({
    
    unique_types <- unique(df$type)
    # Update LogFC values based on PValue conditions
    df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
    df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
    df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
    for (single_type in unique_types) {
      
      # Filter dataframe by current type
      df_type <- df %>% 
        filter(type == single_type) %>%
        arrange(LogFC_difference)
      
      labels_title_1 <- Title1
      labels_title_2 <- Title2
      Length1 <- 1
      Length2 <- 1
      
      annotation_labels <- c(rep(labels_title_1, each=Length1), rep(labels_title_2, each=Length2))
      
      # Function to save pheatmap as PDF
      save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
      }
      
      # Setup conditions for FDR filtering
      conditions <- list(
        list(filter = TRUE, suffix = "_PVALUE_01"),
        list(filter = FALSE, suffix = "")
      )
      
      # Loop through conditions to generate and save heatmaps
      for (condition in conditions) {
        if (condition$filter) {
          df_filtered <- df %>% filter(min_PValue < 0.01)
          # df_filtered <- df %>% filter(Min_PValue_all < 0.05)
        } else {
          df_filtered <- df_type
        }
        df_filtered <- df_filtered %>% arrange(type,LogFC_WT,LogFC_cKO)
        title_suffix <- condition$suffix
        heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
        
        annotation_col_df <- data.frame(Labels = annotation_labels)
        rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
        
        heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                          main = paste0(title, title_suffix, " - ", single_type),
                                          cluster_rows = FALSE, 
                                          cluster_cols = FALSE, 
                                          show_colnames = TRUE, 
                                          show_rownames = TRUE,
                                          breaks = heatmap_breaks,
                                          border_color = "black", 
                                          fontsize = 10,
                                          fontsize_legend = 25,
                                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                          labels_row = df_filtered$lipid,
                                          annotation_col = annotation_col_df)
        
        # Save the heatmap
        save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_newer/", title, title_suffix, "_", single_type,"Limit_",limit, "_LogFC_limit.pdf"), width = 11, height = 10)
      }
    }
    
  }, error = function(e) {
    print(paste("Error processing file:Error message:", e$message))
  })
}
# 
# create_combined_heatmap_logFC <- function(df, col_titles, intensity_columns, title, limit, limit_lower, min_PValue = 0.05) {
#   
#   if (length(col_titles) != 3 || length(intensity_columns) != 3) {
#     stop("The length of col_titles and intensity_columns must be 3.")
#   }
#   
#   # Function to save pheatmap as PDF
#   save_pheatmap_pdf <- function(x, filename, width=6, height=8) {
#     stopifnot(!missing(x))
#     stopifnot(!missing(filename))
#     pdf(filename, width=width, height=height)
#     grid::grid.newpage()
#     grid::grid.draw(x$gtable)
#     dev.off()
#   }
#   
#   tryCatch({
#     
#     # Update LogFC values based on PValue conditions for each column
#     # for (i in 1:3) {
#     #   df[[intensity_columns[i]]] <- ifelse(df[[paste0("PValue_", col_titles[i])]] > min_PValue, 0, df[[intensity_columns[i]]])
#     # }
#     # Update LogFC values based on PValue conditions
#     df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
#     df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
#     df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
#     
#     
#     # Setup conditions for PValue filtering
#     conditions <- list(
#       list(filter = TRUE, suffix = "_PVALUE_05"),
#       list(filter = FALSE, suffix = "")
#     )
#     
#     # Generate the whole heatmap
#     for (condition in conditions) {
#       if (condition$filter) {
#         df_filtered <- df %>% filter(Min_PValue_all < 0.05)
#       } else {
#         df_filtered <- df
#       }
#       
#       df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])
#       
#       title_suffix <- condition$suffix
#       heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
#       
#       annotation_col_df <- data.frame(Labels = col_titles)
#       rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
#       
#       heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
#                                         main = paste0(title, title_suffix),
#                                         cluster_rows = FALSE, 
#                                         cluster_cols = FALSE, 
#                                         show_colnames = TRUE, 
#                                         show_rownames = TRUE,
#                                         breaks = heatmap_breaks,
#                                         border_color = "black", 
#                                         fontsize = 10,
#                                         fontsize_legend = 25,
#                                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
#                                         labels_row = df_filtered$Row_Label,
#                                         annotation_col = annotation_col_df)
#       
#       # Save the heatmap
#       save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix,"Limit_",limit, "_LogFC_all_3_Longer.pdf"), width = 11, height = 10)
#     }
#     
#     # Generate heatmaps for each unique 'type'
#     unique_types <- unique(df$type)
#     for (single_type in unique_types) {
#       
#       df_type <- df %>% filter(type == single_type)
#       
#       for (condition in conditions) {
#         if (condition$filter) {
#           df_filtered <- df_type %>% filter(Min_PValue_all < 0.05)
#         } else {
#           df_filtered <- df_type
#         }
#         
#         df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])
#         
#         title_suffix <- condition$suffix
#         heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
#         
#         annotation_col_df <- data.frame(Labels = col_titles)
#         rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
#         
#         heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
#                                           main = paste0(title, title_suffix, " - ", single_type),
#                                           cluster_rows = FALSE, 
#                                           cluster_cols = FALSE, 
#                                           show_colnames = TRUE, 
#                                           show_rownames = TRUE,
#                                           breaks = heatmap_breaks,
#                                           border_color = "black", 
#                                           fontsize = 10,
#                                           fontsize_legend = 25,
#                                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
#                                           labels_row = df_filtered$lipid,
#                                           annotation_col = annotation_col_df)
#         
#         # Save the heatmap
#         save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, "_", single_type,"Limit_",limit, "_LogFC_Unique_limit_all_3_Longer.pdf"), width = 11, height = 10)
#       }
#     }
#     
#   }, error = function(e) {
#     print(paste("Error processing file: Error message:", e$message))
#   })
# }
create_combined_heatmap_logFC <- function(df, col_titles, intensity_columns, title, limit, limit_lower, min_PValue = 0.05) {
  
  if (length(col_titles) != 3 || length(intensity_columns) != 3) {
    stop("The length of col_titles and intensity_columns must be 3.")
  }
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=6, height=35) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  tryCatch({
    
    # Update LogFC values based on PValue conditions
    df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
    df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
    df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
    
    # Setup conditions for PValue filtering
    conditions <- list(
      list(filter = TRUE, suffix = "_PVALUE_05"),
      list(filter = FALSE, suffix = "")
    )
    
    # Generate the whole heatmap
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df %>% filter(Min_PValue_all < 0.05)
      } else {
        df_filtered <- df
      }
      
      df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])
      first_occurrence <- !duplicated(df_filtered$type)
      df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
      
      title_suffix <- condition$suffix
      heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
      
      annotation_col_df <- data.frame(Labels = col_titles)
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
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$Row_Label,
                                        annotation_col = annotation_col_df)
      
      # Save the heatmap
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix,"Limit_",limit, "_LogFC_all_3_Longest_new_axis.pdf"), width = 11, height = 10)
    }
    
    # Generate heatmaps for each unique 'type'
    unique_types <- unique(df$type)
    for (single_type in unique_types) {
      
      df_type <- df %>% filter(type == single_type)
      
      for (condition in conditions) {
        if (condition$filter) {
          df_filtered <- df_type %>% filter(Min_PValue_all < 0.05)
        } else {
          df_filtered <- df_type
        }
        
        df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]])
        
        title_suffix <- condition$suffix
        heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
        
        annotation_col_df <- data.frame(Labels = col_titles)
        rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
        
        heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                          main = paste0(title, title_suffix, " - ", single_type),
                                          cluster_rows = FALSE, 
                                          cluster_cols = FALSE, 
                                          show_colnames = TRUE, 
                                          show_rownames = TRUE,
                                          breaks = heatmap_breaks,
                                          border_color = "black", 
                                          fontsize = 10,
                                          fontsize_legend = 25,
                                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                          labels_row = df_filtered$lipid,
                                          annotation_col = annotation_col_df)
        
        # Save the heatmap
        save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, "_", single_type,"Limit_",limit, "_LogFC_Unique_limit_all_3_Longest_new_axis.pdf"), width = 11, height = 10)
      }
    }
    
  }, error = function(e) {
    print(paste("Error processing file: Error message:", e$message))
  })
}














df <- read.csv("merged_updated.csv")

# df <- merged_lipid_data

names(df)
filtered_df <- df[(df$PValue_cKO < 0.01) | (df$PValue_WT < 0.01), ]
Title1 <-"WT"
Title2 <-"cKO"
title <- "All Sig Lipids PValue 01 WT LogFC arranged"
intensity_columns <- c("LogFC_WT", "LogFC_cKO")
create_heatmap_logFC_white_non_sig(filtered_df, Title1, Title2, title, intensity_columns,1,-1)
create_heatmap_logFC_Class_non_sig_white(filtered_df, Title1, Title2, title, intensity_columns,1,-1)
create_heatmap_logFC_white_non_sig(filtered_df, Title1, Title2, title, intensity_columns,2,-2)
create_heatmap_logFC_Class_non_sig_white(filtered_df, Title1, Title2, title, intensity_columns,2,-2)



df <- read.csv("merged_updated.csv")
names(df)
filtered_df <- df[(df$PValue_cKO < 0.05) | (df$PValue_WT < 0.05)| (df$PValue_injured_comparison < 0.05), ]
Title1 <-"WT"
Title2 <-"cKO"
Title3 <-"LogFC_Injured_Comparison All"
title <- "All Sig Lipids PValue 05 WT LogFC arranged"
intensity_columns <- c("LogFC_WT", "LogFC_cKO","LogFC_injured_comparison")
col_titles <- c(Title1, Title2,Title3)
create_combined_heatmap_logFC(filtered_df, col_titles, intensity_columns, title, 1, -1, min_PValue = 0.05)



df <- read.csv("merged_updated_with_healthy.csv")
names(df)
filtered_df <- df[(df$PValue_cKO < 0.05) | (df$PValue_WT < 0.05)| (df$PValue_injured_comparison < 0.05)| (df$Pvalue_healthy < 0.05), ]

write.csv(filtered_df, "filtered_df.csv", row.names = FALSE)


Title1 <-"WT"
Title2 <-"cKO"
Title3 <-"LogFC_Injured_Comparison All"
Title4 <-"LogFC_Healthy_comparison"
title <- "All Lipids All Comparisons"
intensity_columns <- c("LogFC_WT", "LogFC_cKO","LogFC_injured_comparison","LogFC_Healthy")
col_titles <- c(Title1, Title2,Title3,Title4)


create_combined_heatmap_logFC_all(filtered_df, col_titles, intensity_columns, title, 1, -1, min_PValue = 0.05)

create_combined_heatmap_logFC_all(df, col_titles, intensity_columns, title, 1, -1, min_PValue = 1)





create_combined_heatmap_logFC_all <- function(df, col_titles, intensity_columns, title, limit, limit_lower, min_PValue = 0.05) {
  
  # if (length(col_titles) < 3 || length(intensity_columns) < 3) {
  #   stop("The length of col_titles and intensity_columns must be 3.")
  # }
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width = 6, height = 35) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  tryCatch({
    
    # Update LogFC values based on PValue conditions
    df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
    df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
    df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
    df$LogFC_Healthy <- ifelse(df$Pvalue_healthy > min_PValue, 0, df$LogFC_Healthy)
    
    # Setup conditions for PValue filtering
    conditions <- list(
      list(filter = TRUE, suffix = "_PVALUE_05"),
      list(filter = FALSE, suffix = "")
    )
    
    # Generate the whole heatmap
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df %>% filter(Min_PValue_all < 0.05)
      } else {
        df_filtered <- df
      }
      
      # Define the custom order
      custom_order <- c("CE", "TAG", "PC", "PE", "SM", "Cer", "CAR", "PS", "PI", "PG", "FA")
      
      # Convert 'type' to a factor with the custom order
      df_filtered$type <- factor(df_filtered$type, levels = custom_order)
      
      # Arrange the dataframe
      # df_filtered <- df_filtered %>%
      #   arrange(type, 
      #           df_filtered[[intensity_columns[1]]], 
      #           df_filtered[[intensity_columns[2]]], 
      #           df_filtered[[intensity_columns[3]]], 
      #           df_filtered[[intensity_columns[4]]])
      df_filtered <- df_filtered %>%
        arrange(type, 
                df_filtered[[intensity_columns[1]]], 
                df_filtered[[intensity_columns[2]]])#, 
                # df_filtered[[intensity_columns[3]]], 
                # df_filtered[[intensity_columns[4]]])
      
      
      # Create a new column 'Row_Label' for labels in the heatmap
      first_occurrence <- !duplicated(df_filtered$type)
      df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
      
      title_suffix <- condition$suffix
      heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
      
      annotation_col_df <- data.frame(Labels = col_titles)
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
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$Row_Label,
                                        annotation_col = annotation_col_df)
      
      # Save the heatmap
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, "Limit_", limit, "_LogFC_All_Comparisons.pdf"), width = 11, height = 10)
    }
    
    # Generate heatmaps for each unique 'type'
    unique_types <- unique(df$type)
    for (single_type in unique_types) {
      
      df_type <- df %>% filter(type == single_type)
      
      for (condition in conditions) {
        if (condition$filter) {
          df_filtered <- df_type %>% filter(Min_PValue_all < 0.05)
        } else {
          df_filtered <- df_type
        }
        
        df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]], df_filtered[[intensity_columns[4]]])
        
        title_suffix <- condition$suffix
        heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
        
        annotation_col_df <- data.frame(Labels = col_titles)
        rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
        
        heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                          main = paste0(title, title_suffix, " - ", single_type),
                                          cluster_rows = FALSE, 
                                          cluster_cols = FALSE, 
                                          show_colnames = TRUE, 
                                          show_rownames = TRUE,
                                          breaks = heatmap_breaks,
                                          border_color = "black", 
                                          fontsize = 10,
                                          fontsize_legend = 25,
                                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                          labels_row = df_filtered$lipid,
                                          annotation_col = annotation_col_df)
        
        # Save the heatmap
        save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, "_", single_type, "Limit_", limit, "_LogFC_All_Comparisons.pdf"), width = 11, height = 10)
      }
    }
    
  }, error = function(e) {
    print(paste("Error processing file: Error message:", e$message))
  })
}











df <- read.csv("merged_updated_with_healthy.csv")
names(df)
filtered_df <- df[(df$PValue_injured_comparison < 0.05)| (df$Pvalue_healthy < 0.05), ]

# write.csv(filtered_df, "filtered_df.csv", row.names = FALSE)


Title1 <-"LogFC_Injured_Comparison"
Title2 <-"LogFC_Healthy_comparison"
# Title3 <-"LogFC_Injured_Comparison All"
# Title4 <-"LogFC_Healthy_comparison"
title <- "Comparison for 3extended"
# intensity_columns <- c("LogFC_WT", "LogFC_cKO","LogFC_injured_comparison","LogFC_Healthy")
intensity_columns <- c("LogFC_injured_comparison","LogFC_Healthy")

col_titles <- c(Title1, Title2)


create_combined_heatmap_logFC_all_updated(filtered_df, col_titles, intensity_columns, title, 1, -1, min_PValue = 0.05)

create_combined_heatmap_logFC_all_updated(df, col_titles, intensity_columns, title, 1, -1, min_PValue = 1)



create_combined_heatmap_logFC_all_updated <- function(df, col_titles, intensity_columns, title, limit, limit_lower, min_PValue = 0.05) {
  
  # if (length(col_titles) < 3 || length(intensity_columns) < 3) {
  #   stop("The length of col_titles and intensity_columns must be 3.")
  # }
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width = 6, height = 35) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  tryCatch({
    
    # Update LogFC values based on PValue conditions
    df$LogFC_cKO <- ifelse(df$PValue_cKO > min_PValue, 0, df$LogFC_cKO)
    df$LogFC_WT <- ifelse(df$PValue_WT > min_PValue, 0, df$LogFC_WT)
    df$LogFC_injured_comparison <- ifelse(df$PValue_injured_comparison > min_PValue, 0, df$LogFC_injured_comparison)
    df$LogFC_Healthy <- ifelse(df$Pvalue_healthy > min_PValue, 0, df$LogFC_Healthy)
    
    # Setup conditions for PValue filtering
    conditions <- list(
      list(filter = TRUE, suffix = "_PVALUE_05"),
      list(filter = FALSE, suffix = "")
    )
    
    # Generate the whole heatmap
    for (condition in conditions) {
      if (condition$filter) {
        df_filtered <- df %>% filter(Min_PValue_all < 0.05)
      } else {
        df_filtered <- df
      }
      
      # Define the custom order
      custom_order <- c("CE", "TAG", "PC", "PE", "SM", "Cer", "CAR", "PS", "PI", "PG", "FA")
      
      # Convert 'type' to a factor with the custom order
      df_filtered$type <- factor(df_filtered$type, levels = custom_order)
      
      # Arrange the dataframe
      # df_filtered <- df_filtered %>%
      #   arrange(type, 
      #           df_filtered[[intensity_columns[1]]], 
      #           df_filtered[[intensity_columns[2]]], 
      #           df_filtered[[intensity_columns[3]]], 
      #           df_filtered[[intensity_columns[4]]])
      df_filtered <- df_filtered %>%
        arrange(type, 
                df_filtered[[intensity_columns[1]]], 
                df_filtered[[intensity_columns[2]]])#, 
      # df_filtered[[intensity_columns[3]]], 
      # df_filtered[[intensity_columns[4]]])
      
      
      # Create a new column 'Row_Label' for labels in the heatmap
      # Create a new column 'Row_Label' for labels in the heatmap
      heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
      # 
      annotation_col_df <- data.frame(Labels = col_titles)
      rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
      # 
      # heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
      #                                   main = paste0(title, title_suffix),
      #                                   cluster_rows = FALSE, 
      #                                   cluster_cols = FALSE, 
      #                                   show_colnames = TRUE, 
      #                                   show_rownames = TRUE,
      #                                   breaks = heatmap_breaks,
      #                                   border_color = "black", 
      #                                   fontsize = 10,
      #                                   fontsize_legend = 25,
      #                                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      #                                   labels_row = df_filtered$Row_Label,
      #                                   annotation_col = annotation_col_df)
      # 
      
      first_occurrence <- !duplicated(df_filtered$type)
      
      df_filtered$Row_Label <- ifelse(first_occurrence, as.character(df_filtered$type), "")
      
      
      # annotation_col_df <- data.frame(Labels = annotation_labels)
      # rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
      
      heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                        main = paste0(title),
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = TRUE, 
                                        show_rownames = TRUE,
                                        breaks = heatmap_breaks,
                                        border_color = "black", 
                                        fontsize = 10,
                                        fontsize_legend = 25,
                                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                        labels_row = df_filtered$Row_Label,
                                        annotation_col = annotation_col_df)
      save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title,"Limit_",limit, "_LogFC_All_Comparisons.pdf"), width = 25, height = 35)
      
      
    }
    
    # Generate heatmaps for each unique 'type'
    unique_types <- unique(df$type)
    for (single_type in unique_types) {
      
      df_type <- df %>% filter(type == single_type)
      
      for (condition in conditions) {
        if (condition$filter) {
          df_filtered <- df_type %>% filter(Min_PValue_all < 0.05)
        } else {
          df_filtered <- df_type
        }
        
        df_filtered <- df_filtered %>% arrange(type, df_filtered[[intensity_columns[1]]], df_filtered[[intensity_columns[2]]], df_filtered[[intensity_columns[3]]], df_filtered[[intensity_columns[4]]])
        
        title_suffix <- condition$suffix
        heatmap_breaks <- seq(limit_lower, limit, length.out = 100)
        
        annotation_col_df <- data.frame(Labels = col_titles)
        rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
        
        heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                          main = paste0(title, title_suffix, " - ", single_type),
                                          cluster_rows = FALSE, 
                                          cluster_cols = FALSE, 
                                          show_colnames = TRUE, 
                                          show_rownames = TRUE,
                                          breaks = heatmap_breaks,
                                          border_color = "black", 
                                          fontsize = 10,
                                          fontsize_legend = 25,
                                          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                          labels_row = df_filtered$lipid,
                                          annotation_col = annotation_col_df)
        
        # Save the heatmap
        save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large_new/", title, title_suffix, "_", single_type, "Limit_", limit, "_LogFC_All_Comparisons.pdf"), width = 11, height = 10)
      }
    }
    
  }, error = function(e) {
    print(paste("Error processing file: Error message:", e$message))
  })
}








