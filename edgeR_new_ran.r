library(forcats)
library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(tidyverse)
library(factoextra)
library(grid)
library(ggsignif)
library(stringr)
library(reshape2)
library(gtools)
library(tidyr)
library(ggbeeswarm) #https://github.com/eclarke/ggbeeswarm
library(ggrepel)
library(scales)
library(svglite)
library(ggplot2)
#library(ggExtra) #https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html

library(cowplot)
library(ggridges)

# library(limma)
# library(writexl)


# install.packages("ggforce")
library(ggforce)
getwd()

# read the variable from the text file
# cwd <- readLines("Variable_Storage/folder_path.txt")[1]
# cwd
# setwd(cwd)



# setwd(cwd)




file_list = list.files(path="Pre_EdgeR", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)

# file_list <- file_list[grepl("^(Original_scaled|Standard_subtracted_scaled|Normalized_scaled)", file_list)]

jj <-file_list[1]
file_list[2]
for (jj in file_list){
  excel_file <- read_csv(paste0("Pre_EdgeR/",jj,sep=""))
  dir.create("plots", F)
  dir.create("results", F)
  dir.create("plots/count", F)
  
  title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub(":", "_",excel_file$Title1[1])
  Title2 <- gsub(":", "_",excel_file$Title2[1])
  
  # Remove "__", " ", and "-" from title_for_plot
  title_for_plot <- gsub("__| |-", "", title_for_plot)
  
  # Remove "__", " ", and "-" from Title1
  Title1 <- gsub("__| |-", "", Title1)
  
  # Remove "__", " ", and "-" from Title2
  Title2 <- gsub("__| |-", "", Title2)
  
  blank_name <- excel_file$Blank_name[1]
  length1 <- excel_file$length1[1]
  length2 <- excel_file$length2[1]
  
  
  
  
  excel_file <- select(excel_file, -c(Title1, Title2, Title, length1, length2,Blank_name, Class2,Class3))
  
  ##Change Class to Type and change Lipid to lipid
  # excel_file <- excel_file %>%
  #   rename(lipid = `Lipid`) %>%
  #   rename(type= `Class`)
  excel_file <- excel_file %>%
    rename(lipid = `Lipid`, type = `Class`) %>%
    filter(!(type == "CAR" & grepl("_QUAL$", lipid))) %>%
    filter(!grepl("STD", type)) %>% filter(!grepl("STD", lipid))# Remove rows where type contains "STD"
  
  # library(dplyr)
  # 
  # excel_file <- excel_file %>%
  #   rename(
  #     lipid = Lipid,
  #     type = Class
  #   )
  
  if (isTRUE(!is.na(length1)) && isTRUE(!is.na(length2))) {
    if (isTRUE(length1 == 1) || isTRUE(length2 == 1)) {
      next
    }
  }
  
  
  
  gr1 = c(rep("GR1",length1))
  gr2 = c(rep("GR2",length2))
  
  
  groups_PCA = c((rep(Title1,length1)),(rep(Title2,length2)))
  
  groups_PCA
  
  
  
  
  blank_name
  
  
  
  
  
  cells_lipid_expr <- excel_file
  
  # Create a dataframe with Length1 and Length2 columns
  df2_other <- data.frame(Length1 = rep(length1, dim(cells_lipid_expr)[1]), 
                          Length2 = rep(length2, dim(cells_lipid_expr)[1]),
                          Title_1 = rep(Title1, dim(cells_lipid_expr)[1]),
                          Title_2 = rep(Title2, dim(cells_lipid_expr)[1]))
  
  #EdgeR groups
  gr_expr = c(gr1,gr2,
              blank_name) %>%
    factor(levels = c(blank_name, "GR1", "GR2"))
  design_expr = model.matrix(~gr_expr)
  
  contrasts_expr = makeContrasts(
    H = gr_exprGR1 - gr_exprGR2,
    levels = design_expr
  )
  
  
  
  ###Functions to do edgeR analysis
  
  perform_analysis_raw <- function(counts, design_mat, gr) {
    
    data.edgeR <- DGEList(counts = counts %>%
                            na.omit %>%
                            mutate(lipid = make.unique(lipid)) %>%
                            select( -type) %>%
                            # select(-Transition) %>%
                            column_to_rownames("lipid"),
                          group = gr
    )
    
    data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
    data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
    data.edgeR
  }
  
  calculate_significance <- function(dge, contrast) {
    dge %>%
      glmFit() %>%
      glmLRT(contrast = contrast)
  }
  
  experiment_helper <- function(df) {
    df %>%
      perform_analysis_raw(design_expr, gr_expr) %>%
      calculate_significance(contrasts_expr) %>%
      topTags(1500000) %>%
      as.data.frame() %>%
      rownames_to_column("lipid") %>%
      as_tibble()
    
  }
  
  cells_lipid_expr <- cells_lipid_expr %>%
    mutate(across(where(is.numeric), ~round(.)))
  cl_e1_tbl <-
    cells_lipid_expr %>%
    experiment_helper
  cl_e1_tbl
  
  names(cells_lipid_expr)
  
  
  
  xx <- merge(cells_lipid_expr,cl_e1_tbl)
  names(cl_e1_tbl)
  names(xx)
  cl_e1_tbl
  # lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG")
  # lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a")
  # 
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  # Define the filtering conditions
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "All_Lipids")
  )
  
  # Loop through each condition
  for (condition in conditions) {
    # Apply the appropriate filtering based on the condition
    if (condition$filter == "FDR") {
      data_to_plot <- xx %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      data_to_plot <- xx %>% filter(PValue < condition$threshold)
    } else {
      data_to_plot <- xx
    }
    # Check if there are at least 3 lipids after filtering
    if (nrow(data_to_plot) < 3) {
      next
    }
    
    # Calculate counts per type for ridge plot labels
    type_counts <- data_to_plot %>%
      group_by(type) %>%
      summarise(count = n(), .groups = 'drop')
    
    # Create labels with counts (just numbers in parentheses, no "n=")
    type_labels <- setNames(paste0(type_counts$type, " (", type_counts$count, ")"), type_counts$type)
    
    # Generate the ridge plot
    ridge_plot <- data_to_plot %>%
      ggplot(aes(x = logFC, y = type, fill = type)) +
      geom_density_ridges2(alpha = 0.5, size = .5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_classic() +
      ggtitle(title_for_plot) +
      xlab(paste0("Fold change, lipids ", ifelse(condition$filter == "NONE", "all", paste0(condition$filter, "<", condition$threshold)))) +
      ylab("") +
      scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
      scale_y_discrete(limits = rev, labels = function(x) type_labels[x]) +
      theme(legend.position = "right")
    
    # Save the ridge plot to both PDF and SVG formats
    ggsave(paste0("plots/Ridge_Plot", condition$suffix, "_", title_for_plot, ".pdf"))
    ggsave(paste0("plots/Ridge_Plot", condition$suffix, "_", title_for_plot, ".svg"))
  }
  
  # Create a directory for saving pie charts
  # if (!dir.exists("plots/count")) {
  #   dir.create("plots/count", recursive = TRUE)
  # }
  
  # Loop through each condition
  for (condition in conditions) {
    # Apply the appropriate filtering based on the condition
    if (condition$filter == "FDR") {
      data_to_plot <- xx %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      data_to_plot <- xx %>% filter(PValue < condition$threshold)
    } else {
      data_to_plot <- xx
    }
    
    # Count the number of significant lipids per class
    lipid_counts <- data_to_plot %>%
      group_by(type) %>%
      summarise(count = n()) %>%
      ungroup()
    
    # Check if there are at least 3 lipids after filtering
    if (nrow(lipid_counts) < 3) {
      next
    }
    
    # Generate the pie chart
    pie_chart <- ggplot(lipid_counts, aes(x = "", y = count, fill = type)) +
      geom_bar(stat = "identity", width = 1, alpha = 0.5) +
      coord_polar(theta = "y") +
      theme_void() +
      ggtitle(title_for_plot) +
      scale_fill_manual(values = lipid_class_colors, name = "Lipid class") +
      theme(legend.position = "right")
    
    # Save the pie chart to PDF and SVG formats
    ggsave(paste0("plots/count/Pie_Chart", condition$suffix, "_", title_for_plot, ".pdf"), plot = pie_chart)
    ggsave(paste0("plots/count/Pie_Chart", condition$suffix, "_", title_for_plot, ".svg"), plot = pie_chart)
  }
  
  # make_volcano_plot <- function(df, title) {
  #   df %>%
  #     mutate(sig = factor(FDR < 0.10)) %>%
  #     ggplot(aes(logFC, -log10(FDR), color = sig)) +
  #     geom_point() +
  #     scale_color_manual(values = c("none" = "black", "TRUE" = "red")) +
  #     guides(color = F) +
  #     ggtitle(title)
  # }
  # 
  # cl_e1_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot,sep=''))
  # ggsave(paste("plots/Volcano_Plot_",title_for_plot,".png",sep=''))
  # 
  # 
  
  
  write_summary_and_results <- function(tbl, df,df2, name) {
    df<-cbind(df2, df)
    tbl %>%
      # merge(df2) %>%
      merge(df) %>%
      as_tibble() %>%
      arrange(FDR) -> results
    
    write_csv(results, paste0("results/", name, "_full.csv"))
    
    results %>%
      group_by(type) %>%
      summarise(Down_FDR = sum(logFC < 0 & FDR < 0.1),
                Up_FDR = sum(logFC > 0 & FDR < 0.1),Down_Pvalue = sum(logFC < 0 & PValue < 0.05),
                Up_Pvalue = sum(logFC > 0 & PValue < 0.05),Down_Pvalue_01 = sum(logFC < 0 & PValue < 0.01),
                Up_Pvalue_01 = sum(logFC > 0 & PValue < 0.01)) %>%
      write_csv(paste0("results/", name, "_summary.csv"))
  }
  
  dir.create("results", F)
  
  cl_e1_tbl %>% write_summary_and_results(cells_lipid_expr,df2_other, title_for_plot)
  # filtered_xxx <- xx %>%
  #   filter(FDR < 0.1)
  # 
  # if (isTRUE(nrow(filtered_xxx) < 2)) {
  #   next
  # }
  
  cl_e1_tbl
  
  # Define all_groups upfront
  all_groups <- select(excel_file, -c(blank_name, lipid, type))
  all_groups <- names(all_groups)
  
  # Setup conditions for FDR and PValue filtering
  conditions <- list(
    list(filter = "FDR", threshold = 0.1, suffix = "_FDR_"),
    list(filter = "PValue", threshold = 0.05, suffix = "_PVALUE_05_"),
    list(filter = "PValue", threshold = 0.01, suffix = "_PVALUE_01_"),
    list(filter = "NONE", threshold = NULL, suffix = "All_Lipids")
  )
  
  # Blank subtraction conditions
  blank_conditions <- list(
    list(subtract = FALSE, suffix = ""),
    list(subtract = TRUE, suffix = "_BLANK_SUBTRACTED")
  )
  
  # Loop through filter conditions
  for (condition in conditions) {
    
    # Filter data based on conditions
    if (condition$filter == "FDR") {
      filtered_xx <- xx %>% filter(FDR < condition$threshold)
    } else if (condition$filter == "PValue") {
      filtered_xx <- xx %>% filter(PValue < condition$threshold)
    } else {
      filtered_xx <- xx
    }
    if (nrow(filtered_xx) < 3) {
      next
    }
    # Loop through blank subtraction conditions
    for (blank_condition in blank_conditions) {
      
      if (blank_condition$subtract) {
        # Subtract the blank column and set values below 0 to 0
        filtered_xx[all_groups] <- pmax(filtered_xx[all_groups] - filtered_xx[[blank_name]], 0)
      }
      # Replace infinite values with NA for data frames
      filtered_xx <- filtered_xx %>%
        mutate(across(all_of(all_groups), ~replace(., is.infinite(.), NA)))
      
      
      
      
      filtered_xx <- filtered_xx %>%
        mutate(across(all_of(all_groups), ~replace(., is.na(.), 0)))
      
      # Remove rows where all values in all_groups are zeros
      filtered_xx <- filtered_xx %>% filter(rowSums(.[all_groups]) > 0)
      
      # Rest of your PCA code
      cells_lipid_expr_subset <- filtered_xx[, all_groups]
      cells_lipid_expr_transposed <- t(cells_lipid_expr_subset)
      
      # Remove columns with zero variance
      # cells_lipid_expr_transposed <- cells_lipid_expr_transposed[, apply(cells_lipid_expr_transposed, 2, var) != 0]
      # 
      # # Check if you have at least 2 columns left after removing zero variance columns
      # if (ncol(cells_lipid_expr_transposed) < 2) {
      #   next
      # }
      
      
      
      pca_result <- prcomp(cells_lipid_expr_transposed, center = TRUE, scale = TRUE)
      
      group_data <- data.frame(
        sample = all_groups,
        group = c(rep(Title1, length1), rep(Title2,length2))
      )
      
      pca_scores <- as.data.frame(pca_result$x[, 1:2])
      
      plot_data <- data.frame(
        PC1 = pca_scores$PC1,
        PC2 = pca_scores$PC2,
        sample = rownames(pca_scores),
        group = group_data$group
      )
      
      pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
        geom_point(size = 3) +
        geom_text_repel(size = 3) +
        stat_ellipse(aes(fill = group), level = 0.95, geom = "polygon", alpha = 0.2) +
        theme_minimal() +
        labs(color = "Groups") +
        ggtitle(title_for_plot) +
        scale_color_manual(values = c("red", "black")) +
        scale_fill_manual(values = c("red", "black")) +
        xlab(paste0("PC1: ", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 2), "% variance")) +
        ylab(paste0("PC2: ", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 2), "% variance"))
      
      
      # Adjust the title of the saved plot to incorporate both the filter and blank subtraction conditions
      ggsave(paste("plots/PCA_", title_for_plot, condition$suffix, blank_condition$suffix, ".pdf", sep = ''), width = 12, height = 10)
    }
  }
  
  
  
  
}

