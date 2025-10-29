
# Loading required libraries
library(ggplot2)
library(ggrepel) # for better label positioning
library(reshape2) # for melting dataframes
library(FactoMineR) # for PCA


labels <- read.csv("labels.csv", stringsAsFactors = FALSE)
df<- read.csv("Blank_subtracted_and_column_normalized.csv", stringsAsFactors = FALSE)

# Transposing the df dataframe as you want dots per column
df_t <- t(df)

# Performing PCA on the transposed data
res.pca <- PCA(df_t, graph = FALSE)

# Extracting the PCA results
pca_data <- as.data.frame(res.pca$ind$coord)

# Merging the PCA results with the labels dataframe
pca_data$Sample.Name <- rownames(pca_data)
final_data <- merge(pca_data, labels, by = "Sample.Name")

# Defining the colors and shapes based on your specifications
shape_mapping <- c("28d ISCI" = 17, "healthy" = 16) # 17 is triangle, 16 is circle in ggplot
color_mapping <- c("cKO" = rgb(167/255, 62/255, 186/255), "WT" = "black")

# Plotting the PCA
ggplot(final_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(shape = `time.point`, color = Genotype), size = 3) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  ggrepel::geom_label_repel(aes(label = Sample.Name), box.padding = 0.5) +
  stat_ellipse(aes(group = interaction(`time.point`, Genotype), shape = `time.point`, color = Genotype), type = "norm", level = 0.68) + # adding ellipses
  labs(
    title = "PCA Plot",
    x = paste("PC1: ", round(res.pca$eig[1, 2], 2), "% variance"),
    y = paste("PC2: ", round(res.pca$eig[2, 2], 2), "% variance")
  ) +
  theme_minimal()

ggsave("Blank_subtracted_and_column_normalized.svg")
# Assuming you have your dataframes loaded as 'df' and 'labels'
# Also assuming you have required libraries installed



####No eclipse
# Transposing the df dataframe as you want dots per column
df_t <- t(df)

# Performing PCA on the transposed data
res.pca <- PCA(df_t, graph = FALSE)

# Extracting the PCA results
pca_data <- as.data.frame(res.pca$ind$coord)

# Merging the PCA results with the labels dataframe
pca_data$Sample.Name <- rownames(pca_data)
final_data <- merge(pca_data, labels, by = "Sample.Name")

# Defining the colors and shapes based on your specifications
shape_mapping <- c("28d ISCI" = 17, "healthy" = 16) # 2 is triangle, 1 is circle in ggplot
color_mapping <- c("cKO" = rgb(167/255, 62/255, 186/255), "WT" = "black")

# Plotting the PCA
ggplot(final_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(shape = `time.point`, color = Genotype), size = 3) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  ggrepel::geom_label_repel(aes(label = Sample.Name), box.padding = 0.5) +
  labs(
    title = "PCA Plot",
    x = paste("PC1: ", round(res.pca$eig[1, 2], 2), "% variance"),
    y = paste("PC2: ", round(res.pca$eig[2, 2], 2), "% variance")
  ) +
  theme_minimal()

