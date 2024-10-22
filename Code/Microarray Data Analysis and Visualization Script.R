#################################################################
# Microarray Data Analysis and Visualization Script
# Author: Luis Fernando Nagano
# Date: 26 May 2022
# Description:
#   This script performs microarray data normalization,
#   differential expression analysis, and heatmap visualization
#   using R. It processes Affymetrix microarray data, identifies
#   differentially expressed genes, and visualizes the results.
#################################################################

#---------------------------
# 1. Preparation
#---------------------------

# Load required packages or install them if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("affy", "simpleaffy", "gcrma", "limma",
                       "ComplexHeatmap", "RColorBrewer", "circlize",
                       "colorRamps", "hgu95av2cdf")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set the working directory to the location of your data files
# Replace the path below with the path to your data directory
setwd("path/to/your/data/directory")

#---------------------------
# 2. Microarray Data Normalization
#---------------------------

# Read Affymetrix CEL files from the working directory
affy_data <- ReadAffy()
print(affy_data)

# Perform RMA normalization
rma_normalized <- rma(affy_data)
exprs_rma <- exprs(rma_normalized)
write.table(exprs_rma, "affymetrix_rma.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

# Perform MAS5 normalization
mas5_normalized <- mas5(affy_data)
exprs_mas5 <- exprs(mas5_normalized)
write.table(exprs_mas5, "affymetrix_mas5.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

# Perform GCRMA normalization
gcrma_normalized <- gcrma(affy_data)
exprs_gcrma <- exprs(gcrma_normalized)
write.table(exprs_gcrma, "affymetrix_gcrma.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

# Visualize the distribution of expression values using boxplots
# Prepare data frames for plotting
exprs_rma_df <- data.frame(exprs_rma)
exprs_mas5_df <- data.frame(log2(exprs_mas5))  # Log2 transformation for MAS5
exprs_gcrma_df <- data.frame(exprs_gcrma)

# Create boxplots for each normalization method
par(mfrow = c(1, 3))  # Set up the plotting area
boxplot(exprs_rma_df, main = "RMA Normalization",
        las = 2, cex.axis = 0.7, col = "lightblue")
boxplot(exprs_mas5_df, main = "MAS5 Normalization (Log2 Transformed)",
        las = 2, cex.axis = 0.7, col = "lightgreen")
boxplot(exprs_gcrma_df, main = "GCRMA Normalization",
        las = 2, cex.axis = 0.7, col = "lightcoral")
par(mfrow = c(1, 1))  # Reset plotting area

#---------------------------
# 3. Differential Expression Analysis
#---------------------------

# Read the targets file containing sample information
targets <- readTargets("targets.txt")
print(targets)

# Read CEL files using filenames from the targets file
abatch <- ReadAffy(filenames = targets$filename)

# Perform RMA normalization on the new data
eset <- rma(abatch)

# Extract phenotype information and create a design matrix
group_factor <- factor(targets$subgroup)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
print(design)

# Fit the linear model to the expression data
fit <- lmFit(eset, design)

# Define contrasts for group comparisons
contrast_matrix <- makeContrasts(
  Diff = "Poor_Prognosis - Good_Prognosis",
  levels = design
)
print(contrast_matrix)

# Apply contrasts and compute statistics with empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract all differentially expressed genes (DEGs)
all_DEGs <- topTable(fit2, coef = "Diff", number = Inf, adjust.method = "BH")
write.table(all_DEGs, file = "ALL_DEGs.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# Extract significant DEGs with adjusted p-value < 0.05 and logFC > 1
significant_DEGs <- topTable(fit2, coef = "Diff", number = Inf,
                             adjust.method = "BH", p.value = 0.05, lfc = 1)
write.table(significant_DEGs, file = "Significant_DEGs.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# Summarize the differential expression results
results <- decideTests(fit2)
summary(results)
vennDiagram(results)

#---------------------------
# 4. Heatmap Visualization
#---------------------------

# Extract expression data for significant DEGs
DEG_genes <- rownames(significant_DEGs)
DEG_exprs <- exprs(eset)[DEG_genes, ]

# Scale the expression data for better visualization
scale_rows <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, sd, na.rm = TRUE)
  scaled_x <- sweep(x, 1, row_means, "-")
  scaled_x <- sweep(scaled_x, 1, row_sds, "/")
  return(scaled_x)
}
scaled_DEG_exprs <- scale_rows(DEG_exprs)

# Prepare sample annotations for the heatmap
sample_conditions <- targets$subgroup
annotation_df <- data.frame(Group = sample_conditions)
rownames(annotation_df) <- colnames(scaled_DEG_exprs)
annotation_colors <- list(Group = c("Poor_Prognosis" = "red",
                                    "Good_Prognosis" = "blue"))

# Create a HeatmapAnnotation object
heatmap_annotation <- HeatmapAnnotation(df = annotation_df,
                                        col = annotation_colors)

# Generate the heatmap
heatmap <- Heatmap(scaled_DEG_exprs,
                   name = "Expression Z-score",
                   top_annotation = heatmap_annotation,
                   show_row_names = FALSE,
                   show_column_names = TRUE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = "pearson",
                   clustering_method_rows = "ward.D2",
                   clustering_distance_columns = "pearson",
                   clustering_method_columns = "ward.D2",
                   col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Save the heatmap as a PNG file
png(filename = "DEG_heatmap.png", width = 10, height = 8, units = "in", res = 300)
draw(heatmap)
dev.off()

# Alternatively, display the heatmap in the R plotting window
draw(heatmap)
