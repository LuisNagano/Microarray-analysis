# Advanced Microarray Data Analysis and Visualization Pipeline in R

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Microarray%20Analysis-brightgreen)](https://www.bioconductor.org/)

*Unlocking Insights from High-Throughput Gene Expression Data*

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installing Required Packages](#installing-required-packages)
- [Usage](#usage)
  - [1. Preparation](#1-preparation)
  - [2. Microarray Data Normalization](#2-microarray-data-normalization)
  - [3. Differential Expression Analysis](#3-differential-expression-analysis)
  - [4. Heatmap Visualization](#4-heatmap-visualization)
- [Results](#results)
- [Visualization](#visualization)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

---

## Introduction

High-throughput microarray technologies have revolutionized genomic research, enabling simultaneous analysis of thousands of genes. This R pipeline provides a comprehensive solution for microarray data normalization, differential expression analysis, and advanced visualization, specifically tailored for Affymetrix microarray platforms.

By leveraging robust statistical methods and high-quality graphics, this pipeline facilitates the identification of differentially expressed genes (DEGs) and generates insightful visualizations to aid researchers and data scientists in interpreting complex gene expression data.

---

## Features

- **Comprehensive Normalization Methods**: Supports RMA, MAS5, and GCRMA normalization techniques to preprocess raw Affymetrix CEL files.
- **Robust Differential Expression Analysis**: Utilizes the `limma` package for reliable identification of DEGs between experimental conditions.
- **Advanced Visualization**: Generates publication-ready heatmaps using `ComplexHeatmap`, with customizable annotations and clustering.
- **Scalability**: Efficiently handles large datasets, making it suitable for high-throughput analyses.
- **Modular and Extensible**: Designed with modular code structure for easy customization and integration with additional analyses like pathway enrichment.

---

## Installation

### Prerequisites

- **R**: Version 4.0.0 or higher.
- **RStudio** (optional but recommended): For an enhanced coding environment.

### Installing Required Packages

Ensure you have an active internet connection for package installation.

```R
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# List of required packages
required_packages <- c("affy", "simpleaffy", "gcrma", "limma",
                       "ComplexHeatmap", "RColorBrewer", "circlize",
                       "colorRamps", "hgu95av2cdf")

# Install and load packages
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg, ask = FALSE)
        library(pkg, character.only = TRUE)
    }
}
```

---

## Usage

### 1. Preparation

- **Set Working Directory**: Replace `"path/to/your/data/directory"` with the path to your data directory containing the CEL files and `targets.txt`.

```R
setwd("path/to/your/data/directory")
```

- **Data Requirements**:
  - **Affymetrix CEL Files**: Raw microarray data files.
  - **Targets File**: A tab-delimited text file named `targets.txt` containing sample information.

**Example `targets.txt`:**

```
filename       subgroup
Sample1.CEL    Poor_Prognosis
Sample2.CEL    Good_Prognosis
Sample3.CEL    Poor_Prognosis
Sample4.CEL    Good_Prognosis
```

### 2. Microarray Data Normalization

The script reads Affymetrix CEL files and performs normalization using RMA, MAS5, and GCRMA methods.

```R
# Read Affymetrix CEL files
affy_data <- ReadAffy()
print(affy_data)

# RMA Normalization
rma_normalized <- rma(affy_data)
exprs_rma <- exprs(rma_normalized)
write.table(exprs_rma, "affymetrix_rma.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

# MAS5 Normalization
mas5_normalized <- mas5(affy_data)
exprs_mas5 <- exprs(mas5_normalized)
write.table(exprs_mas5, "affymetrix_mas5.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

# GCRMA Normalization
gcrma_normalized <- gcrma(affy_data)
exprs_gcrma <- exprs(gcrma_normalized)
write.table(exprs_gcrma, "affymetrix_gcrma.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)
```

- **Boxplot Visualization**: Visualize the distribution of expression values after normalization.

```R
# Prepare data frames
exprs_rma_df <- data.frame(exprs_rma)
exprs_mas5_df <- data.frame(log2(exprs_mas5))  # Log2 transformation
exprs_gcrma_df <- data.frame(exprs_gcrma)

# Plot boxplots
par(mfrow = c(1, 3))  # Layout for multiple plots
boxplot(exprs_rma_df, main = "RMA Normalization",
        las = 2, cex.axis = 0.7, col = "lightblue")
boxplot(exprs_mas5_df, main = "MAS5 Normalization (Log2)",
        las = 2, cex.axis = 0.7, col = "lightgreen")
boxplot(exprs_gcrma_df, main = "GCRMA Normalization",
        las = 2, cex.axis = 0.7, col = "lightcoral")
par(mfrow = c(1, 1))  # Reset layout
```

### 3. Differential Expression Analysis

Identify DEGs between defined groups using linear models and empirical Bayes moderation.

```R
# Read targets file
targets <- readTargets("targets.txt")
print(targets)

# Read CEL files with targets
abatch <- ReadAffy(filenames = targets$filename)

# Perform RMA normalization
eset <- rma(abatch)

# Create design matrix
group_factor <- factor(targets$subgroup)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
print(design)

# Fit linear model
fit <- lmFit(eset, design)

# Define contrasts
contrast_matrix <- makeContrasts(
    Diff = "Poor_Prognosis - Good_Prognosis",
    levels = design
)
print(contrast_matrix)

# Apply contrasts and compute statistics
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract DEGs
all_DEGs <- topTable(fit2, coef = "Diff", number = Inf, adjust.method = "BH")
write.table(all_DEGs, file = "ALL_DEGs.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# Significant DEGs with thresholds
significant_DEGs <- topTable(fit2, coef = "Diff", number = Inf,
                             adjust.method = "BH", p.value = 0.05, lfc = 1)
write.table(significant_DEGs, file = "Significant_DEGs.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# Summarize results
results <- decideTests(fit2)
summary(results)
vennDiagram(results)
```

### 4. Heatmap Visualization

Generate a heatmap of significant DEGs to visualize expression patterns.

```R
# Extract expression data for significant DEGs
DEG_genes <- rownames(significant_DEGs)
DEG_exprs <- exprs(eset)[DEG_genes, ]

# Scale data
scale_rows <- function(x) {
    row_means <- rowMeans(x, na.rm = TRUE)
    row_sds <- apply(x, 1, sd, na.rm = TRUE)
    scaled_x <- sweep(x, 1, row_means, "-")
    scaled_x <- sweep(scaled_x, 1, row_sds, "/")
    return(scaled_x)
}
scaled_DEG_exprs <- scale_rows(DEG_exprs)

# Sample annotations
sample_conditions <- targets$subgroup
annotation_df <- data.frame(Group = sample_conditions)
rownames(annotation_df) <- colnames(scaled_DEG_exprs)
annotation_colors <- list(Group = c("Poor_Prognosis" = "red",
                                    "Good_Prognosis" = "blue"))

# Create heatmap annotation
heatmap_annotation <- HeatmapAnnotation(df = annotation_df,
                                        col = annotation_colors)

# Generate heatmap
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

# Save heatmap
png(filename = "DEG_heatmap.png", width = 10, height = 8, units = "in", res = 300)
draw(heatmap)
dev.off()

# Display heatmap in R
draw(heatmap)
```

---

## Results

- **Normalized Data Files**:
  - `affymetrix_rma.txt`
  - `affymetrix_mas5.txt`
  - `affymetrix_gcrma.txt`
- **Differential Expression Results**:
  - `ALL_DEGs.txt`: Contains all genes with differential expression statistics.
  - `Significant_DEGs.txt`: Contains DEGs meeting the specified thresholds (adjusted p-value < 0.05, logFC > 1).
- **Visualizations**:
  - `DEG_heatmap.png`: Heatmap of significant DEGs across samples.

---

## Visualization

### Heatmap of Significant DEGs

The heatmap illustrates the expression patterns of significant DEGs between the defined groups, providing insights into gene regulation and potential biomarkers.

![DEG Heatmap](DEG_heatmap.png)

- **Annotations**: Samples are annotated by group (e.g., Poor Prognosis vs. Good Prognosis).
- **Clustering**: Both genes and samples are hierarchically clustered to reveal patterns.
- **Color Scale**: Represents the Z-score of expression levels, with blue indicating lower expression and red indicating higher expression.

---

## Contributing

Contributions are highly appreciated! To contribute:

1. **Fork** the repository.
2. **Clone** your fork:

   ```bash
   git clone https://github.com/your-username/microarray-analysis-pipeline.git
   ```

3. **Create a branch** for your feature:

   ```bash
   git checkout -b feature/your-feature
   ```

4. **Commit** your changes:

   ```bash
   git commit -am 'Add your feature'
   ```

5. **Push** to your branch:

   ```bash
   git push origin feature/your-feature
   ```

6. **Open a Pull Request** on GitHub.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Luis Fernando Nagano**

- **Email**: [nagano.luis@gmail.com](mailto:luis.nagano@example.com)
- **LinkedIn**: [www.linkedin.com/in/luis-fernando-nagano-7585b82a8](https://www.linkedin.com/in/luisfernandonagano)
- **GitHub**: [https://github.com/LuisNagano](https://github.com/luisfernandonagano)

---

## Acknowledgements

- **Bioconductor Project**: For providing open-source bioinformatics tools.
- **R Development Core Team**: For the R programming language.
- **Affymetrix**: For microarray technology and data formats.
- **Community Contributors**: For continuous improvements and support.
