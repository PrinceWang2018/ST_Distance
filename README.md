# STDistance - Spatial Transcriptomics Distance Calculation and Visualization

## Description

STDistance is an R package designed for analyzing spatial relationships between cell types in spatial transcriptomics data. It calculates nearest neighbor distances between specified cell types and provides comprehensive visualization tools to explore spatial patterns. The package is particularly useful for studying cell-cell interactions, immune microenvironment characterization, and spatial organization of tissues.

Key features include:
- Distance calculation between reference and target cell types
- Boxplot visualization of distance distributions
- Radial network visualization of spatial relationships
- Spatial mapping of cell type interactions
- Correlation analysis between gene expression and spatial distances

## Installation

```r
# Install from CRAN
install.packages("STDistance")

# Or install the development version from GitHub:
# install.packages("devtools")
library(devtools)
devtools::install_github("PrinceWang2018/STDistance")
```

## Dependencies

STDistance requires R (≥ 4.0.0) and depends on the following packages:
- dplyr
- ggplot2
- Hmisc
- scales
- stats
- RColorBrewer
- tidyr

## Input File Preparation

STDistance requires two input files:

1. **Spatial coordinates file** (e.g., `tissue_positions.csv`):
   
   - Should contain spatial coordinates of spots/cells
   - Must include columns for barcode, x and y coordinates
   - For multiple samples, include "Sample" and "Newbarcode" columns
   - Example format:
     ```
     barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres,Sample,Sampleid,Newbarcode
     AAACCCAAGGCTTTCA-1_1,1,50,102,4950,10020,Sample1,1,AAACCCAAGGCTTTCA-1_1
     ```
   
2. **Metadata file** (e.g., `metadata.csv`):
   
   - Should contain cell type annotations and any expression metrics
   
   - Must include: orig.ident, celltype columns
   
   - The first colume must match the barcode/newbarcode column in tissue_positions.csv
   
   - May include gene expression or splicing index values
   
   - Example format:
   
     ```
     ,orig.ident,nCount_Spatial,nFeature_Spatial,nCount_SCT,nFeature_SCT,integrated_snn_res.0.8,seurat_clusters,celltype_ABCDepi,gen2_SPLIz_numeric
     AAATCGTGTACCACAA-1_6,SP6,5403,2647,6486,2601,5,5,Epithelial_cells_B,0.96565309
     AACCCTACTGTCAATA-1_6,SP6,40683,8876,8578,4328,4,4,Epithelial_cells_A,-0.300446291
     ```
   
   - Can be exported from Seurat object using:
     ```r
     library(Seurat)
     RDS<-readRDS("./Demodata/Demo_SP8_SP6.RDS")
     write.csv(RDS@meta.data,file = "./Demodata/Demo_SP6_SP8_metadata.csv",quote = F)
     ```

## Basic Workflow

Demo data is available in the `Demodata` folder of the R package installed from GitHub. Alternatively, you can download the files directly using the following commands:

```shell
wget https://github.com/PrinceWang2018/ST_Distance/raw/master/Demodata/Demo_SP6_SP8_metadata.csv  
wget https://github.com/PrinceWang2018/ST_Distance/raw/master/Demodata/Demo_SP6_SP8_tissue_positions.csv  
wget https://github.com/PrinceWang2018/ST_Distance/raw/master/Demodata/Demo_SP8_SP6.RDS  
```

Below is a basic workflow demonstrating how to use the demo data for reference:

### 1. Load required packages and data

```r
library(STDistance)
setwd("R package dir or work dir")

# Load spatial coordinates
tissue_posi <- read.csv("./Demodata/Demo_SP6_SP8_tissue_positions.csv", header = TRUE)

# Load metadata
metadata <- read.csv("./Demodata/Demo_SP6_SP8_metadata.csv", header = TRUE, row.names = 1)
```

### 2. Normalize spatial coordinates

```r
tissue_posi_normalized <- normalize_spatial(tissue_posi)
```

### 3. Merge spatial and metadata information

```r
posi <- merge(
  x = tissue_posi_normalized,
  y = metadata,
  by.x = "Newbarcode",
  by.y = "row.names",
  all.y = TRUE
)
```

### 4. Calculate nearest distances between cell types

```r
distance_results <- calculate_nearest_distances(
  posi,
  reference_type = "Macrophage",
  target_types = c("Epithelial_cells_A", "Epithelial_cells_B", "Epithelial_cells_C"),
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  id_col = "Newbarcode",
  type_col = "celltype_ABCDepi"
)
```

### 5. Compare distance among subgroups

```r
plot_distance_boxplot(
  distance_results,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "log10",
  palette = "Dark2"
)
```

### 6. Create radial network visualization

```r
plot_radial_distance(
  distance_results,
  id_col = "Newbarcode",
  reference_type = "Macrophage",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
)
```

### 7. Spatial visualization of interactions

Between two cell types:
```r
visualize_spatial_network(
  posi,
  sample = "SP8",
  reference_type = "Macrophage",
  target_type = "Epithelial_cells_A",
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  type_col = "celltype_ABCDepi",
  color_palette = c("Macrophage" = "#90ee90", "Epithelial_cells_A" = "#377EB8"),
  alpha = 0.7
)
```

Between reference and multiple target types:
```r
visualize_spatial_multinetwork(
  posi,
  sample = "SP8",
  reference_type = "Macrophage",
  target_type = c("Epithelial_cells_A", "Epithelial_cells_B"),
  type_col = "celltype_ABCDepi",
  color_palette = c("Macrophage" = "#90ee90",
                   "Epithelial_cells_A" = "#377EB8",
                   "Epithelial_cells_B" = "#E41A1C"),
  point_alpha = 0.7
)
```

With expression gradient:
```r
visualize_spatial_gradient(
  spatial_data = posi,
  sample = "SP8",
  gradient_type = "Epithelial_cells_A",
  fixed_type = "Macrophage",
  expression_col = "gen2_SPLIz_numeric",
  type_col = "celltype_ABCDepi",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "viridis",
  point_size = 1.5,
  point_alpha = 0.9
)
```

### 8. Correlation analysis

```r
result_correlation <- calculate_correlations(
  spatial_data = posi,
  distance_results = distance_results,
  spatial_feature = "gen2_SPLIz_numeric",
  distance_metric = "Epithelial_cells_A",
  method = "pearson",
  plot = TRUE,
  plot_title = "Correlation between Gene Expression and Distance"
)

print(paste("Correlation coefficient:", result_correlation$estimate))
print(paste("P-value:", result_correlation$p_value))
result_correlation$plot
```

## Applications

STDistance can be used for various spatial transcriptomics analyses:

1. **Cell-cell interaction networks**: Visualize spatial relationships between different cell types
2. **Immune microenvironment characterization**: Study spatial organization of immune cells
3. **Key gene effects on interactions**: Analyze how immune checkpoint genes influence cell proximity
4. **Spatial patterns and disease**: Investigate associations between cell type distributions and clinical outcomes
5. **Tumor microenvironment**: Study spatial relationships between tumor cells and stromal components

## Limitations

1. **Distance is only one aspect of cell communication**: While spatial proximity suggests potential interactions, true cell-cell communication should be validated with receptor-ligand analysis tools like CellChat or NicheNet that consider molecular expression profiles.
2. **Resolution limitations**: Current spatial transcriptomics demo have resolution limits (50μm spot size). Single-cell resolution spatial data (10xVisium HD or other poly(A) / probes captured high resolution data) will provide more accurate distance measurements and increase detection power for spatial relationships.

## Troubleshooting

Common issues and solutions:

1. **Missing columns error**: Ensure your input files contain all required columns:
   - Spatial file: Must have coordinate columns (default "pxl_row_in_fullres" and "pxl_col_in_fullres")
   - Metadata: Must have cell type annotation column

2. **No distances calculated**: Check that:
   - Your reference and target cell types exist in the metadata
   - The type_col parameter correctly specifies your cell type column name
   - There are actually cells of the specified types in your sample

3. **Visualization issues**: 
   - For crowded plots, try adjusting point_size and alpha parameters
   - For color issues, specify a custom color_palette

## Citation

If you use STDistance in your research, please cite:

> Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer. bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282

## Contact

For questions or issues, please contact:
- 970214035yl@gmail.com or wangzixiang@sdu.edu.cn
- GitHub issues: https://github.com/PrinceWang2018/STDistance/issues

## License

GPL-3 © Zixiang Wang, Lei Yang, Zhaojian Liu