# STDistance - Spatial Transcriptomics Distance Calculation and Visualization
[![DOI](https://zenodo.org/badge/999352906.svg)](https://doi.org/10.5281/zenodo.17149296)

## What's New — v0.6.7

- **Direct RDS input**: `calculate_nearest_distances()` now accepts a file path to a Seurat v5 RDS file (or a live Seurat object) as `spatial_data`. Spatial coordinates are extracted automatically — no manual CSV preparation required.
- **Sample-aware distance calculation**: New `sample_col` parameter (e.g. `sample_col = "orig.ident"`) splits multi-sample objects by sample before calculating distances, preventing cross-sample nearest-neighbour matches in pooled datasets.
- **New helper `seurat_to_spatial_df()`**: Converts a Seurat v5 object to the flat data frame format used internally. Exported so users can inspect or modify the extracted data before passing it onward.
- **Nearest-target barcode in output**: Each target type now produces two side-by-side columns — `<target>` (distance) and `<target>_barcode` (barcode of the nearest target spot) — making it straightforward to `merge()` results with metadata.

## Description

STDistance is an R package designed for analyzing spatial relationships between cell types in spatial transcriptomics data. It calculates nearest neighbor distances between specified cell types and provides comprehensive visualization tools to explore spatial patterns. The package is particularly useful for studying cell-cell interactions, immune microenvironment characterization, and spatial organization of tissues.

Key features include:
- Distance calculation between reference and target cell types
- Direct input from Seurat v5 RDS files with automatic coordinate extraction
- Sample-aware distance calculation for multi-sample datasets
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
devtools::install_github("PrinceWang2018/ST_Distance")
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

For RDS / Seurat input (v0.6.7+):
- Seurat (≥ 5.0)
- SeuratObject (≥ 5.0)

## Input File Preparation

STDistance supports two input routes. Use whichever fits your data.

### Option A — Seurat v5 RDS file (recommended, v0.6.7+)

Pass a file path to a Seurat v5 RDS directly as `spatial_data`. The package extracts spot coordinates and cell type labels automatically. No manual CSV export is needed.

```r
# Minimal example — coordinates and labels are extracted from the RDS automatically
distance_results <- calculate_nearest_distances(
  spatial_data   = "Demo_SP6_SP8_v5.RDS",
  reference_type = "Macrophage",
  target_types   = c("Epithelial_cells_A", "Epithelial_cells_B",
                     "Epithelial_cells_C", "Epithelial_cells_D"),
  type_col       = "celltype_ABCDepi",
  sample_col     = "orig.ident"   # recommended for multi-sample objects
)
```

> **Multi-sample objects**: always set `sample_col = "orig.ident"` (or whichever metadata column identifies the sample). Without it, a reference spot in sample A can be matched to the nearest target spot in sample B, producing physically meaningless distances.

If your Seurat object was built with an older version of Seurat, update it first:

```r
library(Seurat)
obj <- readRDS("Demo_SP6_SP8.RDS")
obj <- UpdateSeuratObject(obj)
saveRDS(obj, "Demo_SP6_SP8_v5.RDS")
```

The helper `seurat_to_spatial_df()` is also available if you want to inspect the extracted data frame before running the distance calculation:

```r
df <- seurat_to_spatial_df(obj,
  type_col  = "celltype_ABCDepi",
  meta_cols = "orig.ident")   # carry over any extra metadata columns
head(df)
```

### Option B — CSV files (classic workflow)

STDistance accepts two plain CSV files:

1. **Spatial coordinates file** (e.g., `tissue_positions.csv`):
   
   - Must include columns for barcode, x and y coordinates
   - For multiple samples, include "Sample" and "Newbarcode" columns
   - Example format:
     ```
     barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres,Sample,Sampleid,Newbarcode
     AAACCCAAGGCTTTCA-1_1,1,50,102,4950,10020,Sample1,1,AAACCCAAGGCTTTCA-1_1
     ```
   
2. **Metadata file** (e.g., `metadata.csv`):
   
   - Must include `orig.ident` and a cell type column
   - First column must match the barcode/Newbarcode column in the coordinates file
   - May include gene expression or splicing index values
   - Example format:
     ```
     ,orig.ident,nCount_Spatial,nFeature_Spatial,celltype_ABCDepi,gen2_SPLIz_numeric
     AAATCGTGTACCACAA-1_6,SP6,5403,2647,Epithelial_cells_B,0.96565309
     AACCCTACTGTCAATA-1_6,SP6,40683,8876,Epithelial_cells_A,-0.300446291
     ```
   - Can be exported from a Seurat object using:
     ```r
     library(Seurat)
     RDS <- readRDS("Demo_SP6_SP8_v5.RDS")
     write.csv(RDS@meta.data, file = "Demo_SP6_SP8_metadata.csv", quote = FALSE)
     ```


## Basic Work flow

Demo data is available in the `./inst/extdata/` folder of the R package installed from GitHub.

### Workflow A — RDS input (v0.6.7+, recommended)

Download the Seurat v5 demo RDS from the demo repository:

```shell
wget https://github.com/PrinceWang2018/ST_Distance_demo/raw/master/Demo_SP6_SP8_v5.RDS
```

#### 1. Calculate nearest distances directly from the RDS

```r
library(STDistance)

distance_results <- calculate_nearest_distances(
  spatial_data   = "Demo_SP6_SP8_v5.RDS",
  reference_type = "Macrophage",
  target_types   = c("Epithelial_cells_A", "Epithelial_cells_B",
                     "Epithelial_cells_C", "Epithelial_cells_D"),
  type_col       = "celltype_ABCDepi",
  sample_col     = "orig.ident"   # enables within-sample matching
)
```

The output contains one row per reference spot, with distance and nearest-target barcode columns side by side for each target type:

```
               barcode orig.ident  Epithelial_cells_A  Epithelial_cells_A_barcode  ...
ATAGTGAAGCGTTCTC-1_6      SP6          122.3          ATAAGTTACCGCGACG-1_6       ...
```

The `_barcode` columns make it easy to merge distance results back with any metadata:

```r
# Example: attach metadata for the nearest Epithelial_cells_A spot
merged <- merge(distance_results, metadata,
                by.x = "Epithelial_cells_A_barcode",
                by.y = "row.names")
```

---

### Workflow B — CSV input (classic)

Demo CSVs are available in the `./inst/extdata/` folder of the installed package, or can be downloaded directly:

```shell
wget https://github.com/PrinceWang2018/ST_Distance/raw/master/inst/extdata/Demo_SP6_SP8_metadata.csv
wget https://github.com/PrinceWang2018/ST_Distance/raw/master/inst/extdata/Demo_SP6_SP8_tissue_positions.csv
```

#### 1. Load required packages and data

```r
library(STDistance)
setwd("R package dir or work dir")

# Load spatial coordinates
tissue_posi <- read.csv(system.file("extdata/Demo_SP6_SP8_tissue_positions.csv",
                                    package = "STDistance"), header = TRUE)
# Load metadata
metadata <- read.csv(system.file("extdata/Demo_SP6_SP8_metadata.csv",
                                  package = "STDistance"), header = TRUE, row.names = 1)
```

#### 2. Normalize spatial coordinates

```r
tissue_posi_normalized <- normalize_spatial(tissue_posi)
```

#### 3. Merge spatial and metadata information

```r
posi <- merge(
  x     = tissue_posi_normalized,
  y     = metadata,
  by.x  = "Newbarcode",
  by.y  = "row.names",
  all.y = TRUE
)
```

#### 4. Calculate nearest distances between cell types

```r
distance_results <- calculate_nearest_distances(
  posi,
  reference_type = "Macrophage",
  target_types   = c("Epithelial_cells_A", "Epithelial_cells_B", "Epithelial_cells_C"),
  x_col          = "pxl_row_in_fullres",
  y_col          = "pxl_col_in_fullres",
  id_col         = "Newbarcode",
  type_col       = "celltype_ABCDepi",
  sample_col     = "orig.ident"   # recommended for multi-sample data frames
)
```

#### 5. Compare distance among subgroups

```r
plot_distance_boxplot(
  distance_results,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "log10",
  palette = "Dark2"
)
```

#### 6. Create radial network visualization

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

#### 7. Spatial visualization of interactions

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
  gradient_palette = "C", # Viridis color palette option: "A","B",C"...
  point_size = 1.5,
  point_alpha = 0.9
)
```

#### 8. Correlation analysis

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
   - Spatial file: Must have coordinate columns (default `"pxl_row_in_fullres"` and `"pxl_col_in_fullres"`)
   - Metadata: Must have a cell type annotation column

2. **No distances calculated**: Check that:
   - Your reference and target cell types exist in the metadata
   - The `type_col` parameter correctly specifies your cell type column name
   - There are actually cells of the specified types in your sample

3. **RDS input — Seurat version mismatch**: If you see warnings about `SeuratObject` being built under a different R or Matrix version, update the object before saving:
   ```r
   library(Seurat)
   obj <- UpdateSeuratObject(readRDS("your_file.RDS"))
   saveRDS(obj, "your_file_v5.RDS")
   ```

4. **Cross-sample nearest-neighbour matches**: If distance results show reference spots from one sample matched to target spots in another sample, you have not set `sample_col`. Add `sample_col = "orig.ident"` (or whichever column identifies the sample) to `calculate_nearest_distances()`.

5. **Visualization issues**:
   - For crowded plots, try adjusting `point_size` and `alpha` parameters
   - For color issues, specify a custom `color_palette`

## Citation

If you use STDistance in your research, please cite:

> Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer. bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282

## Contact

For questions or issues, please contact:
- 970214035yl@gmail.com or wangzixiang@sdu.edu.cn
- GitHub issues: https://github.com/PrinceWang2018/ST_Distance/issues

## License

GPL-3 © Zixiang Wang, Lei Yang, Zhaojian Liu
