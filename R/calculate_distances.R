#' Calculate nearest distances between cell types
#'
#' @param spatial_data A data frame containing spatial coordinates and cell type info
#' @param reference_type The reference cell type to calculate distances from
#' @param target_types Vector of target cell types to calculate distances to
#' @param x_col Column name for x-coordinates
#' @param y_col Column name for y-coordinates
#' @param id_col Column name for cell identifiers
#' @param type_col Column name for cell type information
#' @return A data frame with nearest distances for each reference cell
#' @export
#' @examples
#' calculate_nearest_distances(posi, reference_type="Macrophage",
#'                  target_types=c("Epithelial_cells_A","Epithelial_cells_B",
#'                  "Epithelial_cells_C","Epithelial_cells_D"),
#'                  id_col = "Newbarcode",
#'                  type_col = "celltype_ABCDepi")

calculate_nearest_distances <- function(spatial_data, reference_type, target_types,
                                        x_col = "pxl_row_in_fullres",
                                        y_col = "pxl_col_in_fullres",
                                        id_col = "barcode",
                                        type_col = "Epi_strom") {

  # Input validation
  required_cols <- c(x_col, y_col, id_col, type_col)
  if(!all(required_cols %in% names(spatial_data))) {
    stop(paste("Missing required columns:",
               paste(setdiff(required_cols, names(spatial_data)), collapse = ", ")))
  }

  # Split data into reference and target
  reference_cells <- spatial_data[spatial_data[[type_col]] == reference_type, ]
  target_data <- lapply(target_types, function(ttype) {
    spatial_data[spatial_data[[type_col]] == ttype, ]
  })
  names(target_data) <- target_types

  # Initialize result dataframe
  result <- data.frame(matrix(ncol = length(target_types) + 1, nrow = nrow(reference_cells)))
  colnames(result) <- c(id_col, target_types)
  result[[id_col]] <- reference_cells[[id_col]]

  # Calculate nearest distances for each target type
  for(ttype in target_types) {
    target_df <- target_data[[ttype]]

    if(nrow(target_df) == 0) {
      result[[ttype]] <- NA
      next
    }

    distances <- apply(reference_cells, 1, function(ref_row) {
      ref_x <- as.numeric(ref_row[x_col])
      ref_y <- as.numeric(ref_row[y_col])

      # Calculate distances to all target cells
      target_distances <- sqrt(
        (as.numeric(target_df[[x_col]]) - ref_x)^2 +
          (as.numeric(target_df[[y_col]]) - ref_y)^2
      )

      # Return minimum distance (excluding zero for same cell)
      min(target_distances[target_distances > 0])
    })

    result[[ttype]] <- distances
  }

  return(result)
}
