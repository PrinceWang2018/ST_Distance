#' Convert a Seurat v5 spatial object to a flat spatial data frame
#'
#' Extracts spot coordinates from all images in a Seurat object and joins them
#' with selected metadata columns, producing the flat data frame expected by
#' \code{\link{calculate_nearest_distances}}.
#'
#' @param seurat_obj A Seurat v5 object with at least one spatial image.
#' @param type_col Name of the metadata column containing cell/spot type labels.
#' @param id_col Name to use for the barcode column in the output (default \code{"barcode"}).
#' @param meta_cols Additional metadata columns to carry over (character vector, optional).
#'   Pass \code{"orig.ident"} here to retain the sample identifier for use with
#'   \code{sample_col} in \code{\link{calculate_nearest_distances}}.
#' @return A data frame with columns \code{id_col}, \code{"x"}, \code{"y"},
#'   \code{type_col}, and any \code{meta_cols} requested.
#' @export
seurat_to_spatial_df <- function(seurat_obj,
                                  type_col,
                                  id_col    = "barcode",
                                  meta_cols = NULL) {

  if (!inherits(seurat_obj, "Seurat"))
    stop("'seurat_obj' must be a Seurat object.")
  if (!type_col %in% colnames(seurat_obj@meta.data))
    stop("type_col '", type_col, "' not found in Seurat metadata.")

  # Pull coordinates from every image (sample) in the object
  coord_list <- lapply(SeuratObject::Images(seurat_obj), function(img_name) {
    cents <- seurat_obj@images[[img_name]]@boundaries$centroids
    data.frame(
      barcode = cents@cells,
      x       = cents@coords[, "x"],
      y       = cents@coords[, "y"],
      stringsAsFactors = FALSE
    )
  })
  coords_df <- do.call(rbind, coord_list)
  colnames(coords_df)[colnames(coords_df) == "barcode"] <- id_col

  # Attach cell type and any extra metadata columns
  keep_cols <- unique(c(type_col, meta_cols))
  missing   <- setdiff(keep_cols, colnames(seurat_obj@meta.data))
  if (length(missing) > 0)
    stop("Columns not found in Seurat metadata: ", paste(missing, collapse = ", "))

  meta <- seurat_obj@meta.data[coords_df[[id_col]], keep_cols, drop = FALSE]
  cbind(coords_df, meta)
}


# Internal: distance calculation for a single pool of cells (no sample splitting).
.calculate_distances_core <- function(spatial_data, reference_type, target_types,
                                       x_col, y_col, id_col, type_col) {

  reference_cells <- spatial_data[spatial_data[[type_col]] == reference_type, ]
  target_data     <- lapply(target_types, function(ttype) {
    spatial_data[spatial_data[[type_col]] == ttype, ]
  })
  names(target_data) <- target_types

  result           <- data.frame(matrix(ncol = 1, nrow = nrow(reference_cells)),
                                 stringsAsFactors = FALSE)
  colnames(result) <- id_col
  result[[id_col]] <- reference_cells[[id_col]]

  for (ttype in target_types) {
    target_df   <- target_data[[ttype]]
    barcode_col <- paste0(ttype, "_barcode")

    if (nrow(target_df) == 0) {
      result[[ttype]]       <- NA_real_
      result[[barcode_col]] <- NA_character_
      next
    }

    hits <- lapply(seq_len(nrow(reference_cells)), function(i) {
      ref_x <- as.numeric(reference_cells[[x_col]][i])
      ref_y <- as.numeric(reference_cells[[y_col]][i])

      dists <- sqrt(
        (as.numeric(target_df[[x_col]]) - ref_x)^2 +
          (as.numeric(target_df[[y_col]]) - ref_y)^2
      )

      valid <- dists > 0
      if (!any(valid)) return(list(dist = NA_real_, barcode = NA_character_))
      best <- which(valid)[which.min(dists[valid])]
      list(dist = dists[best], barcode = target_df[[id_col]][best])
    })

    result[[ttype]]       <- sapply(hits, `[[`, "dist")
    result[[barcode_col]] <- sapply(hits, `[[`, "barcode")
  }

  result
}


#' Calculate nearest distances between cell types
#'
#' @param spatial_data A data frame, a Seurat v5 object, or a file path to an
#'   RDS file containing a Seurat v5 object. When a Seurat object or RDS path
#'   is supplied, coordinates are extracted automatically via
#'   \code{\link{seurat_to_spatial_df}} and \code{x_col}/\code{y_col} default
#'   to \code{"x"} and \code{"y"}.
#' @param reference_type The reference cell type to calculate distances from.
#' @param target_types Vector of target cell types to calculate distances to.
#' @param x_col Column name for x-coordinates (ignored when input is Seurat/RDS).
#' @param y_col Column name for y-coordinates (ignored when input is Seurat/RDS).
#' @param id_col Column name for cell/spot identifiers.
#' @param type_col Column name for cell type information.
#' @param sample_col Column name identifying the sample each spot belongs to
#'   (e.g. \code{"orig.ident"} in Seurat objects). When provided, distances are
#'   calculated \strong{within each sample independently}, preventing
#'   cross-sample nearest-neighbour matches in multi-sample datasets. The sample
#'   label is retained as a second column in the output to facilitate
#'   per-sample downstream analysis. Default \code{NULL} (no sample splitting).
#' @return A data frame with one row per reference cell. Columns are:
#'   \code{id_col}, optionally \code{sample_col}, then for each target type a
#'   distance column (\code{<target>}) and a nearest-barcode column
#'   (\code{<target>_barcode}).
#' @export
#' @examples
#' # From a plain data frame (no sample awareness)
#' calculate_nearest_distances(posi,
#'   reference_type = "Macrophage",
#'   target_types   = c("Epithelial_cells_A", "Epithelial_cells_B"),
#'   id_col         = "Newbarcode",
#'   type_col       = "celltype_ABCDepi")
#'
#' # From a Seurat v5 RDS file, with sample-aware matching
#' calculate_nearest_distances("Demo_SP6_SP8_v5.RDS",
#'   reference_type = "Macrophage",
#'   target_types   = c("Epithelial_cells_A", "Epithelial_cells_B",
#'                      "Epithelial_cells_C", "Epithelial_cells_D"),
#'   type_col       = "celltype_ABCDepi",
#'   sample_col     = "orig.ident")

calculate_nearest_distances <- function(spatial_data, reference_type, target_types,
                                        x_col      = "pxl_row_in_fullres",
                                        y_col      = "pxl_col_in_fullres",
                                        id_col     = "barcode",
                                        type_col   = "Epi_strom",
                                        sample_col = NULL) {

  # ── Accept RDS file path ────────────────────────────────────────────────────
  if (is.character(spatial_data) && length(spatial_data) == 1) {
    if (!file.exists(spatial_data))
      stop("RDS file not found: ", spatial_data)
    spatial_data <- readRDS(spatial_data)
  }

  # ── Accept Seurat object ────────────────────────────────────────────────────
  if (inherits(spatial_data, "Seurat")) {
    spatial_data <- seurat_to_spatial_df(spatial_data,
                                         type_col  = type_col,
                                         id_col    = id_col,
                                         meta_cols = sample_col)
    x_col <- "x"
    y_col <- "y"
  }

  # ── Input validation ────────────────────────────────────────────────────────
  required_cols <- c(x_col, y_col, id_col, type_col, sample_col)
  if (!all(required_cols %in% names(spatial_data)))
    stop("Missing required columns: ",
         paste(setdiff(required_cols, names(spatial_data)), collapse = ", "))

  # ── Sample-aware: calculate within each sample, then combine ────────────────
  if (!is.null(sample_col)) {
    samples <- unique(spatial_data[[sample_col]])
    per_sample <- lapply(samples, function(s) {
      sub_data   <- spatial_data[spatial_data[[sample_col]] == s, ]
      sub_result <- .calculate_distances_core(sub_data, reference_type, target_types,
                                              x_col, y_col, id_col, type_col)
      sub_result[[sample_col]] <- s
      sub_result
    })
    result <- do.call(rbind, per_sample)
    # Place sample_col right after id_col for readability
    result <- result[, c(id_col, sample_col,
                         setdiff(colnames(result), c(id_col, sample_col)))]
    return(result)
  }

  # ── Single pool (no sample awareness) ──────────────────────────────────────
  .calculate_distances_core(spatial_data, reference_type, target_types,
                             x_col, y_col, id_col, type_col)
}
