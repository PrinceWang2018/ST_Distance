#' Visualize spatial relationships between multiple cell types
#'
#' @param spatial_data Spatial coordinates data frame
#' @param sample Sample name in the spatial transcriptome data
#' @param reference_type Reference cell type (character vector of length 1)
#' @param target_types Target cell type(s) (character vector of 1 or more)
#' @param x_col Column name for x-coordinates
#' @param y_col Column name for y-coordinates
#' @param type_col Column name for cell type information
#' @param color_palette Named vector of colors for cell types
#' @param point_alpha Transparency level for points
#' @param line_alpha Transparency level for connection lines
#' @param point_size Size of points in plot
#' @param line_width Width of connection lines
#' @param show_legend Logical, whether to show legend
#' @return A ggplot object showing the spatial relationships
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom scales hue_pal
#' @examples
#' visualize_spatial_multinetwork(posi, sample="SP8",reference_type="Macrophage",
#'                      target_type=c("Epithelial_cells_A","Epithelial_cells_B"),
#'                      type_col = "celltype_ABCDepi")


visualize_spatial_multinetwork <- function(spatial_data,
                                           sample,
                                      reference_type,
                                      target_types,
                                      x_col = "pxl_row_in_fullres",
                                      y_col = "pxl_col_in_fullres",
                                      type_col = "Epi_strom",
                                      color_palette = NULL,
                                      point_alpha = 0.7,
                                      line_alpha = 0.5,
                                      point_size = 1.5,
                                      line_width = 0.3,
                                      show_legend = TRUE) {

  # Input validation
  if(length(reference_type) != 1) {
    stop("reference_type must be a single cell type")
  }

  if(length(target_types) < 1) {
    stop("At least one target_type must be specified")
  }

  required_cols <- c(x_col, y_col, type_col)
  if(!all(required_cols %in% names(spatial_data))) {
    stop(paste("Missing required columns:",
               paste(setdiff(required_cols, names(spatial_data)), collapse = ", ")))
  }

  #Choose sample
  spatial_data_sample <- spatial_data[spatial_data[["Sample"]] == sample, ]

  # Filter reference cells
  reference_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == reference_type, ]
  if(nrow(reference_cells) == 0) {
    stop(paste("No cells found for reference type:", reference_type))
  }

  # Process each target type
  network_lines_list <- list()

  for(target_type in target_types) {
    target_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == target_type, ]

    if(nrow(target_cells) == 0) {
      warning(paste("No cells found for target type:", target_type))
      next
    }

    # Calculate nearest neighbors
    distances <- apply(reference_cells, 1, function(ref_row) {
      ref_x <- as.numeric(ref_row[x_col])
      ref_y <- as.numeric(ref_row[y_col])

      target_distances <- sqrt(
        (as.numeric(target_cells[[x_col]]) - ref_x)^2 +
          (as.numeric(target_cells[[y_col]]) - ref_y)^2
      )

      list(
        idx = which.min(target_distances),
        distance = min(target_distances)
      )
    })

    # Prepare network lines data for this target type
    network_lines <- data.frame(
      x1 = reference_cells[[x_col]],
      y1 = reference_cells[[y_col]],
      x2 = target_cells[[x_col]][sapply(distances, function(x) x$idx)],
      y2 = target_cells[[y_col]][sapply(distances, function(x) x$idx)],
      group = paste(reference_type, target_type, sep = " to "),
      target_group = target_type,
      stringsAsFactors = FALSE
    )

    network_lines_list[[target_type]] <- network_lines
  }

  # Combine all network lines
  all_network_lines <- do.call(rbind, network_lines_list)

  if(is.null(all_network_lines)) {
    stop("No valid target types with cells found")
  }

  # Prepare plotting data - include reference and all target types
  plot_types <- c(reference_type, target_types)
  plot_data <- spatial_data_sample[spatial_data_sample[[type_col]] %in% plot_types, ]
  plot_data$group <- plot_data[[type_col]]

  # Set default color palette if not provided
  if(is.null(color_palette)) {
    color_palette <- scales::hue_pal()(length(unique(plot_data$group)))
    names(color_palette) <- unique(plot_data$group)

    # Add colors for connection lines (mix of reference and target colors)
    line_colors <- color_palette[target_types]
    names(line_colors) <- paste(reference_type, target_types, sep = " to ")
    color_palette <- c(color_palette, line_colors)
  }

  # Create plot
  p <- ggplot() +
    # Plot all points first
    geom_point(
      data = plot_data,
      aes(x = .data[[x_col]], y = .data[[y_col]], color = group),
      size = point_size,
      alpha = point_alpha
    ) +
    # Plot connection lines
    geom_segment(
      data = all_network_lines,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = group),
      linewidth = line_width,
      alpha = line_alpha
    ) +
    labs(x = "X Position", y = "Y Position", color = "Cell Type") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  # Apply color palette
  if(!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  }

  # Hide legend if requested
  if(!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}
