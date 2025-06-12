#' Visualize spatial relationships between cell types
#'
#' @param spatial_data Spatial coordinates data frame
#' @param sample Sample name in the spatial transcriptome data
#' @param reference_type Reference cell type
#' @param target_type Target cell type
#' @param x_col Column name for x-coordinates
#' @param y_col Column name for y-coordinates
#' @param type_col Column name for cell type information
#' @param color_palette Named vector of colors for cell types
#' @param alpha Transparency level for points and lines
#' @return A ggplot object showing the spatial relationships
#' @export
#' @import ggplot2
#' @import dplyr
#' @examples
#' visualize_spatial_network(posi, sample="SP8", reference_type="Macrophage",
#'                           target_type="Epithelial_cells_A",
#'                           x_col = "pxl_row_in_fullres",
#'                           y_col = "pxl_col_in_fullres",
#'                           type_col = "celltype_ABCDepi",
#'                           color_palette = c("Macrophage" = "#90ee90",
#'                          "Epithelial_cells_A" = "#377EB8"))

visualize_spatial_network <- function(spatial_data, sample, reference_type, target_type,
                                      x_col = "pxl_row_in_fullres",
                                      y_col = "pxl_col_in_fullres",
                                      type_col = "Epi_strom",
                                      color_palette = c("Macrophage" = "#90ee90",
                                                        "Epithelial_cells_A" = "#377EB8"),
                                      alpha = 0.7) {

  #Choose sample
  spatial_data_sample <- spatial_data[spatial_data[["Sample"]] == sample, ]

  # Filter and prepare data
  reference_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == reference_type, ]
  target_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == target_type, ]

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

  # Prepare network lines data
  network_lines <- data.frame(
    x1 = reference_cells[[x_col]],
    y1 = reference_cells[[y_col]],
    x2 = target_cells[[x_col]][sapply(distances, function(x) x$idx)],
    y2 = target_cells[[y_col]][sapply(distances, function(x) x$idx)],
    group = reference_type
  )

  # Combine data for plotting
  plot_data <- rbind(
    reference_cells %>% mutate(group = reference_type),
    target_cells %>% mutate(group = target_type)
  )

  # Create plot
  ggplot() +
    geom_point(data = plot_data,
               aes(x = .data[[x_col]], y = .data[[y_col]], color = group),
               size = 1, alpha = alpha) +
    geom_segment(data = network_lines,
                 aes(x = x1, y = y1, xend = x2, yend = y2, color = group),
                 size = 0.3, alpha = alpha * 0.7) +
    scale_color_manual(values = color_palette) +
    labs(title = sample, x = "X Position", y = "Y Position") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.title = element_blank())
}
