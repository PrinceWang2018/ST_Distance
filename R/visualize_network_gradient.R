#' Visualize spatial network with expression gradient
#'
#' @param spatial_data Spatial coordinates data frame containing cell types and expression values
#' @param sample Sample name in the spatial transcriptome data
#' @param gradient_type Cell type to show with expression gradient coloring
#' @param fixed_type Cell type to show in fixed color (default gray)
#' @param expression_col Column name containing expression values (default "gen2_SPLIz_numeric")
#' @param x_col Column name for x-coordinates (default "pxl_row_in_fullres")
#' @param y_col Column name for y-coordinates (default "pxl_col_in_fullres")
#' @param type_col Column name for cell type information (default "Epi_strom")
#' @param fixed_color Color for the fixed cell type (default "#A9C6D9" - light gray-blue)
#' @param line_color Color for connection lines (default "#666666" - dark gray)
#' @param gradient_palette Color palette for expression gradient (default viridis option "C")
#' @param point_size Size of points (default 1)
#' @param point_alpha Transparency of points (default 0.8)
#' @param line_width Width of connection lines (default 0.3)
#' @param line_alpha Transparency of connection lines (default 0.6)
#' @param show_legend Logical whether to show legend (default TRUE)
#' @param legend_title Title for the legend (default "Expression")
#' @param grid_major_color Color for major grid lines (default "gray90")
#' @param grid_minor_color Color for minor grid lines (default "gray95")
#' @param border_color Color for plot border (default "black")
#' @param background_color Color for plot background (default "white")
#' @return A ggplot object showing the spatial relationships with expression gradient
#' @import dplyr
#' @export
#' @examples
#' visualize_spatial_gradient(spatial_data = posi,
#'                            sample="SP8",
#'                            gradient_type = "Epithelial_cells_A",
#'                            fixed_type = "Macrophage",
#'                            expression_col = "gen2_SPLIz_numeric",
#'                            type_col = "celltype_ABCDepi",
#'                            legend_title = "Expression",
#'                            background_color = "white")

visualize_spatial_gradient <- function(spatial_data,
                                       sample,
                                       gradient_type,
                                       fixed_type,
                                       expression_col = "gen2_SPLIz_numeric",
                                       x_col = "pxl_row_in_fullres",
                                       y_col = "pxl_col_in_fullres",
                                       type_col = "Epi_strom",
                                       fixed_color = "#A9C6D9",
                                       line_color = "#666666",
                                       gradient_palette = "C",
                                       point_size = 1,
                                       point_alpha = 0.8,
                                       line_width = 0.3,
                                       line_alpha = 0.6,
                                       show_legend = TRUE,
                                       legend_title = "Expression",
                                       grid_major_color = "gray90",
                                       grid_minor_color = "gray95",
                                       border_color = "black",
                                       background_color = "white") {

  # Input validation
  required_cols <- c(x_col, y_col, type_col, expression_col)
  missing_cols <- setdiff(required_cols, names(spatial_data))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  if(!gradient_type %in% spatial_data[[type_col]]) {
    stop(paste("gradient_type", gradient_type, "not found in", type_col))
  }

  if(!fixed_type %in% spatial_data[[type_col]]) {
    stop(paste("fixed_type", fixed_type, "not found in", type_col))
  }

  #Choose sample
  spatial_data_sample <- spatial_data[spatial_data[["Sample"]] == sample, ]

  # Filter cells for the two types
  gradient_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == gradient_type, ]
  fixed_cells <- spatial_data_sample[spatial_data_sample[[type_col]] == fixed_type, ]

  if(nrow(gradient_cells) == 0) {
    stop(paste("No cells found for gradient_type:", gradient_type))
  }

  if(nrow(fixed_cells) == 0) {
    stop(paste("No cells found for fixed_type:", fixed_type))
  }

  # Calculate nearest neighbors (gradient to fixed)
  distances <- apply(gradient_cells, 1, function(grad_row) {
    grad_x <- as.numeric(grad_row[x_col])
    grad_y <- as.numeric(grad_row[y_col])

    fixed_distances <- sqrt(
      (as.numeric(fixed_cells[[x_col]]) - grad_x)^2 +
        (as.numeric(fixed_cells[[y_col]]) - grad_y)^2
    )

    list(
      idx = which.min(fixed_distances),
      distance = min(fixed_distances)
    )
  })

  # Prepare network lines data
  network_lines <- data.frame(
    x1 = gradient_cells[[x_col]],
    y1 = gradient_cells[[y_col]],
    x2 = fixed_cells[[x_col]][sapply(distances, function(x) x$idx)],
    y2 = fixed_cells[[y_col]][sapply(distances, function(x) x$idx)]
  )

  # Prepare plotting data
  plot_data <- rbind(
    gradient_cells %>% mutate(cell_group = "gradient"),
    fixed_cells %>% mutate(cell_group = "fixed")
  )

  # Create the plot
  p <- ggplot() +
    # Plot fixed type cells (gray)
    geom_point(
      data = plot_data %>% filter(cell_group == "fixed"),
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      color = fixed_color,
      size = point_size,
      alpha = point_alpha
    ) +
    # Plot gradient type cells (colored by expression)
    geom_point(
      data = plot_data %>% filter(cell_group == "gradient"),
      aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[expression_col]]),
      size = point_size,
      alpha = point_alpha
    ) +
    # Plot connection lines
    geom_segment(
      data = network_lines,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = line_color,
      linewidth = line_width,
      alpha = line_alpha
    ) +
    # Apply viridis color scale
    scale_color_viridis_c(
      option = gradient_palette,
      na.value = fixed_color,
      name = legend_title
    ) +
    labs(x = "Row", y = "Column") +
    theme(
      panel.grid.major = element_line(color = grid_major_color, linewidth = 0.3),
      panel.grid.minor = element_line(color = grid_minor_color, linewidth = 0.3),
      panel.grid = element_blank(),
      panel.border = element_rect(color = border_color, fill = NA, linewidth = 0.8),
      panel.background = element_rect(fill = background_color),
      legend.position = ifelse(show_legend, "right", "none")
    )

  return(p)
}
