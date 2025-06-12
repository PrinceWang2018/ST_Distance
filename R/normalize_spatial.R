#' Normalize spatial coordinates
#'
#' @param data A data frame containing spatial coordinates
#' @param sample_col Column name specifying sample IDs
#' @param x_col Column name for x-coordinates
#' @param y_col Column name for y-coordinates
#' @param min_value Minimum value for normalization range
#' @param max_value Maximum value for normalization range
#' @return A data frame with normalized coordinates
#' @export
#' @import dplyr
#' @examples
#' tissue_posi_normalized<-normalize_spatial(tissue_posi)

normalize_spatial <- function(data, sample_col = "Sample",
                              x_col = "pxl_row_in_fullres",
                              y_col = "pxl_col_in_fullres",
                              min_value = 1, max_value = 10000) {

  # Check inputs
  if(!all(c(sample_col, x_col, y_col) %in% names(data))) {
    stop("One or more specified columns not found in data")
  }

  # Perform normalization
  normalized_data <- data %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::mutate(
      !!x_col := (.data[[x_col]] - min(.data[[x_col]])) /
        (max(.data[[x_col]]) - min(.data[[x_col]])) *
        (max_value - min_value) + min_value,
      !!y_col := (.data[[y_col]] - min(.data[[y_col]])) /
        (max(.data[[y_col]]) - min(.data[[y_col]])) *
        (max_value - min_value) + min_value
    ) %>%
    dplyr::ungroup()

  return(normalized_data)
}
