#' Perform correlation analysis between spatial features and distance metrics with visualization
#'
#' @param spatial_data Spatial data containing feature columns and Newbarcode identifier
#' @param distance_results Distance results containing distance metrics and Newbarcode identifier
#' @param spatial_feature Column name from spatial_data to use for correlation (e.g., "gen2_SPLIz_numeric")
#' @param distance_metric Column name from distance_results to use for correlation (e.g., "Epithelial_cells_A")
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param plot Logical, whether to generate a scatter plot
#' @param plot_title Title for the scatter plot (optional)
#' @return A list containing correlation results and ggplot object (if plot=TRUE)
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom Hmisc rcorr
#' @importFrom stats cor.test
#' @examples
#' calculate_correlations(spatial_data = posi,
#'              distance_results = distance_results,
#'              spatial_feature = "gen2_SPLIz_numeric",
#'              distance_metric = "Epithelial_cells_A",
#'              method = "pearson",
#'              plot = TRUE,
#'              plot_title = "Correlation between Gene Expression and Distance")

calculate_correlations <- function(spatial_data,
                                   distance_results,
                                   spatial_feature,
                                   distance_metric,
                                   method = "pearson",
                                   plot = TRUE,
                                   plot_title = NULL) {

  # Input validation
  if (!spatial_feature %in% colnames(spatial_data)) {
    stop(paste("Spatial feature", spatial_feature, "not found in spatial_data"))
  }

  if (!distance_metric %in% colnames(distance_results)) {
    stop(paste("Distance metric", distance_metric, "not found in distance_results"))
  }

  # Merge data by Newbarcode
  merged_data <- merge(spatial_data, distance_results,
                       by = "Newbarcode", all.x = FALSE)

  # Extract the columns of interest
  x <- merged_data[[spatial_feature]]
  y <- merged_data[[distance_metric]]

  # Calculate correlation
  cor_test <- cor.test(x, y, method = method)
  hrcor <- Hmisc::rcorr(cbind(x, y), type = method)

  # Prepare results
  results <- list(
    estimate = cor_test$estimate,
    p_value = cor_test$p.value,
    method = method,
    n_observations = length(x),
    spatial_feature = spatial_feature,
    distance_metric = distance_metric
  )

  # Generate plot if requested
  if (plot) {
    if (is.null(plot_title)) {
      plot_title <- paste("Correlation between", spatial_feature, "and", distance_metric)
    }

    # Create scatter plot with correlation annotation
    plot_obj <- ggplot(merged_data, aes(x = .data[[spatial_feature]], y = .data[[distance_metric]])) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred") +
      labs(title = plot_title,
           x = spatial_feature,
           y = distance_metric) +
      theme_minimal() +
      annotate("text",
               x = min(x, na.rm = TRUE) + 0.1*diff(range(x, na.rm = TRUE)),
               y = max(y, na.rm = TRUE) - 0.1*diff(range(y, na.rm = TRUE)),
               label = paste0(method, " r = ", round(results$estimate, 3),
                              "\np = ", format.pval(results$p_value, digits = 2)),
               hjust = 0, size = 5)

    results$plot <- plot_obj
  }

  return(results)
}
