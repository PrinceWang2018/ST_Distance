\name{distance_results}
\alias{distance_results}
\docType{data}
\title{
distance_results
}
\description{
A data frame with nearest distances for each reference cell
}
\usage{data("distance_results")}
\format{
  A data frame with 18 observations on the following 5 variables.
  \describe{
    \item{\code{Newbarcode}}{a character vector}
    \item{\code{Epithelial_cells_A}}{a numeric vector}
    \item{\code{Epithelial_cells_B}}{a numeric vector}
    \item{\code{Epithelial_cells_C}}{a numeric vector}
    \item{\code{Epithelial_cells_D}}{a numeric vector}
  }
}
\details{
A data frame with nearest distances for each reference cell.
}
\source{
This study.
}
\references{
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
}
\examples{
data(distance_results)
## maybe str(distance_results) ; plot(distance_results) ...
}
\keyword{datasets}
