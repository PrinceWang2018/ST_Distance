\name{tissue_posi_normalized}
\alias{tissue_posi_normalized}
\docType{data}
\title{
tissue_posi_normalized
}
\description{
Normalized tissue spots position.
}
\usage{data("tissue_posi_normalized")}
\format{
  A data frame with 9984 observations on the following 9 variables.
  \describe{
    \item{\code{barcode}}{a character vector}
    \item{\code{in_tissue}}{a numeric vector}
    \item{\code{array_row}}{a numeric vector}
    \item{\code{array_col}}{a numeric vector}
    \item{\code{pxl_row_in_fullres}}{a numeric vector}
    \item{\code{pxl_col_in_fullres}}{a numeric vector}
    \item{\code{Sample}}{a character vector}
    \item{\code{Sampleid}}{a numeric vector}
    \item{\code{Newbarcode}}{a character vector}
  }
}
\details{
Normalized tissue spots position.
}
\source{
This study.
}
\references{
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
}
\examples{
data(tissue_posi_normalized)
## maybe str(tissue_posi_normalized) ; plot(tissue_posi_normalized) ...
}
\keyword{datasets}
