\name{tissue_posi}
\alias{tissue_posi}
\docType{data}
\title{
tissue_posi
}
\description{
Position of the spots from spatial transcriptome data.
}
\usage{data("tissue_posi")}
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
Please find the tissue_positions.csv from the spaceranger output files: ./out/spatial/tissue_positions.csv
Multiple samples should be merged together, adding "samples" and "Newbarcode" columns.
}
\source{
This study.
}
\references{
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025).
Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
}
\examples{
data(tissue_posi)
## maybe str(tissue_posi) ; plot(tissue_posi) ...
}
\keyword{datasets}
