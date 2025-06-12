\name{posi}
\alias{posi}
\docType{data}
\title{
posi
}
\description{
Merged spatial location and metadata information.
}
\usage{data("posi")}
\format{
  A data frame with 293 observations on the following 18 variables.
  \describe{
    \item{\code{Newbarcode}}{a character vector}
    \item{\code{barcode}}{a character vector}
    \item{\code{in_tissue}}{a numeric vector}
    \item{\code{array_row}}{a numeric vector}
    \item{\code{array_col}}{a numeric vector}
    \item{\code{pxl_row_in_fullres}}{a numeric vector}
    \item{\code{pxl_col_in_fullres}}{a numeric vector}
    \item{\code{Sample}}{a character vector}
    \item{\code{Sampleid}}{a numeric vector}
    \item{\code{orig.ident}}{a character vector}
    \item{\code{nCount_Spatial}}{a numeric vector}
    \item{\code{nFeature_Spatial}}{a numeric vector}
    \item{\code{nCount_SCT}}{a numeric vector}
    \item{\code{nFeature_SCT}}{a numeric vector}
    \item{\code{integrated_snn_res.0.8}}{a numeric vector}
    \item{\code{seurat_clusters}}{a numeric vector}
    \item{\code{celltype_ABCDepi}}{a character vector}
    \item{\code{gen2_SPLIz_numeric}}{a numeric vector}
  }
}
\details{
Merged spatial location and metadata information.
}
\source{
This study.
}
\references{
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
}
\examples{
data(posi)
## maybe str(posi) ; plot(posi) ...
}
\keyword{datasets}
