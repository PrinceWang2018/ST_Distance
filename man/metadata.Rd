\name{metadata}
\alias{metadata}
\docType{data}
\title{
Metadata
}
\description{
Metadata of a Seurat Object.
}
\usage{data("metadata")}
\format{
  A data frame with 293 observations on the following 9 variables.
  \describe{
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
Metadata dataframe can be exported from Seurat Object as follows. The Seurat Object should be fully annotated and the metadata dataframe should included:
orig.ident, celltype, the target gene expression, splicing index (SpliZ), etc.
}
\source{
This study.
}
\references{
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025). Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
}
\examples{
data(metadata)
## maybe str(metadata) ; plot(metadata) ...
}
\keyword{datasets}
