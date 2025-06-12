.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("=========================================================================================
", pkgname, " version ", version, "
Project URL: https://github.com/PrinceWang2018/STDistance
If you use it in published research, please cite:
Wang, Z., Yang, L., Yang, S., Li, G., Xu, M., Kong, B., Shao, C., & Liu, Z. (2025).
Isoform switch of CD47 provokes macrophage-mediated pyroptosis in ovarian cancer.
bioRxiv, 2025.2004.2017.649282. https://doi.org/10.1101/2025.04.17.649282
=========================================================================================
                       --A joyful heart is good medicine! ^_^--")
  base::packageStartupMessage(msg)
}
