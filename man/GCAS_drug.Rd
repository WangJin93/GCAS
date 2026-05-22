\name{GCAS_drug}
\alias{GCAS_drug}
\docType{data}
\title{Anti-tumor drug sensitivity of all samples}
\description{Anti-tumor drug sensitivity of all samples calculated based on oncoPredict package and GDSC database.}
\usage{data("GCAS_drug")}
\format{
  A data frame containing drug sensitivity data for cancer samples.
  \describe{
    \item{\code{ID}}{Sample ID}
    \item{\code{Camptothecin_1003}}{Drug sensitivity score}
    \item{\code{Vinblastine_1004}}{Drug sensitivity score}
    \item{\code{Cisplatin_1005}}{Drug sensitivity score}
    \item{\code{Docetaxel_1007}}{Drug sensitivity score}
    \item{\code{Gefitinib_1010}}{Drug sensitivity score}
  }
}
\references{
Garnett, M. J., et al. (2012). Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature, 483(7391), 570-575.
}
\examples{
data(GCAS_drug)
str(GCAS_drug)
}