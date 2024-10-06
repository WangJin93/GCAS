\name{drug_info}
\alias{drug_info}
\docType{data}
\title{
Drug informations
}
\description{
Anti-tumor drugs informations obtained from GDSC2.0 datasets.
}
\usage{data("drug_info")}
\format{
  A data frame with 198 observations on the following 6 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Name}}{a character vector}
    \item{\code{Synonyms}}{a character vector}
    \item{\code{Targets}}{a character vector}
    \item{\code{Target.pathway}}{a character vector}
    \item{\code{PubCHEM}}{a character vector}
  }
}
\source{
Genomics of Drug Sensitivity in Cancer (https://www.cancerrxgene.org/)
}
\examples{
data(drug_info)
str(drug_info)
}
