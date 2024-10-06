\name{sample_subtype}
\alias{sample_subtype}
\docType{data}
\title{
Subtype information of all samples
}
\description{
Subtype information of all samples
}
\usage{data("sample_subtype")}
\format{
  A data frame with 28032 observations on the following 5 variables.
  \describe{
    \item{\code{tissue}}{a character vector}
    \item{\code{ID}}{a character vector}
    \item{\code{dataset}}{a character vector}
    \item{\code{subtype}}{a character vector}
    \item{\code{Patient.ID}}{a character vector}
  }
}
\examples{
data(sample_subtype)
str(sample_subtype)
}
