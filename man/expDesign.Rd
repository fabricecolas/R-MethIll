\name{expDesign}
\alias{expDesign}
\title{ Characterize the experimental design of an analysis}
\description{ The two by two joint distribution of the factors provided as parameters to the method are reported.
This may enable to identify experimental design artifacts, such as non-randomization. Indeed, batch effect removal
is not possible when batches and covariates to preserve are singularly associated. }
\usage{ expDesign(x, gr1 = NULL, gr2 = NULL, latex = FALSE, exclude = NULL) }
\arguments{
  \item{x}{ a matrix or data frame }
  \item{gr1}{ eventually null, a first group of factors to compare; either one to another when gr2 is null, or with
  respect to gr2}
  \item{gr2}{ a second group of factors to compare against gr1 }
  \item{latex}{ whether the output should be reported as latex }
  \item{exclude}{ whether missing values should be included as levels of the joint distribution; by default TRUE}
}
\author{ Fabrice Colas}
