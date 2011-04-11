\name{MethIll}
\alias{MethIll}
\alias{summary.MethIll}
\alias{print.MethIll}
\alias{plot.MethIll}
\alias{MethIll-package}
\title{Main wrapper method of Methill package}
\description{ Accepts a MethyLumiSet and can be called from to perform batch effect cancelation, sample network
analyses.  }
\usage{ MethIll(x, batch = NULL, covariate = NULL, strata = NULL, asfactor = NULL, filter = NULL, latex = FALSE,
impute = TRUE) }
\arguments{
  \item{x}{ a methylumiset object }
  \item{batch}{ batches whose effect should be canceled out from the assay data}
  \item{covariate}{ covariates whose signal should be preserved while canceling batch effects with ComBat}
  \item{strata}{ not use }
  \item{asfactor}{ which variables are to be considered as factor }
  \item{filter}{ whether to filter out samples having scaled intra array connectivity lower than a provided threshold (\code{filter})}
  \item{latex}{ whether to report output as latex }
  \item{impute}{ whether to perform imputation, by default TRUE}
}
\author{ Fabrice Colas }
