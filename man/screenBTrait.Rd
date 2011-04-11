\name{screenBTrait}
\alias{screenBTrait}
\alias{plot.screenBTrait}
\alias{print.screenBTrait}
\alias{summary.screenBTrait}
\title{Standard screening binary trait (WGCNA)}
\description{ Association testing over the whole set of probes }
\usage{ screenBTrait(x, btrait = "SampleLabel", nnodes = 1, csv = TRUE) }
\arguments{
  \item{x}{ a \code{MethIll} object}
  \item{btrait}{ name of the binary trait covariate}
  \item{nnodes}{ number of computing nodes for parallel computing by snow }
  \item{csv}{ whether to save as csv the output}
}
\references{For WGCNA: \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/} }
\author{Steve Horvath et al. (WGCNA), Fabrice Colas (packaging, adpatation and extension) }
