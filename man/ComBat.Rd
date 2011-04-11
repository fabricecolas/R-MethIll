\name{ComBat}
\alias{ComBat}
\alias{plot.ComBat}
\title{Cancels batch effects by combining batches of gene expression microarray data}
\description{ Non-biological experimental variation or batch effects are commonly observed across multiple batches
of microarray experiments, often rendering the task of combining data from these batches difficult. The ability to
combine microarray data sets is advantageous to researchers to increase statistical power to detect biological
phenomena from studies where logistical considerations restrict sample size or in studies that require the
sequential hybridization of arrays. In general, it is inappropriate to combine data sets without adjusting for
batch effects.  Methods have been proposed to filter batch effects from data, but these are often complicated and
require large batch sizes (>25) to implement. Because the majority of microarray studies are conducted using much
smaller sample sizes, existing methods are not sufficient. We propose parametric and nonparametric empirical Bayes
frameworks for adjusting data for batch effects that is robust to outliers in small sample sizes and performs
comparable to existing methods for large samples. We illustrate our methods using two example data sets and show
that our methods are justifiable, easy to apply, and useful in practice.  }
\usage{ ComBat(x, pdat, batch, covariate, par.prior = TRUE) }
\arguments{
  \item{x}{ the expression data }
  \item{pdat}{ data frame of the pheno data describing each batch }
  \item{batch}{ the series of batches to adjust for }
  \item{covariate}{ the signal to preserve while canceling batch effects }
  \item{par.prior}{ type of prior probabilities (cf. code)}
}
\references{ Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data
using Empirical Bayes methods. Biostatistics 8(1):118-127. 

\url{http://jlab.byu.edu/ComBat/Abstract.html} }
\author{ Johnson WE, Rabinovic A, and Li C}
