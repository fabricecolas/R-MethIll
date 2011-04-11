\name{MLStack}
\alias{MLStack}
\alias{MLearnExp}
\alias{print.MLearnExp}
\title{Create a stack of machine learning experiments}
\description{ The stack of machine learning experiments to a wrapper function of MLearn}  
\usage{ MLStack(Class, DataSet, FeatSel = "rand", NFeat = as.character(as.integer(2^(0.5 + 3:9))), Algo =
   c("svmI"), DataHyp = c("data", "H0"), DataH0Seed = 6013:6017) }
\arguments{
  \item{Class}{ the different levels of the two or more target classes }
  \item{DataSet}{ dataset to perform experiment on, e.g. meth, unmethylated, beta }
  \item{FeatSel}{ feature selection method to call }
  \item{NFeat}{ number of features to be selected }
  \item{Algo}{ MLInt classifiers }
  \item{DataHyp}{ data hypotheses to be tested }
  \item{DataH0Seed}{ if H0, which random seeds}
}
\author{ Fabrice Colas }
