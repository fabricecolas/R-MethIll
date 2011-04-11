\name{confuMat2}
\alias{confuMat2}
\alias{summary.confuMat2}
\title{Confusion matrix from a machine learning experiment}
\description{ConfuMat2 is an extension of confuMat from MLInterfaces, that returns a confusion matrix with
all class-levels, even those rows or columns with zero counts and which calculates the precision, recall and F1
measure for class. The macroF1 is also calculated, as an unweighted mean of the F1 scores for each class. }
\usage{ confuMat2(x, clSet = NULL) }
\arguments{
  \item{x}{ an MLInt classifier (in the same way as MLInterfaces)}
  \item{clSet}{ the levels of the class factor }
}
\author{ Fabrie Colas}
