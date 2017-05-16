\name{tbmetaphase}
\alias{tbmetaphase}
\docType{data}
\title{
Genetic damage from drugs used to treat TB
}
\description{
This is a matched comparison of the effects of two drug sequences, namely HRZ and H2R2Z2, for the treatment of tuberculosis.  HRZ is a higher dose sequence than H2R2Z2.  The outcome is a measure of genetic damage, namely the frequency of aberrant metaphases two months after treatment.  Individuals were matched for the frequency of aberrant metaphases before treatment.  15 individuals treated with HRZ are matched to 1 or 2 controls treated with H2R2Z2.  Each row is one matched set.  If a set is a pair, the third element in a row is NA.  The data are originally from Rao, Gupta and Thomas (1991) and were used as an example in Rosenbaum (2007, Table 3).  Data are used to illustrate the senmv function in the sensitivitymv package.
}
\usage{data(tbmetaphase)}
\format{
  A data frame with 15 observations on the following 3 variables.
  \describe{
    \item{\code{HRZ}}{Aberrant metaphases for individual treated with HRZ.}
    \item{\code{H2R2Z2.1}}{Aberrant metaphases for first matched individual treated with H2R2Z2.}
    \item{\code{H2R2Z2.2}}{Aberrant metaphases for second matched individual treated with H2R2Z2.  For matched pairs, this is NA.}
  }
}

\references{
Rao, V. V. N. G., Gupta, E. V. V and Thomas, I. M. Chromosomal aberrations in tuberculosis patients before and after treatment with short-term chemotherapy.  Mutation Research 1991, 259, 13-19.

Rosenbaum, P. R. Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464.
}
\examples{
data(tbmetaphase)
}
\keyword{datasets}
