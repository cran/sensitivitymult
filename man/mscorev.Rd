\name{mscorev}
\alias{mscorev}

\title{
Computes M-scores for M-tests and estimates.
}
\description{
Of limited interest to most users, function mscorev() computes the M-scores used by
functions senm(), senmCI(), comparison(), and principal() that perform Huber-Maritz M-tests.
The function is also in the package sensitivitymv.
}
\usage{
mscorev(ymat, inner = 0, trim = 2.5, qu = 0.5, TonT = FALSE)
}

\arguments{
  \item{ymat}{
ymat is a matrix as described in the documentation for senm().
}
  \item{inner}{
inner is the parameter described in the documentation for senm().
The inner parameter is discussed in Rosenbaum (2013).
}
  \item{trim}{
trim is the parameter described in the documentation for senm().  Note that the default
for mscorev is trim = 2.5, but other functions in this package reset the default to
trim = 3.  This is for consistency with the sensitivitymv package which also contains
the mscorev function.
}
  \item{qu}{
qu is the lambda parameter described in the documentation for senm().
}
\item{TonT}{
If TonT=FALSE, then the total score in set (row) i is divided by the number ni of individuals in row i, as in expression (8) in Rosenbaum (2007).  If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets.  See the discussion of TonT in the documentation for senm().
}
}

\value{
Generally, a matrix with the same dimensions as ymat containing the M-scores.  Exception: if a matched set does not contain at least one treated subject and at least one control, then that set will not appear in the result, and the result will have fewer rows than ymat.  However, if a matched set has several controls but no treated subject, then these controls will contribute to the estimate of the scale parameter, typically the median absolute pair difference.
}
\references{
Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464.  <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)
}
\author{
Paul R. Rosenbaum
}

\note{
The example reproduces Table 3 in Rosenbaum (2007).
}

\examples{
# The example reproduces Table 3 in Rosenbaum (2007).
data(tbmetaphase)
mscorev(tbmetaphase,trim=1)
}
