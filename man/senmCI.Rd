\name{senmCI}
\alias{senmCI}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Confidence interval}
\title{
Sensitivity Analysis for a Confidence Interval in a Full Match.
}
\description{
Each matched set contains one treated
individual and one or more controls.
Uses Huber's M-statistic as the basis for
a confidence interval for an additive constant treatment effect,
\eqn{\tau}.  Produces either a randomization
based confidence interval or an analysis of sensitivity to
departures from random assignment.  Also produces a point estimate
for randomization inference or an interval of point estimates
for a sensitivity analysis.  For tests, use function senm().
The method is described in Rosenbaum (2007,2013).
}
\usage{
senmCI(y, z, mset, gamma=1, inner=0, trim=3, lambda=1/2,
                  alpha=0.05, twosided=TRUE, upper=TRUE, TonT=FALSE)
}
\arguments{
\item{y}{ A vector of responses with no missing data.}
  \item{z}{ Treatment indicators, z=1 for treated, z=0 for control with
  length(z)==length(y).
}
  \item{mset}{ Matched set indicators, 1, 2, ..., sum(z) with
  length(mset)==length(y).  Matched set indicators should be either integers or a factor.
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming ignorable treatment assignment given
the matched sets, and it returns a randomization-based confidence interval.
}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
The default values yield a version of Huber's \eqn{\psi}-function, while
setting inner = 0 and trim = Inf uses the mean within each matched set.
The \eqn{\psi}-function is an odd function, so \eqn{\psi(w) = -\psi(-w)}.  For
\eqn{w \ge 0}, the \eqn{\psi}-function is \eqn{\psi(w)=0} for \eqn{0 \le w \le } inner,
is \eqn{\psi(w)= } trim for \eqn{w \ge } trim, and rises linearly from 0 to trim
for inner < w < trim.

An error will result unless \eqn{0 \le } inner \eqn{ \le } trim.

Taking trim < Inf limits the influence of outliers; see Huber (1981).
Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
See inner.
}
  \item{lambda}{
Before applying the \eqn{\psi}-function to treated-minus-control differences,
the differences are scaled by dividing by the lambda quantile of all
within set absolute differences.  Typically, lambda = 1/2 for the median.
The value of lambda has no effect if trim=Inf and inner=0.  See Maritz (1979)
for the paired case and Rosenbaum (2007) for matched sets.

An error will result unless 0 < lambda < 1.
}
  \item{alpha}{
The coverage rate of the confidence interval is 1-alpha. If the bias in
treatment assignment is at most \eqn{\Gamma}, then the confidence interval
will cover the true \eqn{\tau} with probability at least \eqn{1-\alpha}.
}
  \item{twosided}{
If twosided==TRUE, then a two-sided \eqn{1-\alpha} confidence interval is
constructed.  If twosided==FALSE, then a one-sided \eqn{1-\alpha}
confidence interval is constructed.  The two sided interval is the
intersection of two one-sided \eqn{1-\alpha/2} intervals.
  }
  \item{upper}{
If twosided==TRUE, then upper is ignored. If twosided==FALSE and
upper=TRUE, then an upper \eqn{1-\alpha} confidence interval
is constructed.  If twosided==FALSE and upper=FALSE, then
a lower \eqn{1-\alpha} confidence interval is constructed.
  }
\item{TonT}{
TonT refers to the effect of the treatment on the treated; see Rosenbaum and Rubin (1985, equation 1.1.1)  The default is TonT=FALSE.  If TonT=FALSE, then the total score in matched set i is divided by the number ni of individuals in set i, as in expression (8) in Rosenbaum (2007).  This division by ni has few consequences when every matched set has the same number of individuals, but when set sizes vary, dividing by ni is intended to increase efficiency by weighting inversely as the variance; see the discussion in section 4.2 of Rosenbaum (2007).  If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets to make it a type of mean.  If TonT=TRUE and trim=Inf, then the statistic is the mean over matched sets of the treated minus mean-control response, so it is weighted to estimate the average effect of the treatment on the treated.
}
}
\details{
For the given \eqn{\Gamma}, senmCI() inverts the test in the function
senm() to produce the confidence interval.  That is, it tests every
\eqn{\tau} and retains the values not rejected at level \eqn{\alpha}.

The test is a large sample approximation
based on a Normal approximation to the null distribution;
see Rosenbaum (2007).

If TonT=FALSE, matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).
}
\value{
\item{PointEstimates }{The interval of point estimates of \eqn{\tau}.
If gamma=1, then the interval is a single point estimate.}
\item{ConfidenceInterval }{The confidence interval for \eqn{\tau}.}
\item{description }{Reminder of the coverage rate, gamma, and
type of interval.}
}


\references{

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (1993). Hodges-Lehmann point estimates of
treatment effect in observational studies.
Journal of the American Statistical Association, 88, 1250-1253.
(Introduces sensitivity analysis for point estimates.)

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv) <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. DOI: 10.1214/16-AOAS942.

Rosenbaum, P. R., & Rubin, D. B. (1985). The bias due to incomplete
matching. Biometrics, 41, 103-116.

}
\author{
Paul R. Rosenbaum.
}
\note{
In a sensitivity analysis, a one-sided confidence interval is
not conservative; however, two-sided intervals formed
as the intersection of two one-sided \eqn{1-\alpha/2}
intervals are somewhat conservative.  See the
discussion of two-sided tests in the documentation
for senm().
}


\examples{
data(teeth)
attach(teeth)
#
# Note: Computing confidence intervals takes a few moments
# The calls to senmCI() are commented to meet time requirements
# for cran examples.  Remove the comment symbol to run them.
#
# The calculations that follow reproduce the intervals from
# section 5.1 of Rosenbaum (2016, p. 1466)
# senmCI(either4low,smoker,mset,trim=2.5,gamma=1.5)
# senmCI(either4up,smoker,mset,trim=2.5,gamma=1.5)
# Confidence interval using the mean by inverting the
# permuational t-test.
# senmCI(either4low,smoker,mset,trim=Inf,TonT=TRUE)
dif<-either4low[smoker==1]-either4low[smoker==0] # Matched pair differences
mean(dif) # Equals the point estimate above
t.test(dif) # But permutational t-interval and t-interval differ
# Sensitivity analysis using the mean difference
# senmCI(either4low,smoker,mset,gamma=1.5,trim=Inf,TonT=TRUE)
detach(teeth)
}

\keyword{ htest }
\keyword{ robust }
