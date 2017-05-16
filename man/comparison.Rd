\name{comparison}
\alias{comparison}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Sheffe projection}
\title{
Sensitivity Analysis for a Comparison Involving Several Outcomes in an Observational Study.
}
\description{
For multiple outcomes in an observation study, computes
a weighted combination of M-statistics, one for each outcome,
and performs either a one-sided randomization test or an analysis
of sensitivity to departures from random assignment.
Each matched set contains one treated individual and one or more controls.
The method is described in Rosenbaum (2016).  For one outcome,
use the function senm().
}
\usage{
comparison(y, z, mset, w, gamma = 1, inner = 0, trim = 3, lambda = 0.5,
     TonT = FALSE, apriori = FALSE, Scheffe = FALSE)
}
\arguments{
  \item{y}{ A matrix of responses with no missing data.  Different
  columns of y are different variables.  If present, the column names
  of y are used to label output.}
  \item{z}{ Treatment indicators, z=1 for treated, z=0 for control with
  length(z)==dim(y)[1].
}
  \item{mset}{ Matched set indicators, 1, 2, ..., sum(z) with
  length(mset)==dim(y)[1].  The vector mset may contain integers
  or may be a factor.
}
  \item{w}{ Vector of weights to be applied to the M-statistics
  for the several outcomes with length(w)==dim(y)[2].  At least
  one weight must be nonzero.  See Details for discussion of
  scaling.
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming ignorable treatment assignment given
the matched sets, and it performs a within-set randomization test.
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
Taking trim < Inf and inner = 0 uses Huber's psi function.
Taking trim = Inf does no trimming and is similar to a weighted mean; see TonT.
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
\item{TonT}{
TonT refers to the effect of the treatment on the treated; see Rosenbaum and Rubin (1985, equation 1.1.1)  The default is TonT=FALSE.  If TonT=FALSE, then the total score in matched set i is divided by the number ni of individuals in set i, as in expression (8) in Rosenbaum (2007).  This division by ni has few consequences when every matched set has the same number of individuals, but when set sizes vary, dividing by ni is intended to increase efficiency by weighting inversely as the variance; see the discussion in section 4.2 of Rosenbaum (2007).  If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets to make it a type of mean.  If TonT=TRUE and trim=Inf, then the statistic is the mean over matched sets of the treated minus mean-control response, so it is weighted to estimate the average effect of the treatment on the treated.
}
  \item{apriori}{
  If Scheffe=FALSE and apriori=TRUE, then the weights are assumed to have been chosen
  a priori, and a one-sided, uncorrected P-value is reported for gamma=1 or an
  upper bound on the one-sided, uncorrected P-value is reported for gamma>1. In
  either case, this is a Normal approximation based on the central limit
  theorem and equals 1-pnorm(deviate).
}
  \item{Scheffe}{
  If Scheffe=TRUE, then the weights are assumed to have been chosen after
  looking at the data.  In this case, the P-value or P-value bound is
  corrected using Scheffe projections.  The approximate corrected P-value
  or P-value bound is 1-pchisq(max(0,deviate)^2,dim(y)[2]).  If Scheffe=FALSE
  and apriori=FALSE, then the deviate is returned, but no P-value is given.
  See Rosenbaum (2016).  A Scheffe correction entitles you to look in
  both tails, which you do by considering both w and -w.
  See the planScheffe() function for a combination
  of an apriori and Scheffe comparisons.
}
}
\details{
If y has k columns for k outcomes, then comparison computes k M-statistics,
one for each outcome, combines them into a single statistic using weights
w, and computes a one-sided, upper-tailed deviate for a randomization test
or a sensitivity analysis, as described in Rosenbaum (2016).  The k individual
statistics for the k outcomes separately are as described in Rosenbaum (2007).

When trim<Inf, outcomes are scaled using by the \eqn{\lambda} quantile of
the absolute differences before applying the \eqn{\psi}-function.  In this
sense, when trim<Inf, the test statistics share a common scaling before they
are combined using the weights in w.  Weights w=c(1,1) give equal emphasis
to two outcomes that may be recorded in different units, inches or pounds
or whatever.

When inner=0 and trim=Inf, the \eqn{\psi}-function is the identity, and no
scaling is done.  In this case, the weights refer to the original variables
in their unscaled original units, inches or pounds or whatever.  Because of
this, the weights w have a different meaning with trim=Inf or trim<Inf.
In the comparison() function, inner>0 is not permitted if trim=Inf.

If one has a single a priori choice of weights, w, then the one-sided
P-value (for gamma=1) or the one-sided upper bound on the P-value
(for gamma>1) is approximately 1-pnorm(deviate).

If one considers every possible choice of weights, w, then the
P-value (for gamma=1) or the upper bound on the P-value
(for gamma>1) uses a Scheffe correction and is approximately
1-pchisq(max(0,deviate^2),dim(y)[2]); see Rosenbaum (2016).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).
}
\value{
\item{deviate }{The upper bound on the standardized deviate
that is used to approximate P-values using the Normal or
chi-square distribution; see apriori and Scheffe in the
arguments.}
\item{aprioriPval }{If Scheffe=FALSE and apriori=TRUE, the
deviate is compared to the upper tail of the Normal
distribution to produce either a P-value for gamma=1
or an upper bound on the P-value for gamma>1.}
\item{ScheffePval }{If Scheffe=TRUE, the
deviate is compared to the upper tail of the
chi-square distribution to produce either
a P-value for gamma=1
or an upper bound on the P-value for gamma>1.}
\item{weights }{The weights possibly relabeled with
colnames of y.}
}
\references{

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv) <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming, inner>0.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Rosenbaum, P. R., & Rubin, D. B. (1985). The bias due to incomplete
matching. Biometrics, 41, 103-116.

Scheffe, H. (1953) A method for judging all contrasts in the analysis
of variance.  Biometrika, 40, 87-104.

}
\author{
Paul R. Rosenbaum.
}
\note{
For confidence intervals for individual outcomes, use function senmCI().
}


\examples{
data(teeth)
attach(teeth)
# The following calculation reproduces the comparison in
# expression (5.1) of Rosenbaum (2016, p. 1464)
comparison(cbind(teeth$either4low,teeth$either4up),teeth$smoker,
   teeth$mset,c(.714,.286),gamma=2.2,trim=2.5,Scheffe=TRUE)
#
# The following example reproduces the deviate for lower teeth
# mentioned on line 4 of Rosenbaum (2016, p. 1466) as a one-sided
# test with w picked a priori as w=c(1,0):
comparison(cbind(either4low,either4up),smoker,mset,c(1,0),
   trim=2.5,gamma=2.2,apriori=TRUE)
# Because the previous comparison implicitly involves just one outcome, it
# could be done more simply with senm() as follows:
senm(either4low,smoker,mset,trim=2.5,gamma=2.2)
# Had one done all comparisons including the comparison for lower teeth,
# then one would need to adjust for multiple testing:
comparison(cbind(either4low,either4up),smoker,mset,c(1,0),
   trim=2.5,gamma=2.2,Scheffe=TRUE)
detach(teeth)
}

\keyword{ htest }
\keyword{ robust }
