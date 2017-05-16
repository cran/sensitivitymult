\name{senm}
\alias{senm}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matching}
\concept{rbounds}
\title{
Sensitivity Analysis for a Matched Comparison in an Observational Study.
}
\description{
Each matched set contains one treated
individual and one or more controls.
Uses Huber's M-statistic as the basis for
the test, for instance, a mean.  Performs either a randomization
test or an analysis of sensitivity to departures from random
assignment.  For confidence intervals, use function senmCI().
The method is described in Rosenbaum (2007,2013).
The senm() function is intended as a convenience for a user
of the comparison() function in the sensitivitymult package.
The function senm() in the sensitivitymult package is essentially
the same as senmv() in the sensitivitymv package, except the
format of the input to senm() resembles the format of the
input to comparison().  In particular, in the sensitivitymv
package, the rows of y are matched sets, not people, whereas
in sensitivitymult the rows of y are people with treated
and control people identified by z and matched sets
identified by mset.
}
\usage{
senm(y, z, mset, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
              tau = 0, alternative="greater", TonT = FALSE)
}
\arguments{
  \item{y}{ A vector of responses with no missing data.}
  \item{z}{ Treatment indicator, z=1 for treated, z=0 for control with
  length(z)==length(y).
}
  \item{mset}{ Matched set indicator, 1, 2, ..., sum(z) with
  length(mset)==length(y).  Matched set indicators should be either integers or a factor.
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
  \item{tau}{
The null hypothesis asserts that the treatment has an additive effect, tau.
By default, tau=0, so by default the null hypothesis is Fisher's sharp null
hypothesis of no treatment effect.
}
  \item{alternative}{
If alternative="greater", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect larger than tau.  If
alternative="less", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect smaller than tau.
In particular, alternative="less" is equivalent to: (i) alternative="greater",
(ii) y replaced by -y, and (iii) tau replaced by -tau.  See the note for
discussion of two-sided sensitivity analyses.
  }
\item{TonT}{
TonT refers to the effect of the treatment on the treated; see Rosenbaum and Rubin (1985, equation 1.1.1)  The default is TonT=FALSE.  If TonT=FALSE, then the total score in matched set i is divided by the number ni of individuals in set i, as in expression (8) in Rosenbaum (2007).  This division by ni has few consequences when every matched set has the same number of individuals, but when set sizes vary, dividing by ni is intended to increase efficiency by weighting inversely as the variance; see the discussion in section
4.2 of Rosenbaum (2007).  If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets to make it a type of mean.  If TonT=TRUE and trim=Inf, then the statistic is the mean over matched sets of the treated minus mean-control response, so it is weighted to estimate the average effect of the treatment on the treated.  See the examples.
}
}
\details{
For the given \eqn{\Gamma}, senm() computes the upper bound on the 1-sided
P-value testing the null hypothesis
of an additive treatment effect tau against the alternative hypothesis of
a treatment effect larger than tau.  By default, senm() tests the null hypothesis of
no treatment effect against the alternative of a positive treatment effect.
The P-value is an approximate P-value
based on a Normal approximation to the null distribution; see Rosenbaum (2007).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).
}
\value{
\item{pval }{The upper bound on the 1-sided P-value.}
\item{deviate }{The deviate that was compared to the Normal
distribution to produce pval.}
\item{statistic }{The value of the M-statistic.}
\item{expectation }{The maximum expectation of the
M-statistic for the given \eqn{\Gamma}.}
\item{variance }{The maximum variance of the
M-statistic among treatment assignments
that achieve the maximum expectation.
Part of the separable approximation.}
}
\references{
Cox, D. R. (1977).  The role of signficance tests (with Discussion).
Scand. J. Statist. 4, 49-70.

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
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2014). Weighted M-statistics with superior design
sensitivity in matched observational studies with multiple controls.
J. Am. Statist. Assoc. 109 1145-1158. (R package sensitivitymw)
<doi:10.1080/01621459.2013.879261>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Rosenbaum, P. R., & Rubin, D. B. (1985). The bias due to incomplete
matching. Biometrics, 41, 103-116.

}
\author{
Paul R. Rosenbaum.
}
\note{
The function senm() performs 1-sided tests.  One approach
to a 2-sided, \eqn{\alpha}-level test does both 1-sided tests
at level \eqn{\alpha/2}, and rejects the null hypothesis if either
1-sided
test rejects.  Equivalently, a bound on the two sided
P-value is the smaller of 1 and twice the smaller of the two 1-sided
P-values.  This approach views a 2-sided test as two 1-sided tests
with a Bonferroni correction; see Cox (1977, Section 4.2).  In all
cases, this approach is a valid large sample test: a true
null hypothesis is falsely
rejected with probability at most \eqn{\alpha} if the bias in
treatment assignment is at most \eqn{\Gamma}; so, this procedure
is entirely safe to use.  For a randomization test, \eqn{\Gamma=1}, this
Bonferroni procedure is not typically conservative. For large \eqn{\Gamma},
this Bonferroni procedure tends to be somewhat conservative.

Related packages are sensitivitymv, sensitivitymw, sensitivityfull
and sensitivity2x2xk.
}


\examples{
data(teeth)
attach(teeth)
# The following example reproduces the deviate for lower teeth
# mentioned on line 4 of Rosenbaum (2016, p. 1466).
senm(either4low,smoker,mset,trim=2.5,gamma=2.2)
# The calculation above is equivalent to using comparison()
# with weights w=c(0,1) so upper teeth are ignored.
comparison(cbind(either4up,either4low),smoker,mset,c(0,1),trim=2.5,gamma=2.2)
# The following example illustrates the permutational t-test
# which uses the mean of the pair differences as a test statistic.
senm(either4low,smoker,mset,trim=Inf,TonT=TRUE)
dif<-either4low[smoker==1]-either4low[smoker==0] # Matched pair differences
mean(dif) # Equals the test statistic above
detach(teeth)
}

\keyword{ htest }
\keyword{ robust }
