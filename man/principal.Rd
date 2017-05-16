\name{principal}
\alias{principal}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Sheffe projection}
\concept{Principal Components}
\title{
Sensitivity Analysis for Principal Components of M-Scores for Several Outcomes in an Observational Study.
}
\description{
For k outcomes in a matched observational study, principal() computes
the M-scores for the outcomes one at a time, computes the principal
components of the M-scores, and uses some of the larger principal
components as the outcomes in a sensitivity analysis.  The user
controls the number of components using w: (i) if is.null(w)==TRUE,
then the first principal component is used in a one-sided test,
(ii) if length(w)==1, then w=1 and w=-1 both use the first
principal component, but direct attention to the upper or
lower tails, respectively, (iii) if length(w)>1, then
the first length(w) principal components are used with
weights w; e.g., w=c(1,1) adds the first two principal
components together.  Setting Scheffe=TRUE with length(w)=2
permits the user to test every linear combination of the
first two principal components -- that is, every with
length(w)=2 -- while controlling the family-wise error
rate.  Every matched set contains one treated subject and
one or more controls.
}
\usage{
principal(y,z,mset,w=NULL,gamma=1,inner=0,trim=3,lambda=0.5,
                     TonT=FALSE,apriori=FALSE,Scheffe=FALSE,detail=FALSE,
                     cor=FALSE)
}
\arguments{
  \item{y}{ A matrix of responses with no missing data.  Different
  columns of y are different variables, and there are k=dim(y)[2] variables.
  If present, the column names of y are used to label output.}
  \item{z}{ Treatment indicators, z=1 for treated, z=0 for control with
  length(z)==dim(y)[1].
}
  \item{mset}{ Matched set indicators, 1, 2, ..., sum(z) with
  length(mset)==dim(y)[1].  The vector mset may contain integers
  or may be a factor.
}
  \item{w}{ A vector of weights to be applied to principal
  components.  There are k=dim(y)[2] variables in y.
  (i) If is.null(w)=TRUE or if length(w)=1 with w!=0, then
  the test is applied to the first principal component of the nvars
  M-test scores, and no adjustment for multiple testing is needed.
  (ii) If length(w)>1, then w determines a comparison among the first
  length(w) principal components of the k M-scores.  At least
  one weight must be nonzero.  The meaning of the weights is affected
  by whether cor=FALSE (the default) or cor=TRUE.  If Scheffe=TRUE,
  the dimensionality of the Scheffe correction is length(w),
  so w=c(1,1) is different from w=(1,1,0), because the latter
  implies the investigator may consider the third principal
  component in some comparison, even though w=(1,1,0) ignores it.
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
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
See inner.  Unlike other functions in this package, principal() requires
trim<Inf.  When trim<Inf, the M-statistics for different outcomes have
been scaled so their magnitudes are comparable, but this would not be
true for trim=Inf.  If you would like to do analogous calculations without
trimming, then give the comparison() function principal component scores rather than
data for y, set trim=Inf and inner=0, and the w in that function will be
applied to the principal component scores that you supplied.
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
  or P-value bound is 1-pchisq(max(0,deviate)^2,length(w)).  If Scheffe=FALSE
  and apriori=FALSE, then the deviate is returned, but no P-value is given.
  See Rosenbaum (2016).  Note carefully that length(w) determines the extent of
  the correction for multiplicity, so that, for example, w=c(1,0) focuses
  attention on the first principal component but allows for consideration of all
  comparisons using the first two principal components.  See the discussion of
  w above and the examples below.  A Scheffe correction entitles you to look in
  both tails, which you do by considering both w and -w.
  See the planScheffe() function for a combination
  of an apriori and Scheffe comparisons.
}
\item{detail}{
  If detail=TRUE, then some detail from the princomp() function in the
  stats package is returned.
}
\item{cor}{
  If cor=FALSE, the principal components of the M-scores are computed
  from the covariance matrix, but if cor=TRUE, then they are computed
  from the correlation matrix.  Because the columns of y were scaled
  using lambda and scored by the same \eqn{\psi}-function, they are
  on a common scale, and it is reasonable to compute the principal
  components from the covariance matrix with cor=FALSE, the default.
  Setting cor=TRUE standardizes the columns of y twice, so perhaps it
  is pointless.  In any event, because the weights w are applied to
  the principal components, and the latter are affected by cor, it
  follows that the meaning of w is affected by the value of cor.
}
}
\details{
If y has k columns for k outcomes, then comparison computes k M-scores,
one for each outcome, computes principal components from these scores,
combines the scores into a single comparison using w, and
computes a one-sided, upper-tailed deviate for a randomization test
or a sensitivity analysis, as described in Rosenbaum (2007, 2016).

Outcomes are scaled using by the \eqn{\lambda} quantile of
the absolute differences before applying the \eqn{\psi}-function.  In this
sense, when trim<Inf, the M-scores share a common scaling before
principal components are computed.

Taking Scheffe=TRUE and w=(w1,w2) for all w1 and w2 considers all
comparisons based on the first two principal components.

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
\item{weights }{The weights are returned.}
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

Rosenbaum, P. R. (2017) Combining planned and discovered comparisons
in observational studies.  Manuscript.  (The artcog example below
is from this manuscript.)

Rosenbaum, P. R., & Rubin, D. B. (1985). The bias due to incomplete
matching. Biometrics, 41, 103-116.

}
\author{
Paul R. Rosenbaum.
}
\note{
For confidence intervals for individual outcomes, use function senmCI().

Under Fisher's hypothesis of no treatment effect, the principal
components of the outcomes are unaffected by the treatment, so
they may be used in randomization tests of no effect.  However, this logic
does not permit confidence intervals for the magnitude of effect
on a principal component.

The principal() function computes principal components of M-scores, not
of the outcomes themselves.  This has various implications.  The
M-scores share a common, resistant scaling, so it is reasonable to
consider principal components of the covariance matrix of the M-scores.
In contrast, M-tests computed from principal components of outcomes
are not resistant to outliers because the components themselves are
not resistant to outliers.  M-scores add to zero within each matched
set; see the example for the mscorev() function.  In this specific
and limited sense, variation among M-scores reflects variation within
matched sets rather than variation between matched sets.  For example,
if the matched sets had been exactly matched for age, then the M-scores
would be uncorrelated with age.  In contrast,
principal components of outcomes reflect both variation within and
variation between matched sets.  For instance, principal components of
outcomes might be correlated with age even if the sets had been matched
exactly for age.  When the matched set size is
variable, the M-scores incorporate variable weights, and the
principal components are affected by these weights.  For this reason,
principal components of M-scores are more interpretable when every
matched set has the same size, say matched pairs or matching 1-to-2,
and they may be difficult to interpret if the set sizes vary
widely, say 1-1 mixed with 1-5.  In thinking about the relationship
between outcomes and their M-scores, it can be helpful to examine
the small, univariate example for the mscorev() function.
}


\examples{
data(artcog)
attach(artcog)
# A randomization test using the first principal component for the three memory scores.
# The loadings show that the first component gives positive weight to each memory score.
principal(cbind(words,wordsdelay,animals),arthritis,mset,w=1,apriori=TRUE,detail=TRUE)
#
# The comparison above is insensitive to a bias of gamma=1.45
principal(cbind(words,wordsdelay,animals),arthritis,mset,w=1,gamma=1.45,apriori=TRUE,detail=TRUE)
#
# gamma=1.45 is an unobserved covariate that more than triples the odd of a poor memory score
# and more than doubles the odds of arthritis.
amplify(1.45,c(2,3,4))
#
# Although the first principal component is insensitive to a bias of gamma=1.45, each
# of the three individual variables is sensitive to a bias of gamma=1.45
senm(words,arthritis,mset,gamma=1.45)
senm(wordsdelay,arthritis,mset,gamma=1.45)
senm(animals,arthritis,mset,gamma=1.45)
#
# Although not particularly useful or enlightening in this one example, we can
# explore all weighted combinations of the first two principle components,
# correcting for multiple testing using Scheffe projections for dimension 2.
# This would be more interesting in an example with 50 outcomes, where we
# might want to reduce the dimensionality to 2 or 3 from 50, rather than to 1.
# We will do calculations for gamma=1.25.  A gamma=1.25 is an unobserved
# covariate that doubles the odds of arthritis and doubles
# the odds of a worse memory score.
amplify(1.25,2)
# The deviate is the same but the corrected P-value is different if w=1 or w=c(1,0),
# because the former is doing a single one-sided test, while the latter is anticipating
# consideration of all possible combinations of the first two components.
principal(cbind(words,wordsdelay,animals),arthritis,mset,w=1,gamma=1.25,apriori=TRUE,detail=TRUE)
principal(cbind(words,wordsdelay,animals),arthritis,mset,w=c(1,0),gamma=1.25,Scheffe =TRUE)
# A weighted combination of the first two principal components, w=c(1,-.1), is ever so
# slightly less sensitive than using the first component alone.
principal(cbind(words,wordsdelay,animals),arthritis,mset,w=c(1,-.1),gamma=1.25,Scheffe =TRUE)
detach(artcog)
}

\keyword{ htest }
\keyword{ robust }
