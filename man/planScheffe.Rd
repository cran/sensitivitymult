\name{planScheffe}
\alias{planScheffe}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matching}
\concept{Sheffe projection}
\concept{Multiple comparisons}
\title{
Combining One Planned Comparison and a Scheffe Correction For All Comparisons.
}
\description{
The function planScheffe() computes the critical values for a level
alpha test that combines one planned linear combination of
a K-dimensional multivariate Normal outcome and consideration of
all possible combinations correcting for multiple testing using
a Scheffe projection.
}
\usage{
planScheffe(K, alpha = 0.05)
}
\arguments{
  \item{K}{ An integer >=2 giving the number of outcomes to be compared.}
  \item{alpha}{ The level of the test, with 0 < alpha < 1.}
}
\details{
Although the calculation uses the multivariate Normal distribution,
a typical application uses K test statistics that are asymptotically
Normal.

The method is based on Rosenbaum (2017).  The example below
reproduces some of the comparisons in that manuscript.
}
\value{
\item{critical }{ critical is a vector with two elements, a and c.
The null hypothesis is rejected at level alpha if either the
Normal deviate for the planned comparison is >= a or if the square
of the Normal deviate for any comparison is >= b.  Then the
probability of a false rejection is <= alpha.}
\item{alpha }{ alpha is a vector with three elements, a, c and
joint.  The value of joint should equal the input value of alpha
aside from numerical errors of computation: it is the probability
of a false rejection using the joint test that rejects if
either of the two critical values in critical is exceeded.
In contrast, a is the probability that the planned deviate
will be >= critical[1] when the null hypothesis is true.
Also, c is the probability that at least one comparison
will have a squared deviate >= critical[2] when the
null hypothesis is true.}
}
\references{

Miller, R. G., Jr. (1981) Simultaneous Statistical Inference (2nd edition).
New York: Springer.  Section 2.2, pp. 48-67 discusses Scheffe projections.

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Rosenbaum, P. R. (2017) Combining planned and discovered comparisons
in observational studies.  Manuscript.

Scheffe, H. (1953) A method for judging all contrasts in the analysis
of variance.  Biometrika, 40, 87-104.

}
\author{
Paul R. Rosenbaum.
}
\note{
The method is based on Rosenbaum (2017).

The functions comparison() and principal() may be used to calculate
the standardized deviates that are compared to the critical values
from planScheffe.  Those functions have options for an a priori
comparison or consideration of all possible comparisons with a
Scheffe correction.  The function planScheffe provides a third
option: one planned comparison plus all possible comparisons.

}


\examples{
# Please READ the documentation for artcog, and in particular
# the distinction between simulated and actual data.
# The dontrun section refers to the acutal data and
# reproduces results in Rosenbaum (2017).

planScheffe(2,alpha=0.05)
# Interpretation of this output follows.
# Suppose there is a bivariate Normal outcome.  We specify
# one a priori linear combination of its two coordinates.
# We test that the expectation is (0,0) with known covariance
# matrix.  We compute the standardized difference for the
# a priori contrast, rejecting if it is >=1.895.  We also
# reject if we can find any linear combination of the two
# coordinates whose squared standardized difference is
# >=7.077.  The chance that we falsely reject a true
# null hypothesis is 0.05.  The chance of a false rejection
# using the a priori comparison is 0.029.  The chance of
# false rejection using any linear combination is 0.029.
#
# The a priori comparison could be the first principal
# component.  Using the principal() function with
# w=c(1,0) gives the deviate for the first principal
# component.  Exploring every w=c(w1,w2) gives the
# deviates that are squared for comparison with 7.077.

\dontrun{
# For this illustration, obtain the actual data,
# as described in the documentation for artcog.
# An illustration from Rosenbaum (2017) follows.
data(artcog)
attach(artcog)
# The comparison using the first principal component:
principal(cbind(words,wordsdelay,animals),arthritis,mset,
     w=c(1,0),gamma=1.396,detail=TRUE)
# The resulting deviate, 1.900 is slightly greated than 1.895,
# so the hypothesis of no effect would be rejected at 0.05 even if
# we allow for a bias in treatment assignment of gamma=1.396.
principal(cbind(words,wordsdelay,animals),arthritis,mset,
     w=c(1,-.075),gamma=1.396,detail=TRUE)
# The comparison w=c(1,-.075) yields a slightly larger
# deviate, 1.907, but 1.907^2 < 7.077, so this ad hoc
# comparison would not lead to rejection.
detach(artcog)
# Interpret gamma:
amplify(1.396,c(2,3))
amplify(1.4,c(2,3))
  }
}

\keyword{ htest }
\keyword{ robust }
