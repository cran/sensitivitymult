\name{artcog}
\alias{artcog}
\docType{data}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matched sets}
\concept{Matching with multiple controls}
\concept{Alzheimers disease}
\title{
Arthritis and cognition in the elderly.
}
\description{
Data from the 2009-2011 Irish Longitudinal Study of Aging relating
arthritis and cognitive performance; see the Details for
motivation.  There are 219 matched triples containing
one individual with arthritis (arthritis=1) and two without
(arthritis=0).  There are three measures of cognitive
performance, words, wordsdelay and animals.
}
\usage{data("artcog")}
\format{
  A data frame with 657 observations on the following 7 variables.
  \describe{
    \item{\code{arthritis}}{1 if osteoarthritis, 0 if no arthritis}
    \item{\code{arthritisAge}}{Decade in which arthritis was diagnosed.
    The value is -99 for people without arthritis.}
    \item{\code{words}}{Individuals are read a list of words and are
    immediately asked to recall as many as they can.}
    \item{\code{wordsdelay}}{After a delay and another task, individuals
    are again asked to recall as many of the words as they can.}
    \item{\code{animals}}{Individuals are asked to name all of the
    animals they can think of.}
    \item{\code{mset}}{Indicator of the matched set: 1, 2, ..., 219.}
  }
}
\details{
A theory that NSAIDs reduce the risk of Alzheimer's disease has often been
examined by comparing elderly people with and without arthritis, reasoning
that many people with arthritis have consumed NSAIDs in quantity for a
long period; see McGeer et al. (1996).  This comparison does not
ask a person with Alzheimer's disease to recall past use of NSAIDs.

Triples were matched for age, sex, education and mother's education.

Everyone is 75 years old or older.

Presumably, many individuals whose arthritis was diagnosed earlier
have taken NSAIDs for a longer time.
}
\source{
Irish Longitudinal Study of Aging 2009-2011.
}
\references{
Irish Longitudinal Study of Aging.  http://tilda.tcd.ie/

McGeer, P. L., Schulzer, M., and McGeer, E. G. (1996). Arthritis and anti-
inflammatory agents as possible protective factors for alzheimer's disease.
Neurology, 47, 425-432.

Rosenbaum, P. R. (2017) Combining planned and discovered comparisons
in observational studies.  Manuscript. (The artcog example is discussed
in this manuscript.)
}
\examples{
data(artcog)
# Three correlated outcomes.
cor(artcog[,4:6])
# See documentation for principal() for use of this example.
}
\keyword{datasets}
