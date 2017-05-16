\name{teeth}
\alias{teeth}
\docType{data}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matched pairs}
\title{
Smoking and Periodontal Disease.
}
\description{
Data from NHANES 2011-2012 concerning periodonal disease in 441 matched
pairs of smokers and nonsmokers.
}
\usage{data("teeth")}
\format{
  A data frame with 882 observations on the following 4 variables.
  \describe{
    \item{\code{mset}}{Matched pair indicator: 1, 2, ..., 441.}
    \item{\code{smoker}}{Treatment indicator: 1 if current smoker, 0 if never smoker}
    \item{\code{either4up}}{Measure of periodontal disease for upper teeth; see
    Details.}
    \item{\code{either4low}}{Measure of periodontal disease for lower teeth; see
    Details.}
    \item{\code{cigsperday}}{Cigarettes smoked per day.  Zero for nonsmokers.}
  }
}
\details{
Smoking is believed to cause periodontal disease; see Tomar and Asma (2000).
Using more recent data from NHANES 2011-2012, the data describe 441 matched
pairs of a daily smoker and a never smoker.  Daily smokers smoked every day
of the last 30 days.  Never smokers smoked fewer than 100 cigarettes in their
lives, do not smoke now, and had no tobacco use in the previous five days.

Pairs are matched for education, income, age, gender and black race.

Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding
4 wisdom teeth.  Pocket depth and loss of attachment are two
complementary measures of the degree to which the gums have separated
from the teeth; see Wei, Barker and Eke (2013).  Pocket depth and loss of
attachment are measured at six locations on each tooth, providing the
tooth is present.  A measurement at a location was taken to exhibit
disease if it had either a loss of attachement >=4mm or a pocked depth
>=4mm, so each tooth contributes a score from 0 to 6.  Upper and lower
are the number of measurements exhibiting disease for upper and
lower teeth.

This example is from Rosenbaum (2016) where more information may be
found.
}
\source{
National Health and Nutrition Examination Survey (NHANES), 2011-2012.
https://www.cdc.gov/nchs/nhanes/
}
\references{
Rosenbaum, P. R. (2016). Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Tomar, S. L. and Asma, S. (2000). Smoking attributable periodontitis in the United States: Findings from NHANES III. J. Periodont. 71, 743-751.

Wei, L., Barker, L. and Eke, P. (2013). Array applications in determining periodontal disease
measurement. SouthEast SAS User's Group. (SESUG2013) Paper CC-15, analytics.ncsu.edu/
sesug/2013/CC-15.pdf.
}
\examples{
data(teeth)
attach(teeth)
# The following calculation reproduces the comparison
# in expression (5.1) of Rosenbaum (2016, p. 1464)
comparison(cbind(either4low,either4up),smoker,
   mset,c(.714,.286),gamma=2.2,trim=2.5,Scheffe=TRUE)
# Note that Rosenbaum (2016) used trim=2.5, but comparison()
# has a default of trim=3.
# The parameter gamma=2.2 is given alternative interpretations
# in Rosenbaum (2016, p. 1465) as follows:
amplify(2.2,c(3,4,5,6,7))
# The calculation (Rosenbaum 2016, p. 1465) for lower teeth alone is:
comparison(cbind(either4low,either4up),
  smoker,mset,c(1,0),gamma=2.2,trim=2.5,apriori=TRUE)
# Because w = c(1,0) ignores upper teeth, it may also be done as follows.
# Remove the comment sign to execute senmCI which is a little slow.
# senm(either4low,smoker,mset,gamma=2.2,trim=2.5)
# The calculations that follow reproduce the intervals from
# section 5.1 of Rosenbaum (2016, p. 1466)
# Remove the comment sign to execute senmCI which is a little slow.
# senmCI(teeth$either4low,teeth$smoker,teeth$mset,trim=2.5,gamma=1.5)
# The example that follows uses inner=0.5 as in Rosenbaum (2016,
# p. 1466, section 5.2):
comparison(cbind(either4low,either4up),smoker,
  mset,c(.714,.286),gamma=2.2,inner=.5,trim=2.5,Scheffe=TRUE)
comparison(cbind(either4low,either4up),smoker,
  mset,c(.714,.286),gamma=2.37,inner=.5,trim=2.5,Scheffe=TRUE)
detach(teeth)
}
\keyword{datasets}
