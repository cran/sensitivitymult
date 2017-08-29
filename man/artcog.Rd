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
The R package contains a simulated data set similar to actual
data from 2009-2011 Irish Longitudinal Study of Aging (TILDA) used
in Rosenbaum (2017).  Additionally, in the documentation
below, instructions are given for constructing the actual
data set after downloading a file from ICPSR.  The simulated
data may be used to try the methods in this package.  The
actual data may be used to replicate the calculations in
Rosenbaum (2017).  Please be careful to distinguish
the simulated data (with continuous outcomes) and the
actual data (with integer outcomes), as scientific
conclusions should not be based on the simulated data.

Instructions for creating the actual data are in
the example section, but are not executed because you
must obtain a file from ICPSR.

The simulated data were built from the actual data by
calculating the trivariate within group means and the
pooled within group covariance matrix.  Then a data
set of the same size was sampled from the trivariate
Normal distribution, using the actual means and
covariance matrix as the population parameters for
the simulation.  The simulation used the mvtnorm
package.  Although the data set consists of matched
triples, in the simulated version, the matched
sets are independent of the outcomes.

There are 219 matched triples containing
one individual with arthritis (arthritis=1) and two without
(arthritis=0).  There are three measures of cognitive
performance, words, wordsdelay and animals.
}
\usage{data("artcog")}
\format{
  A data frame with 657 observations on the following 5 variables.
  \describe{
    \item{\code{arthritis}}{1 if osteoarthritis, 0 if no arthritis}
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

}
\source{
Simulated data with a script for obtaining the actual data
from the Irish Longitudinal Study of Aging 2009-2011.
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
# data(artcog) returns the simulated example.
data(artcog)
# Three correlated outcomes.
cor(artcog[,2:4])
# See documentation for principal() for use of this example.

# The code below constructs the actual data, as distinct
# from the simulated example.  The lengthy list of numbers
# assembles the 219 matched triples, or 657 = 3*219 rows,
# from the larger TILDA data set.

\dontrun{
# Obtain from ICPSR the R data file
# ICPSR_34315-1IrishAging/34315-0001-Data.rda
# http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/34315?q=34315

# The data should have 8504 rows and 1992 columns

d<-da34315.0001
attach(d)

wordsC<-PH118
wordsI<-PH119
wordsC[wordsC<0]<-0
wordsI[wordsI<0]<-0
words<-wordsC+wordsI

wordsdelayC<-PH712
wordsdelayC[is.na(wordsdelayC)]<-0
wordsdelayC[wordsdelayC<=-1]<-0
wordsdelayI<-PH713
wordsdelayI[is.na(wordsdelayI)]<-0
wordsdelayI[wordsdelayI<=-1]<-0
wordsdelay<-wordsdelayC+wordsdelayI

animals<-PH125

arthritis<-PH301_03
osteoA<-PH304_1
z<-rep(NA,dim(d)[1])
z[arthritis==0]<-0
z[(arthritis==1)&(osteoA==1)]<-1

detach(d)

artcog<-data.frame(z,words,wordsdelay,animals)

who <- c(91, 4408, 7754, 129, 4716, 8383, 135, 6066,
 8028, 280, 894, 5300, 288, 151, 667, 298, 4889, 5977,
 333, 1100, 3707, 480, 696, 8148, 568, 372,
 7578, 584, 1852, 8057, 589, 3799, 6567, 590, 7422,
 8419, 609, 2825, 8272, 669, 1197, 8471, 684, 141,
 1847, 687, 2416, 7591, 771, 5239, 6986, 782,
 4857, 7654, 850, 885, 2239, 892, 2717, 7788, 929, 248,
 4740, 975, 1965, 8242, 1036, 6459, 7973, 1059, 1541,
 5901, 1103, 6518, 8264, 1160, 4798, 7330,
 1168, 4678, 7319, 1180, 152, 2735, 1191, 3740,
 7260, 1199, 26, 5209, 1252, 2615, 3251, 1444, 4790,
 7298, 1549, 898, 7630, 1587, 4418, 7122, 1596, 5875,
 8489, 1604, 3594, 7246, 1614, 3189, 7052, 1646,
 5415, 6828, 1708, 1634, 7029, 1760, 1950, 7815, 1840,
 5860, 8334, 1843, 6054, 7331, 1849, 5617, 8046, 1854,
 2890, 7703, 1885, 5846, 7247, 1896, 4365, 7803,
 1898, 3952, 4187, 1977, 544, 940, 1987, 768, 960,
 2029, 5363, 6293, 2161, 10, 4432, 2270, 5620, 7132,
 2330, 445, 1301, 2372, 1014, 1138, 2379, 3906,
 6183, 2386, 6226, 7203, 2417, 2458, 6616, 2437, 6262,
 7178, 2442, 3840, 8024, 2443, 4955, 5834, 2455, 1969,
 5967, 2457, 6962, 7560, 2466, 986, 2895, 2498, 2461,
 5876, 2522, 1837, 4803, 2618, 7279, 7764, 2734, 4005,
 4477, 2747, 221, 3837, 2763, 4440, 7863, 2765,
 6173, 7377, 2799, 7711, 7822, 2820, 2676, 7288, 2853,
 3035, 7518, 2914, 3142, 6891, 2952, 3081, 4908, 2969,
 3077, 6837, 3013, 747, 7614, 3107, 1754, 6564,
 3178, 2242, 4377, 3192, 260, 4530, 3246, 3019, 6478,
 3313, 4710, 7271, 3389, 356, 1796, 3481, 99, 491,
 3571, 658, 1410, 3693, 4341, 7624, 3694, 522,
 7702, 3704, 6532, 7171, 3705, 4973, 7131, 3806, 2163,
 5400, 3848, 4811, 7097, 3850, 2154, 5773, 3851, 3547,
 7613, 3862, 3357, 3370, 3877, 6186, 7990, 3913,
 455, 2883, 3931, 3548, 3699, 3933, 3210, 6164, 3935,
 4712, 7813, 3940, 5598, 7826, 3964, 2129, 8005, 3997,
 49, 1537, 4000, 3915, 5392, 4044, 3014, 6130,
 4052, 5208, 7213, 4186, 1586, 4249, 4264, 7058, 7182,
 4324, 3950, 7507, 4343, 3701, 6359, 4358, 567, 1020,
 4387, 2919, 4011, 4389, 5851, 7125, 4409, 3310,
 8100, 4427, 767, 2108, 4439, 1263, 6024, 4447, 3814,
 8373, 4478, 3493, 6743, 4479, 939, 2621, 4537, 1264,
 7942, 4608, 1797, 2987, 4633, 976, 1814, 4641,
 274, 1116, 4697, 4718, 7008, 4750, 2842, 5787, 4791,
 4386, 6966, 4812, 2817, 5640, 4815, 845, 5430, 4856,
 2288, 2289, 4887, 2182, 4874, 4942, 460, 4300,
 4945, 565, 3644, 4946, 487, 3369, 4953, 4352, 6709,
 4956, 2731, 3387, 4958, 4436, 6460, 4964, 3388, 6692,
 5078, 278, 963, 5110, 842, 4842, 5166, 600,
 1530, 5199, 1775, 6210, 5204, 4993, 8477, 5210, 2646,
 5563, 5291, 2957, 7777, 5325, 4881, 7053, 5342, 4385,
 6444, 5377, 3957, 4319, 5384, 3144, 7757, 5385,
 2813, 3054, 5386, 3636, 6185, 5474, 2507, 5085, 5488,
 4278, 5675, 5584, 2606, 5359, 5599, 3180, 7037, 5605,
 2459, 5304, 5637, 2581, 3621, 5641, 2781, 4302,
 5805, 6424, 7227, 5870, 492, 5847, 5909, 1750, 5158,
 5923, 3199, 6492, 6039, 4347, 4762, 6048, 5332, 7320,
 6080, 1992, 2830, 6091, 5213, 7045, 6099, 5167,
 6511, 6135, 5177, 6944, 6172, 2983, 6455, 6176, 2319,
 3737, 6189, 5525, 7257, 6196, 4423, 6893, 6256, 2639,
 5740, 6322, 1427, 2435, 6370, 7321, 7385, 6371,
 1212, 2423, 6417, 205, 1674, 6462, 2393, 2882, 6463,
 2170, 4765, 6496, 1630, 5048, 6519, 3058, 7498, 6901,
 5237, 7508, 6984, 3819, 6548, 7042, 2961, 6445,
 7057, 5457, 7984, 7061, 5401, 6049, 7093, 502, 3847,
 7094, 4717, 6348, 7096, 825, 7844, 7099, 2188, 8337,
 7251, 7423, 7576, 7269, 2616, 6401, 7270, 2394,
 5039, 7273, 2337, 4941, 7300, 2241, 4934, 7316, 2604,
 6369, 7355, 2113, 3880, 7402, 1456, 2378, 7473, 6368,
 7243, 7592, 5583, 7892, 7615, 855, 6924, 7684,
 6412, 6822, 7852, 405, 7077, 7862, 3623, 3990, 7879,
 2447, 6334, 7913, 3927, 5299, 7930, 5289, 5844, 7983,
 297, 7772, 8006, 3869, 6930, 8009, 2729, 6480,
 8081, 4700, 6560, 8109, 62, 8061, 8130, 3351, 4381,
 8149, 2854, 6513, 8157, 5220, 7559, 8184, 1523, 4195,
 8185, 1459, 3820, 8188, 117, 1050, 8206, 2513,
 3954, 8335, 2352, 4435, 8346, 5109, 8207)

artcog<-artcog[who,]
mset<-as.numeric(gl(219,3))
artcog<-cbind(artcog,mset)

rm(z,words,wordsdelay,animals,mset)
	}
}
\keyword{datasets}
