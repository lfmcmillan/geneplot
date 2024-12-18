% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extended_geneplot.R
\name{calc_geneplot_spdfs}
\alias{calc_geneplot_spdfs}
\title{Lower-level function for running the calculation step of \link{extended_geneplot}.
Carries out the GenePlot calculations and the saddlepoint distribution
calculations.}
\usage{
calc_geneplot_spdfs(
  dat,
  refpopnames,
  locnames,
  includepopnames = NULL,
  quantiles_vec,
  prior,
  logten = T,
  saddlepoint = T,
  leave_one_out = F,
  rel_tol = NULL,
  abs_tol = NULL,
  npts = 1000
)
}
\arguments{
\item{dat}{The data, in a data frame, with two columns labelled as 'id' and
'pop', and with two additional columns per locus. Missing data at any locus
should be marked as '0' for each allele. The locus columns must be labelled
in the format Loc1.a1, Loc1.a2, Loc2.a1, Loc2.a2, etc.
Missing data must be for BOTH alleles at any locus.
See \link{read_genepop_format} for details of how to import Genepop
format data files into the appropriate format.}

\item{refpopnames}{Character vector of reference population names, that must
match two values in the 'pop' column of \code{dat}. The SPDF methods
currently only work for a pair of baseline populations, so
\code{refpopnames} must be length 2.}

\item{locnames}{Character vector, names of the loci, which must match the
column names in the data so e.g. if dat has columns id, pop, EV1.a1,
EV1.a2, EV14.a1, EV14.a2, etc. then you could use 'locnames =
c("EV1","EV14") etc. The locnames do not need to be in any particular order
but all of them must be in \code{dat}.}

\item{includepopnames}{Character vector (default NULL) of population names to
be included in the calculations as comparison populations. The reference
populations are automatically used as comparison populations for each
other, but you can also add additional comparison populations using
\code{includepopnames}. For example, if the reference pops are Pop1 and
Pop2, and you have some new individuals which you have labelled as PopNew,
then use \code{includepopnames=c("PopNew")} to compare those individuals to
Pop1 and Pop2. You can specify the populations in any order, provided that
they are all in \code{dat}.}

\item{quantiles_vec}{Specify which quantiles to show on the plots, as a
vector of numbers between 0 and 1. They do not have to be ordered. If
NULL, quantiles will not be plotted.}

\item{prior}{(default "Rannala") String, either "Rannala" or "Baudouin",
giving the choice of prior parameter for the Dirichlet priors for the
allele frequency estimates. Both options define parameter values that
depend on the number of alleles at each locus, k.
"Baudouin" gives slightly more weight to rare alleles than "Rannala"
does, or less weight to the data, so Baudouin may be more suitable for
small reference samples, but there is no major difference between them.
For more details, see McMillan and Fewster (2017), Biometrics.
Additional options are "Half" or "Quarter" which specify parameters 1/2
or 1/4, respectively. These options have priors whose parameters do not
depend on the number of alleles at each locus, and so may be more suitable
for microsatellite data with varying numbers of alleles at each locus.}

\item{logten}{(default TRUE) Boolean, indicates whether to use base 10 for the
logarithms, or base e (i.e. natural logarithms). logten=TRUE is default
because it's easier to recalculate the original non-log numbers in your
head when looking at the plots. Use FALSE for natural logarithms.}

\item{saddlepoint}{(default TRUE) If TRUE, use saddlepoint approximation to
impute Log-Genotype Probability for individual genotypes with missing
data. If not, use an empirical approximation to impute the LGPs.
Defaults to TRUE because the side plots in the extended GenePlot use the
saddlepoint approximation process.}

\item{leave_one_out}{(default TRUE) Boolean, indicates whether or not to
calculate leave-one-out results for any individual from the reference
pops. If TRUE, any individual from a reference population will have
their Log-Genotype-Probability with respect to their own reference
population after temporarily removing the individual's genotype from the
sample data for that reference population. The individual's
Log-Genotype-Probabilities with respect to all populations they are not
a member of will be calculated as normal.
We STRONGLY RECOMMEND using leave-one-out=TRUE for any small reference
samples (<30).}

\item{rel_tol}{(default NULL) Specify the relative tolerance for the
numerical integration function that is used to calculate the overlap area
and also the normalization constants for the various distributions. The
default value corresponds to the \code{integrate} default i.e.
\code{.Machine$double.eps^0.25}.}

\item{abs_tol}{(default NULL) Specify the absolute tolerance for the
numerical integration function that is used to calculate the overlap area
and also the normalization constants for the various distributions. The
default value corresponds to the \code{integrate} default i.e.
\code{.Machine$double.eps^0.25}.}

\item{npts}{(default 1000) Number of values to use when calculating numerical
integrals (for the overlap area measures, and for the normalization
constants of the distributions). Increasing this value will increase the
precision of the numerical integrals but will also increase the
computational cost. Reducing this below 1000 may save some computation
time if you are not too concerned with the precision of the results.}
}
\value{
A list with the following components:

   \code{logprob} GenePlot calculation results: Log-Genotype Probability
   values for all individuals with respect to all of the reference populations.
   If there are individuals with missing values, their raw LGPs are shown
   which are based on the loci that are present, and also the imputed LGPs
   for the full set of loci. This output is the same as the output from
   \link{calc_logprob}.

   \code{spdf_vals1} Saddlepoint distribution approximations with the first
   reference population as the baseline. List contains \code{xvals},
   \code{yvals}, \code{wvals} and \code{zvals}, the raw distribution curves
   for plotting the distributions.
   This list also contains \code{oavals}, \code{probvals} and
   \code{diffvals}, which are the Overlap Area, Incumbent Selection
   Probability and Home Assignment Probability values for the given baseline
   population with all the other populations as comparisons. The sub-list
   also records the quantile values requested, the name of the given baseline
   pop for this sub-list as \code{refpopA}, the indices of the baseline pop
   and comparison pop in the reference pops, the name of the other reference
   population and the number of other reference populations (always equal to
   one).

   \code{spdf_vals2} As for \code{spdf_vals1}, but with the second reference
   population as the baseline.
}
\description{
Use before \link{prepare_plot_params} and \link{plot_geneplot_spdfs}.
}
