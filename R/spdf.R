#' Run GenePlot directional analysis of the connectivity of two populations.
#'
#' This function calculates three measures of genetic differentiation between
#' populations that relate to assignment probabilities. It also calculates the
#' values for the auxiliary curves in the extended GenePlot (see also
#' \code{extended_geneplot}.
#'
#' All three measures are \strong{directional},
#' meaning that for two populations, A and B, there is an A to B value and a
#' separate B to A value. The function calculates all three measures in both
#' directions, i.e. with each of the reference populations as the baseline in
#' turn.
#'
#' 1) Overlap Area (OA): the area of overlap between the two auxiliary curves in
#' one side of the extended GenePlot.
#'
#' 2) Incumbent Selection Probability (ISP): what is the probability that if you
#' take a random individual from baseline population B and another from
#' comparison population A, the individual from B has a better fit to B than the
#' individual from A has to population B?
#' As an analogy, if you asked a random Dutch child and a random English child
#' to take an English language test, what is the probability that the English
#' child gets a higher mark than the Dutch child? (England is the baseline
#' population B in this analogy.)
#'
#' 3) Home Assignment Probability (HAP): what is the probability that if you
#' take a random individual from baseline population B, the individul has
#' better fit to population B than population A?
#' As an analogy, if you pick a random English child and give them an English
#' test and a Dutch test, what is the probability that they will do better at
#' the English test? (England is again the baseline population in this analogy.)
#'
#' The three measures are based on saddlepoint approximations to various
#' distributions relating to the fit of potential individual genotypes to two
#' candidate source populations. "spdf" stands for Saddlepoint Probability
#' Density Function, because the values calculated are based on probability
#' density functions.
#'
#' By default the three directional measures are shown as matrices with the
#' baseline populations as the \strong{columns} and the comparison populations
#' as the \strong{rows}.
#'
#' By default, the outputs also include the full genetic distributions of the
#' two reference populations, approximated using the saddlepoint approximation.
#' The genetic distribution of a population B is the distribution of
#' log-genotype probabilities for all possible genotypes that could arise from
#' the population, given its estimated allele frequencies. It shows the range of
#' fits to the population that are possible for individuals that could arise
#' from it.
#'
#' @param dat The data, in a data frame, with two columns labelled as 'id' and
#'    'pop', and with two additional columns per locus. Missing data at any
#'    locus should be marked as '0' for each allele.
#'    The locus columns must be labelled in the format Loc1.a1, Loc1.a2,
#'    Loc2.a1, Loc2.a2, etc.
#'    Missing data must be for BOTH alleles at any locus.
#'    See \code{\link{read_genepop_format}} for details of how to import Genepop
#'    format data files into the appropriate format.
#'    The user must supply either an input data frame \code{dat} and a vector of
#'    loci names \code{locnames}, or a list of allele frequencies at all relevant
#'    loci for all the populations, \code{allele_freqs}. The \code{allele_freqs}
#'    list can be obtained as an attribute of the output of GenePlot \code{logprob}
#'    calculations, and may cover more populations than the two used in this SPDF
#'    analysis (the SPDF analysis populations are indicated in the argument
#'    \code{refpopnames})
#'
#' @param refpopnames Character vector of reference population names, that must
#'    match two values in the 'pop' column of \code{dat}. The SPDF methods
#'    currently only work for a pair of baseline populations, so
#'    \code{refpopnames} must be length 2.
#'
#' @param locnames Character vector, names of the loci, which must match the
#'      column names in the data so e.g. if dat has columns
#'      id, pop, EV1.a1, EV1.a2, EV14.a1, EV14.a2, etc.
#'      then you could use 'locnames = c("EV1","EV14") etc.
#'      The locnames do not need to be in any particular order but all of them
#'      must be in \code{dat}.
#'
#' @param includepopnames Character vector (default NULL) of population names to
#'      be included in the calculations as comparison populations. The reference
#'      populations are automatically used as comparison populations for each
#'      other, but you can also add additional comparison populations using
#'      \code{includepopnames}.
#'      For example, if the reference pops are Pop1 and Pop2, and you have some
#'      new individuals which you have labelled as PopNew, then use
#'      \code{includepopnames=c("PopNew")} to compare those individuals to Pop1
#'      and Pop2.
#'      You can specify the populations in any order, provided that they are all
#'      in \code{dat}.
#'
#' @param allele_freqs (default=NULL) Alternative input format, which you can
#'      supply instead of \code{dat} and \code{locnames}. You can calculate the
#'      allele frequencies object using the \code{calc_logprob} function.
#'
#' @param display_names (default=refpopnames) Use this to supply alternative
#'      display names for the populations. The refpopnames, as columns in the
#'      dataset, cannot have spaces, for example, whereas the display names can
#'      have spaces.
#'
#' @param plot_spdfs (default=FALSE) If true, display plots of the genetic
#'      distributions as part of the output. Not needed if calculating the
#'      values for an extended GenePlot.
#'      The first two plots have the first reference pop as the baseline, then
#'      the last two have the second reference pop as the baseline. Within each
#'      pair, one plot shows the baseline population genetic distribution and
#'      the distribution of the comparison population relative to it, and the
#'      second shows the distribution of the differences between fit to the
#'      baseline and the comparison, for all genotypes that could arise from
#'      the baseline population. In each "differences" plot, the values above 0
#'      are from genotypes that have a better fit to their own population, the
#'      baseline, than to the comparison population, and the values below 0 are
#'      from genotypes that have a better fit to the comparison population than
#'      their own baseline population.
#'
#' @param difference_threshold (default=0) When this is zero, the Home Assignment
#'      Probability is the probability of a random individual from the baseline
#'      B having a better fit to B than A. If you want to instead calculate the
#'      probability that the individual from B has 10x better fit to B than A,
#'      then set \code{difference_threshold} to 1, because log10(10) = 1
#'      (for \code{logten = TRUE}) and the probabilities are on a log scale.
#'      Positive values make the measure calculation more conservative,
#'      negative values make the measure calculation less conservative.
#'      Only applies to the Home Assignment Probability, not the other measures.
#'
#' @param show_quantiles (default=TRUE) If TRUE, show quantiles on the
#'      distribution plots as vertical lines. For example, on the distribution
#'      of the baseline population, if one of the quantiles is 0.01 (1%) then
#'      the vertical line will show that 1% of the genotypes that could arise
#'      from the population will have a worse fit than the quantile value and
#'      99% of them will have a better fit than the quantile value.
#'      Ignored if \code{plot_spdfs = FALSE}.
#'
#' @param quantiles_vec (default=0.01) Specify which quantiles to show on the
#'      distribution plots, as a vector of numbers between 0 and 1. They do not
#'      have to be ordered.
#'
#' @param prior (default="Rannala") String, either "Rannala" or "Baudouin",
#'      giving the choice of prior parameter for the Dirichlet priors for the
#'      allele frequency estimates. Both options define parameter values that
#'      depend on the number of alleles at each locus, k.
#'      "Baudouin" gives slightly more weight to rare alleles than "Rannala"
#'      does, or less weight to the data, so Baudouin may be more suitable for
#'      small reference samples, but there is no major difference between them.
#'      For more details, see McMillan and Fewster (2017), Biometrics.
#'      Additional options are "Half" or "Quarter" which specify parameters 1/2
#'      or 1/4, respectively. These options have priors whose parameters do not
#'      depend on the number of alleles at each locus, and so may be more suitable
#'      for microsatellite data with varying numbers of alleles at each locus.
#'
#' @param leave_one_out Boolean (default TRUE), indicates whether or not to
#'      calculate leave-one-out results for any individual from the reference
#'      pops. If TRUE, any individual from a reference population will have
#'      their Log-Genotype-Probability with respect to their own reference
#'      population after temporarily removing the individual's genotype from the
#'      sample data for that reference population. The individual's
#'      Log-Genotype-Probabilities with respect to all populations they are not
#'      a member of will be calculated as normal.
#'      We STRONGLY RECOMMEND using leave-one-out=TRUE for any small reference
#'      samples (<30).
#'
#' @param logten (default TRUE) Boolean, indicates whether to use base 10 for the
#'      logarithms, or base e (i.e. natural logarithms). logten=TRUE is default
#'      because it's easier to recalculate the original non-log numbers in your
#'      head when looking at the plots. Use FALSE for natural logarithms.
#'
#' @param calc_positive_stats (default=TRUE) Which form of the genetic measures
#'      to calculate. If TRUE, the Incumbent Selection Probability calculates
#'      the probability that for two random individuals from baseline B and
#'      comparison A, the one from B has the best fit to B. The Home Assignment
#'      Probability calculates the probability that a random individual from
#'      baseline B has a better fit to B than A.
#'      If FALSE, the Incumbent Selection Probability becomes instead the
#'      Interloper Selection Probability, ie. the probability that the individual
#'      from A has a better fit to B than the individual from B. And the Home
#'      Assignment Probability becomes the Away Assignment Probability, i.e. the
#'      probability that the individual from B has a better fit to A than to B.
#'
#' @param calc_details (default=FALSE) If TRUE, the function displays additional
#'      statistics relating to the genetic distribution of each of the reference
#'      populations.
#'
#' @param calc_vecs (default=TRUE) If TRUE, the function calculates the three
#'      measures in both directions, and calculates the full distribution curves
#'      required for the extended GenePlot. Therefore leave this as the default
#'      (TRUE) if using this function prior to plotting the extended GenePlot.
#'
#' @param output_as_vectors (default=FALSE) By default, the output measures are
#'      for the two reference populations and any additional included populations
#'      are collated into a matrix of values, where each column represents one
#'      of the reference populations as the baseline, and the comparison
#'      populations are in the rows. If \code{output_as_vectors} is TRUE, the
#'      outputs are instead in vector format, where the columns indicate the
#'      baseline and comparison populations for each value (e.g. PopB.PopA is
#'      the value for baseline B and comparison A).
#'
#' @param rel_tol (default=NULL) Specify the relative tolerance for the
#'      numerical integration function that is used to calculate the overlap area
#'      and also the normalization constants for the various distributions. The
#'      default value corresponds to the \code{integrate} default i.e.
#'      \code{.Machine$double.eps^0.25}.
#'
#' @param abs_tol (default=NULL) Specify the absolute tolerance for the
#'      numerical integration function that is used to calculate the overlap area
#'      and also the normalization constants for the various distributions. The
#'      default value corresponds to the \code{integrate} default i.e.
#'      \code{.Machine$double.eps^0.25}.
#'
#' @param npts (default=1000) Number of values to use when calculating numerical
#'      integrals (for the overlap area measures, and for the normalization
#'      constants of the distributions). Increasing this value will increase the
#'      precision of the numerical integrals but will also increase the
#'      computational cost. Reducing this below 1000 may save some computation
#'      time if you are not too concerned with the precision of the results.
#'
#' @param only_plot_baseline_pop (default=FALSE) Ignored if \code{plot_spdfs=FALSE}.
#'      By default, the first and third plots in the output show the genetic
#'      distribution of the baseline pop, and also the genetic distribution of
#'      the comparison pop against the baseline pop. If TRUE, only plot the
#'      distribution of the baseline pop.
#'
#' @param show_statistics_on_plot (default=TRUE) Ignored if \code{plot_spdfs=FALSE}.
#'      If TRUE, state the incumbent selection probabilities and home assignment
#'      probabilities beneath their respective plots.
#'
#' @param line_cols (default=NULL) Ignored if \code{plot_spdfs=FALSE}. Vector of
#'      line colours for the populations, starting with the two reference
#'      populations in the same order as in \code{refpopnames}. You can use any
#'      R colour specification, including named colours.
#'
#' @param line_widths (default=NULL) Ignored if \code{plot_spdfs=FALSE}. Vector
#'      of line widths for the populations, starting with the two reference
#'      populations in the same order as in \code{refpopnames}. Standard R line
#'      widths.
#'
#' @param title_text (default=NULL) Ignored if \code{plot_spdfs=FALSE}. Title
#'      text for the Incumbent Selection Probability plots.
#' @param title_text_difference (default=NULL) Ignored if \code{plot_spdfs=FALSE}.
#'      Title text for the Home Assignment Probability plots.
#'
#' @param ISP_xlim (default=NULL) Ignored if \code{plot_spdfs=FALSE}. x-axis
#'      limits for the Incumbent Selection Probability plots.
#' @param ISP_ylim (default=NULL) Ignored if \code{plot_spdfs=FALSE}. y-axis
#'      limits for the Incumbent Selection Probability plots.
#'
#' @param HAP_xlim (default=NULL) Ignored if \code{plot_spdfs=FALSE}. x-axis
#'      limits for the Home Assignment Probability plots.
#' @param HAP_ylim (default=NULL) Ignored if \code{plot_spdfs=FALSE}. y-axis
#'      limits for the Home Assignment Probability plots.
#'
#' @param ISP_legend_xy (default=NULL) Ignored if \code{plot_spdfs=FALSE}.
#'      x-y position for the legend in Incumbent Selection Probability plots.
#' @param HAP_legend_xy (default=NULL) Ignored if \code{plot_spdfs=FALSE}.
#'      x-y position for the legend in Home Assignment Probability plots.
#'
#' @param axis_labels (default=TRUE) Ignored if \code{plot_spdfs=FALSE}. If TRUE,
#'      include axis labels on all the SPDF plots.
#'
#'
#' @returns A list with components:
#'
#'     \code{overlap_results}: Overlap Area values calculated using each of the
#'     two reference populations as the baseline in turn. Comparison populations
#'     are the other reference population and any additional included
#'     populations.
#'
#'     \code{incumbent_results}: Incumbent Selection Probability values
#'     calculated using each of the two reference populations as the baseline in
#'     turn. Comparison populations are the other reference population and any
#'     additional included populations.
#'
#'     \code{home_assignment_results}: Home Assignment Probability values
#'     calculated using each of the two reference populations as the baseline in
#'     turn. Comparison populations are the other reference population and any
#'     additional included populations.
#'
#'     \code{spdf_results}: List with two sub-lists, one for each reference
#'     population as the baseline. Within a sublist, the entries \code{xvals},
#'     \code{yvals}, \code{wvals} and \code{zvals} are blank unless argument
#'     \code{plot_sdfs = TRUE} or \code{calc_vecs = TRUE}. They are the raw
#'     distribution curves for plotting the ISP and HAP values.
#'     Each sub-list also contains \code{oavals}, \code{probvals} and
#'     \code{diffvals}, which are the Overlap Area, Incumbent Selection
#'     Probability and Home Assignment Probability values for the given baseline
#'     population with all the other populations as comparisons.
#'     The sub-list also records the quantile values requested, the name of the
#'     given baseline pop for this sub-list as \code{refpopA}, the indices of
#'     the baseline pop and comparison pop in the reference pops, the name of
#'     the other reference population and the number of other reference populations
#'     (always equal to one).
#'
#' @references McMillan, L. F., "Concepts of statistical analysis,
#'   visualization, and communication in population genetics" (2019). Doctoral
#'   thesis, https://researchspace.auckland.ac.nz/handle/2292/47358.
#'
#' @references McMillan, L. and Fewster, R. "Visualizations for genetic assignment analyses
#'  using the saddlepoint approximation method" (2017) \emph{Biometrics}.
#'  Rannala, B., and Mountain, J. L. (1997). Detecting immigration by using
#'  multilocus genotypes. \emph{Proceedings of the National Academy of Sciences}
#'  \strong{94}, 9197--9201.
#'
#' @references Piry, S., Alapetite, A., Cornuet, J.-M., Paetkau, D., Baudouin, L., and
#'  Estoup, A. (2004). GENECLASS2: A software for genetic assignment and
#'  first-generation migrant detection. \emph{Journal of Heredity} \strong{95},
#'  536--539.
#'
#' @author Saddlepoint approximation to distributions of genetic fit to
#'   populations developed by McMillan and Fewster,  based on calculations of
#'   Log-Genotype-Probability from the method of Rannala and Mountain (1997) as
#'   implemented in GeneClass2, updated to allow for individuals with missing
#'   data and to enable accurate calculations of quantiles of the
#'   Log-Genotype-Probability distributions of the reference populations. See
#'   McMillan and Fewster (2017) for details.
#'
#' @export
calc_spdfs <- function(dat=NULL,
                       refpopnames,
                       locnames=NULL,
                       includepopnames=NULL,
                       allele_freqs=NULL,
                       display_names=refpopnames,
                       plot_spdfs=FALSE,
                       difference_threshold=0,
                       show_quantiles=TRUE,
                       quantiles_vec=c(0.01,1.00),
                       prior='Rannala',
                       logten=T,
                       leave_one_out=T,
                       calc_positive_stats=TRUE,
                       calc_details=FALSE,
                       calc_vecs=TRUE,
                       output_as_vectors=FALSE,
                       rel_tol=NULL, abs_tol=NULL, npts=1000,
                       only_plot_baseline_pop=FALSE,
                       show_statistics_on_plot=TRUE,
                       line_cols=NULL, line_widths=NULL,
                       title_text=NULL, title_text_difference=NULL,
                       ISP_xlim=NULL, ISP_ylim=NULL,
                       HAP_xlim=NULL, HAP_ylim=NULL,
                       ISP_legend_xy=NULL, HAP_legend_xy=NULL,
                       axis_labels=TRUE)
{
    if (is.null(rel_tol)) rel_tol <- .Machine$double.eps^0.25
    if (is.null(abs_tol)) abs_tol <- .Machine$double.eps^0.5

    ## A note on tolerances: the 'integrate' function seems to work better if
    ## abs_tol is smaller than or equal to rel_tol, e.g.
    ##      rel_tol <- .Machine$double.eps^0.25
    ##      abs_tol <- .Machine$double.eps^0.5
    ## Often, even with rel_tol and abs_tol both equal to .Machine$double.eps^0.35
    ## (i.e. a larger rel_tol, but also a larger abs_tol), the integration fails,
    ## but it seems to be fine with a range of different powers for rel_tol,
    ## e.g. from 0.15 all the way down to 0.4, provided that abs_tol is strictly
    ## smaller e.g. power 0.5.

    npop <- length(refpopnames)

    if ((is.null(dat) | is.null(locnames)) && is.null(allele_freqs)) stop("You must provide either a data frame with allelic data and a list of loci names, or a list of allele_freqs for each locus, in the specified format.")

    if (!is.null(dat))
    {
        if (!all(refpopnames %in% unique(dat$pop))) stop("Some of the named reference pops are not in the data set.")
        if (npop < 2) stop("This function only works for more than one ref pop.")
    }

    if (is.null(line_cols)) line_cols <- c("black","blue3","forestgreen","brown","darkorchid4","gold4","dodgerblue2")[1:npop]
    if (is.null(line_widths)) line_widths <- rep(4,times=npop)
    if (npts < 500)
    {
        # npts <- 500 ## npts needs to be at least 500 to ensure the accuracy of the integral of the overlap area
        # warning("npts too small, resetting npts to 500.")
    }

    if (!is.null(quantiles_vec) & (!is.vector(quantiles_vec) || !is.numeric(quantiles_vec) ||
                                   any(quantiles_vec <0) || any(quantiles_vec > 1))) stop("quantiles_vec must be a vector of numeric values between 0 and 1 inclusive, indicating which percentiles to plot.")

    if (is.null(allele_freqs)) posterior_nu_list <- calc_posterior_nu(dat,refpopnames,locnames,prior,includepopnames=includepopnames)
    else posterior_nu_list <- allele_freqs

    ## If includepopnames is not empty, need to run the calculations part of GenePlot
    ## to calculate the logprobs for the individuals in the includepopnames
    if (!is.null(includepopnames)) logprob <- calc_logprob(dat, refpopnames, includepopnames,
                                                           locnames, prior, logten,
                                                           leave_one_out=F, saddlepoint=T,
                                                           quantiles=NULL)

    ## Create matrices to store the statistic results
    overlap_results <- matrix(rep(0,npop*npop),nrow=npop)
    incumbent_results <- matrix(rep(0,npop*npop),nrow=npop)
    home_assignment_results <- matrix(rep(0,npop*npop),nrow=npop)
    spdf_results <- list()

    if (!calc_positive_stats) print("Calculating negative statistics e.g. probability of the OTHER pop individual having higher LGP.")

    for (refpopA in refpopnames)
    {
        other_refpops <- refpopnames[-which(refpopnames == refpopA)]
        n_other_pops <- length(other_refpops)
        base_pop_idxs <- match(c(refpopA,other_refpops),refpopnames)

        spdf_vals <- calc_spdf_single_refpop_usingRcpp(posterior_nu_list=posterior_nu_list,
                                                       refpopA=refpopA,
                                                       base_pop_idxs=base_pop_idxs,
                                                       other_refpops=other_refpops,
                                                       n_other_pops=n_other_pops,
                                                       npts=npts, plot_spdfs=plot_spdfs,
                                                       calc_vecs=calc_vecs,
                                                       calc_positive_stats=calc_positive_stats,
                                                       logten=logten,
                                                       leave_one_out=leave_one_out,
                                                       quantiles_vec=quantiles_vec,
                                                       rel_tol=rel_tol,abs_tol=abs_tol,
                                                       difference_threshold=difference_threshold)

        ## Store the results for this baseline pop (popA)
        overlap_results[base_pop_idxs,base_pop_idxs[1]] <- spdf_vals$oavals ## Columns are the baseline pops, and rows are the comparison pops
        incumbent_results[base_pop_idxs,base_pop_idxs[1]] <- spdf_vals$probvals ## Columns are the baseline pops, and rows are the comparison pops
        home_assignment_results[base_pop_idxs,base_pop_idxs[1]] <- spdf_vals$diffvals ## Columns are the baseline pops, and rows are the comparison pops

        if (plot_spdfs)
        {
            plot_spdf_single_refpop_grid(spdf_vals, display_names=display_names,
                                         includepopnames=includepopnames,
                                         only_plot_baseline_pop=only_plot_baseline_pop,
                                         calc_positive_stats=calc_positive_stats,
                                         show_statistics_on_plot=show_statistics_on_plot,
                                         show_quantiles=show_quantiles,
                                         line_cols=line_cols, line_widths=line_widths,
                                         title_text=title_text,
                                         ISP_xlim=ISP_xlim, ISP_ylim=ISP_ylim,
                                         ISP_legend_xy=ISP_legend_xy, axis_labels=axis_labels)

            plot_difference_spdf_single_refpop(spdf_vals, display_names=display_names,
                                               calc_positive_stats=calc_positive_stats,
                                               show_statistics_on_plot=show_statistics_on_plot,
                                               line_cols=line_cols, line_widths=line_widths,
                                               title_text=title_text_difference,
                                               HAP_xlim=HAP_xlim, HAP_ylim=HAP_ylim,
                                               HAP_legend_xy=HAP_legend_xy, axis_labels=axis_labels)
        }

        spdf_results[[length(spdf_results)+1]] <- spdf_vals
    }

    overlap_results <- prepare_outputs(overlap_results, refpopnames, output_as_vectors=output_as_vectors)
    incumbent_results <- prepare_outputs(incumbent_results, refpopnames, output_as_vectors=output_as_vectors)
    home_assignment_results <- prepare_outputs(home_assignment_results, refpopnames, output_as_vectors=output_as_vectors)

    list(overlap_results=overlap_results,
         incumbent_results=incumbent_results,
         home_assignment_results=home_assignment_results,
         spdf_results=spdf_results)
}

prepare_outputs <- function(raw_outputs, refpopnames, output_as_vectors)
{
    if (output_as_vectors)
    {
        ## as.vector outputs all the entries in the first column, then the
        ## second column, etc. which is fine for us because we want to output
        ## all the results for the first Baseline/column, then the second
        ## Baseline/column, etc. (the comparison pops are the rows).
        ## When the "outer" results are converted using as.vector, then this
        ## also goes down the first column of outer, then the second column,
        ## etc. But that is equivalent to looping through the FIRST set of
        ## entries to "outer" within the loop of the SECOND set, but we want to
        ## loop through the SECOND entries within the FIRST, because that
        ## matches with going down the comparison pops in the first
        ## column/baseline, THEN going down the comparison pops in the second
        ## column/baseline, etc.
        final_outputs <- as.vector(raw_outputs)
        names(final_outputs) <- as.vector(outer(refpopnames,refpopnames,Vectorize(function(pop1,pop2) paste0(pop2,".",pop1))))

        ## Now remove the outputs corresponding to each population with itself
        npops <- length(refpopnames)
        pop_self_idxs <- sapply(0:(npops-1), function(idx) idx*npops + idx + 1)
        final_outputs <- final_outputs[-pop_self_idxs]
    }
    else
    {
        final_outputs <- raw_outputs
        colnames(final_outputs) <- refpopnames
        rownames(final_outputs) <- refpopnames
    }

    final_outputs
}

calc_posterior_nu <- function(dat, refpopnames, locnames, prior, includepopnames=NULL)
{
    ## Make sure that dat$pop is stored as a character vector:
    dat$pop <- as.character(dat$pop)
    if (is.null(includepopnames)) allpopnames <- unique(refpopnames)
    else allpopnames <- unique(c(refpopnames,includepopnames))

    if (!all(allpopnames %in% unique(dat$pop))) stop("Some of the pops listed are not in the dataset. Please check refpopnames and includepopnames.")

    ## Reduce the data frame "dat" to include only the populations in ref pops
    dat <- dat[!is.na(match(dat$pop, allpopnames)),]

    ## popx2 is useful so we can concatenate the vector for the first allele followed by
    ## the vector for the second allele, keeping the population information intact:
    popx2 <- rep(dat$pop, 2)

    ## Calculate the posterior allele frequencies for all populations at each locus
    ## This function returns an array with rows matching the refpopnames, and with
    ## entries being the vector (nu_1, ..., nu_k) for this row. For example the
    ## output could be (for a locus with 8 alleles):
    ##           150 154 156 164 172 178 182 184
    ## Pop1       1    1   19     6    12    1   18    10
    ## Pop2       1    1   16     1    42    1    1     1
    ## meaning that the posterior for allele frequencies for population Pop2 is
    ## (p1, ..., p8)~Dirichlet(1, 1, 16, ..., 1).
    posterior_nu_list <- lapply(locnames, calc_posterior_locus, dat, popx2, prior, refpopnames)
    names(posterior_nu_list) <- locnames

    posterior_nu_list
}

#' @importFrom stats integrate
#' @importFrom utils head tail
calc_spdf_single_refpop_usingRcpp <- function(posterior_nu_list, refpopA, base_pop_idxs,
                                              other_refpops, n_other_pops, npts,
                                              plot_spdfs=TRUE, calc_vecs=TRUE,
                                              calc_details=FALSE, calc_positive_stats=TRUE,
                                              quantiles_vec=c(0.01),
                                              logten=T, leave_one_out=T,
                                              rel_tol=NULL, abs_tol=NULL,
                                              difference_threshold=0)
{
    if (is.null(rel_tol)) rel_tol <- .Machine$double.eps^0.25
    if (is.null(abs_tol)) abs_tol <- .Machine$double.eps^0.5

    yvals <- NULL

    ## Calculate the standard distribution for population A
    postA <- lapply(posterior_nu_list, function(x) x[refpopA,, drop=FALSE])

    postA_info <- rcpp_calc_multi_locus_dist(postA, leave_one_out=leave_one_out, differenceLGPs=FALSE)

    min_postA <- postA_info$min
    max_postA <- postA_info$max
    mean_postA <- rcpp_calc_mu(postA_info$dist)
    loge_mean_postA <- mean_postA
    if (logten) mean_postA <- mean_postA/log(10)

    ## If min_postA is infinite, need to correct it to a really large number
    if (is.infinite(min_postA)) min_postA <- 1E-6

    ## Find points along the distribution, remove the end points where
    ## f_hat is typically NA and calculate f_hat at each point, and normalize
    xvals <- seq(min_postA, max_postA, length=npts)
    xvals <- xvals[-c(1,npts)]
    if (logten) xvals <- xvals/log(10)

    fhA_area <- tryCatch(
        integrate(rcpp_calc_fhat, lower=xvals[1], upper=xvals[length(xvals)],
                  postA_info$dist, min_postA, max_postA,
                  loge_mean_postA, logten=logten, rel.tol=rel_tol, abs.tol=abs_tol,
                  subdivisions=npts)$value,
        error=function(err){
            message(paste("My error: ",err))
            # browser()
        })

    if (plot_spdfs | calc_vecs) {
        fhA <- rcpp_calc_fhat(xvals, postA_info$dist, postA_info$min, postA_info$max,
                              loge_mean_postA, logten=logten)

        fhA_norm <- fhA/fhA_area
        yvals <- fhA_norm
    } else yvals <- rep(NA,times=npts-2)

    if (calc_details)
    {
        var_postA <- rcpp_calc_multi_locus_K2(postA_info$dist, 0)
        skew_postA <- rcpp_calc_multi_locus_K3(postA_info$dist, 0)/(rcpp_calc_multi_locus_K2(postA_info$dist, 0)^1.5)

        if (logten)
        {
            min_postA <- min_postA/log(10)
            max_postA <- max_postA/log(10)
            var_postA <- var_postA/(log(10)^2)
            skew_postA <- skew_postA/(log(10)^3)
        }
        range_postA <- max_postA - min_postA
        Fh_at_mean_postA <- Fhat(mean_postA,postA_info$dist, postA_info$min, postA_info$max,
                                 loge_mean_postA,logten=logten)

        print(paste0("Summary statistics for the standard LGP distribution of ",refpopA))
        print(paste0("Mean: ",mean_postA))
        print(paste0("Fh at mean: ",Fh_at_mean_postA))
        print(paste0("Variance: ",var_postA))
        print(paste0("Normalization constant: ",fhA_area))
    }

    if (is.null(quantiles_vec)) quantiles_vec <- c(0.01)
    quantiles <- matrix(rep(NA,length(quantiles_vec)), nrow=1)
    spdf_qsearch_params_pop <- rcpp_calc_qsearch_params(postA_info$dist, postA_info$min,
                                                        postA_info$max, logten=logten)
    for (i in 1:length(quantiles_vec))
    {
        quantiles[1,i] <- rcpp_calc_Qhat(quantiles_vec[i], spdf_qsearch_params_pop,
                                         postA, logten=logten)
    }
    rownames(quantiles) <- refpopA
    attributes(quantiles)$quantiles_vec <- quantiles_vec

    oavals <- 1
    probvals <- 1
    diffvals <- 1

    ## Loop through all other pops ---------------------------------------------
    for (refpopB in other_refpops)
    {
        ## Calculate the distribution with values corresponding to LGPs from
        ## one pop and probabilities corresponding to the probabilities of LGPs
        ## in the other pop. For example, the individuals with genotypes that occur
        ## commonly in pop A may have very low LGPs with respect to pop B
        ## and this would indicate that a typical individual taken from pop A
        ## would be expected to have a poor fit to pop B.
        ## The resulting distribution will have probs from B and values being LGPs
        ## from pop A.

        post_both_pops <- lapply(posterior_nu_list, function(x) rbind(x[refpopA,],x[refpopB,]))

        ## This is a combined distribution from both populations so CANNOT be leave-one-out
        post_both_pops_info <- rcpp_calc_multi_locus_dist(post_both_pops, leave_one_out=FALSE, differenceLGPs=FALSE)

        mean_post_both_pops <- rcpp_calc_mu(post_both_pops_info$dist)
        loge_mean_post_both_pops <- mean_post_both_pops
        if (logten) mean_post_both_pops <- mean_post_both_pops/log(10)

        fh_both_pops_area <- integrate(rcpp_calc_fhat,
                                       lower=head(xvals,1), upper=tail(xvals,1),
                                       post_both_pops_info$dist, post_both_pops_info$min,
                                       post_both_pops_info$max, loge_mean_post_both_pops,
                                       logten=logten, rel.tol=rel_tol, abs.tol=abs_tol)$value

        if (plot_spdfs | calc_vecs) {
            fh_both_pops <- rcpp_calc_fhat(xvals, post_both_pops_info$dist,
                                           post_both_pops_info$min, post_both_pops_info$max,
                                           loge_mean_post_both_pops, logten=logten)
            fh_both_pops_norm <- fh_both_pops/fh_both_pops_area
            yvals <- rbind(yvals, fh_both_pops_norm)
        } else yvals <- rbind(yvals, rep(NA,npts-2))

        print(paste0("refpopA = ",refpopA,", refpopB = ",refpopB))

        ## Calculate the area of overlap between the 2 distributions, the
        ## baseline pop distribution and the combined bothPops distribution.
        overlap_function <- function(x) min(rcpp_calc_fhat(x,postA_info$dist, postA_info$min,
                                                           postA_info$max, loge_mean_postA,
                                                           logten=logten)/fhA_area,
                                            rcpp_calc_fhat(x,post_both_pops_info$dist,
                                                           post_both_pops_info$min,
                                                           post_both_pops_info$max,
                                                           loge_mean_post_both_pops,
                                                           logten=logten)/fh_both_pops_area)
        overlap_area <- integrate(Vectorize(overlap_function),
                                  lower=xvals[1], upper=xvals[length(xvals)],
                                  rel.tol=rel_tol, abs.tol=abs_tol)$value
        # rel.tol=rel_tol, abs.tol=abs_tol, subdivisions=npts)$value
        ### NOTE: Overlap Area is always calculated the same way even under calc_positive_stats=FALSE
        oavals <- rbind(oavals, overlap_area)

        ## Calculate the integral of the PDF of bothPops distribution multiplied
        ## by the CDF of the baseline pop distribution. This integral gives the
        ## probability that for a given pair of genotypes from popA and popB,
        ## the LGP with respect to popA will be higher for the genotype from popB.
        incumbent_function <- function(x) rcpp_calc_fhat(x,post_both_pops_info$dist,
                                                         post_both_pops_info$min,
                                                         post_both_pops_info$max,
                                                         loge_mean_post_both_pops,
                                                         logten=logten)*
            rcpp_calc_Fhat(x,postA_info$dist,
                           postA_info$min, postA_info$max,
                           loge_mean_postA,logten=logten)

        incumbent_value <- integrate(Vectorize(incumbent_function),
                                     lower=xvals[1], upper=xvals[length(xvals)],
                                     rel.tol=rel_tol, abs.tol=abs_tol)$value/fh_both_pops_area
        # rel.tol=rel_tol, abs.tol=abs_tol, subdivisions=npts)$value/fh_both_pops_area
        if (calc_positive_stats) incumbent_value <- 1-incumbent_value

        probvals <- rbind(probvals, incumbent_value)

        if (calc_positive_stats) {
            print(paste0("Overlap area: of refpopA SPDF and SPDF of refpopB into refPopA = ",round(overlap_area,6)))
            print(paste0("ISP: Prob. of refpopA LGP being greater than refpopB LGP = ",round(incumbent_value,6)))
        } else {
            print(paste0("Overlap area: of refPopA SPDF and SPDF of refpopB into refpopA = ",round(overlap_area,6)))
            print(paste0("Inverse-ISP: Prob. of refpopB LGP being greater than refpopA LGP = ",round(incumbent_value,6)))
        }

        ## Now calculate the SPDF for the differences between the LGPs with
        ## respect to the two populations, which produces a distribution
        ## over wvals not xvals

        post_diff_pops_info <- rcpp_calc_multi_locus_dist(post_both_pops,
                                                          leave_one_out=leave_one_out,
                                                          differenceLGPs=TRUE)

        mean_post_diff_pops <- rcpp_calc_mu(post_diff_pops_info$dist)
        loge_mean_post_diff_pops <- mean_post_diff_pops
        if (logten) mean_post_diff_pops <- mean_post_diff_pops/log(10)

        ## Find points along the distribution, remove the end points where
        ## f_hat is typically NA and calculate f_hat at each point, and normalize
        wvals <- seq(post_diff_pops_info$min, post_diff_pops_info$max, length=npts)
        wvals <- wvals[-c(1,npts)]
        if (logten) wvals <- wvals/log(10)

        fh_diff_pops_area <- integrate(rcpp_calc_fhat,
                                       lower=head(wvals,1), upper=tail(wvals,1),
                                       post_diff_pops_info$dist, post_diff_pops_info$min,
                                       post_diff_pops_info$max, loge_mean_post_diff_pops,
                                       logten=logten, rel.tol=rel_tol, abs.tol=abs_tol)$value

        if (plot_spdfs | calc_vecs) {
            fh_diff_pops <- rcpp_calc_fhat(wvals, post_diff_pops_info$dist,
                                           post_diff_pops_info$min,
                                           post_diff_pops_info$max,
                                           loge_mean_post_diff_pops, logten=logten)
            fh_diff_popsNorm <- fh_diff_pops/fh_diff_pops_area
            if (refpopB == other_refpops[1]) zvals <- fh_diff_popsNorm
            else zvals <- rbind(zvals, fh_diff_popsNorm)
        } else zvals <- rep(NA,npts-2)

        home_assignment_value <- rcpp_calc_Fhat(difference_threshold,
                                                post_diff_pops_info$dist,
                                                post_diff_pops_info$min,
                                                post_diff_pops_info$max,
                                                loge_mean_post_diff_pops,
                                                logten=logten)
        if (calc_positive_stats) home_assignment_value <- 1-home_assignment_value

        diffvals <- rbind(diffvals, home_assignment_value)

        if (difference_threshold==0) {
            if (calc_positive_stats) print(paste0("HAP: Prob. of LGP for random refpopA genotype being greater w.r.t. refpopA = ",round(home_assignment_value,6)))
            else print(paste0("Inverse-HAP: Prob. of LGP for random refpopA genotype being greater w.r.t. refpopB = ",round(home_assignment_value,6)))
        } else {
            if (logten) timesGreater <- 10^difference_threshold
            else timesGreater <- exp(difference_threshold)

            if (calc_positive_stats) print(paste0("HAP: Prob. of LGP for random refpopA genotype being ",timesGreater," times greater w.r.t. refpopA = ",round(home_assignment_value,6)))
            else print(paste0("Inverse-HAP: Prob. of LGP for random refpopA genotype being ",timesGreater," times greater w.r.t. refpopB = ",round(home_assignment_value,6)))
        }

        print("")

        ## Calc details for other pops -----------------------------------------
        if (calc_details)
        {
            min_post_both_pops <- post_both_pops_info$min
            max_post_both_pops <- post_both_pops_info$max

            var_post_both_pops <- rcpp_calc_multi_locus_K2(post_both_pops_info$dist, 0)
            skew_post_both_pops <- rcpp_calc_multi_locus_K3(post_both_pops_info$dist, 0)/
                (rcpp_calc_multi_locus_K2(post_both_pops_info$dist, 0)^1.5)

            Fh_at_mean_post_both_pops <- rcpp_calc_Fhat(mean_post_both_pops,
                                                        post_both_pops_info$dist,
                                                        post_both_pops_info$min,
                                                        post_both_pops_info$max,
                                                        loge_mean_post_both_pops,
                                                        logten=logten)

            ## Calculate the CDF value for the SPDF distribution at the mean of the SDPF* distribution
            FhAt_both_popsMean <- Fhat(mean_post_both_pops, postA_info$dist,
                                       postA_info$min, postA_info$max,
                                       loge_mean_postA,logten=logten)

            if (logten)
            {
                min_post_both_pops <- min_post_both_pops/log(10)
                max_post_both_pops <- max_post_both_pops/log(10)
                mean_post_both_pops <- loge_mean_post_both_pops/log(10)
                var_post_both_pops <- var_post_both_pops/(log(10)^2)
                skew_post_both_pops <- skew_post_both_pops/(log(10)^3)
            }

            range_post_both_pops <- max_post_both_pops - min_post_both_pops

            print(paste0("Summary statistics for the LGP distribution of ",refpopA," in ",refpopB))
            print(paste0("Mean: ",mean_post_both_pops))
            print(paste0("Fh at mean: ",Fh_at_mean_post_both_pops))
            print(paste0("Variance: ",var_post_both_pops))
            print(paste0("Normalization constant: ",fh_both_pops_area))
            print(paste0("Difference of means of SPDF and SPDF** = ",round(mean_postA-mean_post_both_pops,4)))
            print(paste0("Log10 of SCDF percentile at SPDF** mean = ",log10(FhAt_both_popsMean)))
            print("")
        }
    }

    if (is.vector(zvals)) zvals <- matrix(zvals, ncol=length(wvals))

    list(xvals=xvals, yvals=yvals, wvals=wvals, zvals=zvals,
         oavals=oavals, probvals=probvals, diffvals=diffvals,
         quantiles=quantiles, refpopA=refpopA, base_pop_idxs=base_pop_idxs,
         other_refpops=other_refpops, n_other_pops=n_other_pops)
}

plot_difference_spdf_single_refpop <- function(spdf_vals, display_names,
                                               calc_positive_stats=TRUE,
                                               show_statistics_on_plot=TRUE,
                                               line_cols=NULL, line_widths=NULL,
                                               title_text=NULL, axis_labels=TRUE,
                                               HAP_xlim=NULL, HAP_ylim=NULL,
                                               HAP_legend_xy=NULL, show_legend=TRUE)
{
    base_pop_idxs <- spdf_vals$base_pop_idxs
    n_other_pops <- spdf_vals$n_other_pops

    ### Setup plot params ------------------------------------------------------
    if (is.null(HAP_ylim))
    {
        zvals_vec <- as.vector(spdf_vals$zvals)
        zlim_actual <- range(zvals_vec[which(!is.nan(zvals_vec) & !is.infinite(zvals_vec))])
    }
    else zlim_actual <- HAP_ylim
    if (is.null(HAP_xlim)) wlim_actual <- range(spdf_vals$wvals[is.finite(spdf_vals$wvals)])
    else wlim_actual <- HAP_xlim
    # wlim_actual <- c(-0.1,0.1)

    if (show_statistics_on_plot)
    {
        plot_matrix <- matrix(1:(n_other_pops+1),ncol=1)
        layout(plot_matrix, widths=1, heights=c(8,rep(1,n_other_pops)))
    }

    if (axis_labels)
    {
        wlab <- paste('Difference of log-probabilities for multilocus genotypes')
        zlab <- paste('Probability density within',display_names[base_pop_idxs[1]],'dist')
    }
    else
    {
        wlab <- ""
        zlab <- ""
    }

    if (is.null(line_cols)) line_cols <- grDevices::rainbow(n_other_pops+1, s=0.5, start=0.625, end=0.42)
    if (is.null(line_widths)) line_widths <- rep(4,(n_other_pops+1))

    ## Create main plot --------------------------------------------------------
    par(mar=c(2, 2, 2, 2) + 0.1)
    plot(spdf_vals$wvals, spdf_vals$zvals[1,], type="l", lty=1, col=line_cols[base_pop_idxs[2]],
         lwd=line_widths[base_pop_idxs[1]], cex.lab=1.2, cex.axis=1.2,
         xlim=wlim_actual, ylim=zlim_actual,xlab=wlab,ylab=zlab)
    if (is.null(title_text)) graphics::title(paste('LGPs in',display_names[base_pop_idxs[1]],'pop minus LGPs in other pops'))
    else graphics::title(title_text)

    graphics::polygon(c(min(spdf_vals$wvals), spdf_vals$wvals, max(spdf_vals$wvals)),
                      c(min(spdf_vals$zvals[1,]), spdf_vals$zvals[1,], min(spdf_vals$zvals[1,])),
                      border=NA, col=grDevices::adjustcolor(line_cols[base_pop_idxs[2]],alpha.f=0.3))

    ## Plot other pops ---------------------------------------------------------
    if (is.null(HAP_legend_xy)) HAP_legend_xy <- c(spdf_vals$wvals[1]+3, 0.25)
    if (n_other_pops > 1) {
        for (pop_idx in 2:n_other_pops)
        {
            # first row of yvals is popA in popA, so ignore that one
            graphics::lines(spdf_vals$wvals,spdf_vals$zvals[pop_idx,],
                            col=line_cols[base_pop_idxs[pop_idx+1]],
                            lwd=line_widths[base_pop_idxs[pop_idx+1]])

            graphics::polygon(c(min(spdf_vals$wvals), spdf_vals$wvals, max(spdf_vals$wvals)),
                              c(min(spdf_vals$zvals[pop_idx,]), spdf_vals$zvals[pop_idx,], min(spdf_vals$zvals[pop_idx,])),
                              border=NA, col=grDevices::adjustcolor(line_cols[base_pop_idxs[pop_idx+1]],alpha.f=0.3))
        }
    }
    if (show_legend) legend(HAP_legend_xy[1],HAP_legend_xy[2], display_names[base_pop_idxs], cex=1.3,
                            lty=rep(1,times=n_other_pops), lwd=line_widths[base_pop_idxs],
                            col=line_cols[base_pop_idxs])

    ## Show proportion of the graph that is below/above 0 ----------------------
    if (show_statistics_on_plot)
    {
        ### Home assignment probabilities --------------------------------------
        for (pop_idx in 1:n_other_pops)
        {
            par(mar=c(1,1,1,1))
            graphics::plot.new()
            if (calc_positive_stats) mtext(paste("Probability of home assignment vs.",
                                                 display_names[base_pop_idxs[pop_idx+1]],"=",
                                                 round(spdf_vals$diffvals[pop_idx+1],6)), side=3, font=2)
            else mtext(paste("Probability of misclassification vs.",
                             display_names[base_pop_idxs[pop_idx+1]],"=",
                             round(spdf_vals$diffvals[pop_idx+1],6)), side=3, font=2)
        }
    }

    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
}

plot_spdf_single_refpop <- function(spdf_vals, display_names, includepopnames=NULL,
                                    only_plot_baseline_pop=FALSE,
                                    calc_positive_stats=TRUE,
                                    show_statistics_on_plot=TRUE,
                                    show_quantiles=TRUE, title_text=NULL,
                                    line_cols=NULL, line_widths=NULL,
                                    ISP_xlim=NULL, ISP_ylim=NULL,
                                    ISP_legend_xy=NULL, axis_labels=TRUE)
{
    base_pop_idxs <- spdf_vals$base_pop_idxs
    n_other_pops <- spdf_vals$n_other_pops

    ### Setup plot params ------------------------------------------------------
    if (is.null(ISP_ylim))
    {
        yvals_vec <- as.vector(spdf_vals$yvals)
        ylim_actual <- range(yvals_vec[which(!is.nan(yvals_vec) & !is.infinite(yvals_vec))])
    }
    else ylim_actual <- ISP_ylim
    if (is.null(ISP_xlim)) xlim_actual <- range(spdf_vals$xvals[is.finite(spdf_vals$xvals)])
    else xlim_actual <- ISP_xlim

    if (show_statistics_on_plot)
    {
        plot_matrix <- matrix(1:(n_other_pops+1),ncol=1)
        layout(plot_matrix, widths=1, heights=c(8,rep(1,n_other_pops)))
    }

    if (axis_labels)
    {
        xlab <- paste('Log-probabilities for multilocus genotypes')
        ylab <- paste('Probability density within',display_names[base_pop_idxs[1]],'dist')
    }
    else
    {
        xlab <- ""
        ylab <- ""
    }

    if (is.null(line_cols)) line_cols <- grDevices::rainbow(n_other_pops+1, s=0.5, start=0.625, end=0.42)
    if (is.null(line_widths)) line_widths <- rep(4,(n_other_pops+1))

    ## Create main plot --------------------------------------------------------
    par(mar=c(2, 2, 2, 2) + 0.1)
    plot(spdf_vals$xvals, spdf_vals$yvals[1,], type="l", lty=2, col=line_cols[base_pop_idxs[1]],
         lwd=line_widths[base_pop_idxs[1]], cex.lab=1.3, cex.axis=1.3,
         xlim=xlim_actual, ylim=ylim_actual,xlab=xlab,ylab=ylab)
    if (is.null(title_text)) {
        if (only_plot_baseline_pop) graphics::title(paste('Log-probabilities in',display_names[base_pop_idxs[1]],'distribution'))
        else if (calc_positive_stats) graphics::title(paste('Incumbent Selection Probability in',display_names[base_pop_idxs[1]],'distribution'))
        else graphics::title(paste('Inverse Incumbent Selection Prob. in',display_names[base_pop_idxs[1]],'distribution'))
    }
    else graphics::title(title_text)

    graphics::polygon(c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                      c(min(spdf_vals$yvals[1,]), spdf_vals$yvals[1,], min(spdf_vals$yvals[1,])),
                      border=NA, col=grDevices::adjustcolor(line_cols[base_pop_idxs[1]],alpha.f=0.3))

    ## Plot the individuals from includepopnames, if required
    if (!is.null(includepopnames))
    {
        include_logprobs <- spdf_vals$logprob[which(spdf_vals$logprob$pop %in% includepopnames),spdf_vals$refpopA]
        points(include_logprobs, rep(0,length(include_logprobs)), pch=18, cex=2)
    }

    if (show_quantiles)
    {
        ## Add column headings to the quantile table
        quantiles_vec <- attributes(spdf_vals$quantiles)$quantiles_vec
        for(idx in 1:length(quantiles_vec)) abline(v=spdf_vals$quantiles[spdf_vals$refpopA,idx],lwd=2,lty=2)
        axis(side=1, at=spdf_vals$quantiles[spdf_vals$refpopA,idx],
             labels=paste0(round(100*quantiles_vec,0),"%"),
             col.axis="black", tick=F, cex.axis=1, line=-0.8)
    }

    ## Plot other pops ---------------------------------------------------------
    if (is.null(ISP_legend_xy)) ISP_legend_xy <- c(spdf_vals$xvals[1]+3, 0.25)
    if (only_plot_baseline_pop)
    {
        legend(ISP_legend_xy[1],ISP_legend_xy[2], display_names[base_pop_idxs[1]], cex=1.3, lty=2, lwd=line_widths[base_pop_idxs[1]], col=line_cols[base_pop_idxs[1]])
    }
    else
    {
        for (pop_idx in 1:n_other_pops)
        {
            # first row of yvals is popA in popA, so ignore that one
            graphics::lines(spdf_vals$xvals,spdf_vals$yvals[pop_idx+1,],
                            col=line_cols[base_pop_idxs[pop_idx+1]],
                            lwd=line_widths[base_pop_idxs[pop_idx+1]])

            graphics::polygon(c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                              c(min(spdf_vals$yvals[pop_idx+1,]), spdf_vals$yvals[pop_idx+1,], min(spdf_vals$yvals[pop_idx+1,])),
                              border=NA, col=grDevices::adjustcolor(line_cols[base_pop_idxs[pop_idx+1]],alpha.f=0.3))
        }
        legend(ISP_legend_xy[1],ISP_legend_xy[2], display_names[base_pop_idxs], cex=1.3,
               lty=c(2,rep(1,times=n_other_pops)),
               lwd=line_widths[base_pop_idxs], col=line_cols[base_pop_idxs])
    }

    if (show_statistics_on_plot)
    {
        ### Incumbent selection probabilities ----------------------------------
        for (pop_idx in 1:n_other_pops)
        {
            par(mar=c(1,1,1,1))
            graphics::plot.new()
            mtext(paste("Incumbent Selection Prob for",
                        display_names[base_pop_idxs[pop_idx+1]],"=",
                        round(spdf_vals$probvals[pop_idx+1],6)), side=3, font=2)
        }

        par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    }
}

plot_spdf_single_refpop_grid <- function(spdf_vals, display_names,
                                         includepopnames=NULL,
                                         show_title=TRUE, title_text=NULL,
                                         only_plot_baseline_pop=FALSE,
                                         calc_positive_stats=TRUE,
                                         show_statistics_on_plot=TRUE,
                                         show_quantiles=TRUE,
                                         line_cols=NULL, line_widths=NULL,
                                         ISP_xlim=NULL, ISP_ylim=NULL,
                                         ISP_legend_xy=NULL, axis_labels=TRUE)
{
    ### Setup plot params ------------------------------------------------------
    if (is.null(spdf_vals$xvals) || is.null(spdf_vals$yvals)) stop("Missing inputs. Need xvals and yvals within spdf_vals object.")

    base_pop_idxs <- spdf_vals$base_pop_idxs
    n_other_pops <- spdf_vals$n_other_pops

    if (is.null(display_names)) display_names <- c(spdf_vals$refpopA, spdf_vals$otherRefPops)

    if (is.null(ISP_ylim))
    {
        yvals_vec <- as.vector(spdf_vals$yvals)
        ylim_actual <- range(yvals_vec[which(!is.nan(yvals_vec) & !is.infinite(yvals_vec))])
        ylim_actual[1] <- 0
    }
    else ylim_actual <- ISP_ylim
    if (is.null(ISP_xlim)) xlim_actual <- range(spdf_vals$xvals[is.finite(spdf_vals$xvals)])
    else xlim_actual <- ISP_xlim

    if (is.null(line_widths)) line_widths <- rep(4,times=(1+n_other_pops))
    if (is.null(line_cols)) line_cols <- grDevices::rainbow(1+n_other_pops, s=0.5, start=0.625, end=0.42)

    if (axis_labels)
    {
        xlab <- paste('Log-probabilities for multilocus genotypes')
        ylab <- paste('Probability density within',display_names[base_pop_idxs[1]],'dist')
    } else {
        xlab <- ""
        ylab <- ""
    }

    ### Setup overall plot -----------------------------------------------------

    element_height_units <- c("cm","null","cm",rep("cm",n_other_pops))
    element_heights <- c(show_title,1,1.5, rep(1*show_statistics_on_plot,n_other_pops))
    element_height_units <- element_height_units[(element_heights != 0)]
    element_heights <- element_heights[(element_heights != 0)]

    ## Now convert the central plot height to 1 if all the rest of the elements
    ## are missing (i.e. not shown)
    if (length(element_heights) == 1) element_heights <- 1

    lay <- grid::grid.layout(nrow=(show_title+1+1+show_statistics_on_plot*n_other_pops), ncol=3,
                             widths=grid::unit(c(2,1,0.5),c("cm","null","cm")),
                             heights=grid::unit(element_heights,element_height_units))

    vp1 <- grid::viewport(layout = lay)

    grid::grid.newpage()
    grid::pushViewport(vp1)

    ### Part 1: the title ------------------------------------------------------
    if (show_title){
        if (is.null(title_text) && only_plot_baseline_pop) title_text <- paste('Log-probabilities in',display_names[base_pop_idxs[1]],'distribution')
        else if (is.null(title_text)) {
            if (only_plot_baseline_pop) title_text <- paste('Log-probabilities in',display_names[base_pop_idxs[1]],'distribution')
            else if (calc_positive_stats) title_text <- paste('Incumbent Selection Probability in',display_names[base_pop_idxs[1]],'distribution')
            else title_text <- paste('Inverse Incumbent Selection Prob. in',display_names[base_pop_idxs[1]],'distribution')
        }
        vp2 <- grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3)
        grid::pushViewport(vp2)
        grid::grid.text(title_text,gp=grid::gpar(fontface="bold"))
        grid::popViewport()  ## pop the title viewport i.e. navigate back up to the level of the overall viewport
    }

    ### Part 2: the y-axis label -----------------------------------------------
    vp3 <- grid::viewport(layout.pos.row = show_title+1, layout.pos.col = 1)
    grid::pushViewport(vp3)

    # grid::grid.yaxis(main=FALSE)
    grid::grid.text(ylab, x = grid::unit(0.25, "cm"), rot = 90)

    grid::popViewport()  ## pop the title viewport i.e. navigate back up to the level of the overall viewport

    ### Part 3: the main plot --------------------------------------------------
    vp4 <- grid::viewport(layout.pos.row = show_title+1, layout.pos.col = 2,
                          xscale = xlim_actual, yscale = ylim_actual, gp=grid::gpar(cex=1))
    grid::pushViewport(vp4)

    grid::grid.yaxis(main=TRUE)
    grid::grid.xaxis(main=TRUE)

    grid::grid.lines(spdf_vals$xvals, spdf_vals$yvals[1,], default.units="native",
                     gp=grid::gpar(col=line_cols[base_pop_idxs[1]], lty=2, cex=1.3,
                             lwd=line_widths[base_pop_idxs[1]]))

    grid::grid.polygon(c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                       c(min(spdf_vals$yvals[1,]), spdf_vals$yvals[1,], min(spdf_vals$yvals[1,])),
                       default.units="native",
                       gp=grid::gpar(col=NA, fill=line_cols[base_pop_idxs[1]], alpha=0.3))

    ## Plot the individuals from includepopnames, if required
    if (!is.null(includepopnames))
    {
        include_logprobs <- spdf_vals$logprob[which(spdf_vals$logprob$pop %in% includepopnames),spdf_vals$refpopA]
        grid::grid.points(include_logprobs, rep(0,length(include_logprobs)),
                          pch=18, gp=grid::gpar(cex=2))
    }

    if (show_quantiles)
    {
        ## Add column headings to the quantile table
        quantiles_vec <- attributes(spdf_vals$quantiles)$quantiles_vec
        for(idx in 1:length(quantiles_vec))
        {
            grid::grid.function(function(z) list(x=spdf_vals$quantiles[spdf_vals$refpopA,idx],y=z),
                                range="y", gp=grid::gpar(col="black", lty=2, lwd=2))
            grid::grid.text(paste0(round(100*quantiles_vec[idx],0),"%"),
                            x=grid::unit(spdf_vals$quantiles[spdf_vals$refpopA,idx],"native"),
                            y=grid::unit(-0.5, "lines"), gp=grid::gpar(col="black"))
        }
    }

    ### Part 3b: the other pops ------------------------------------------------
    if (is.null(ISP_legend_xy)) ISP_legend_xy <- c(spdf_vals$xvals[1]+3, 0.25)
    if (only_plot_baseline_pop)
    {
        # legend(ISP_legend_xy[1],ISP_legend_xy[2], display_names[base_pop_idxs[1]], cex=1.3, lty=2, lwd=line_widths[base_pop_idxs[1]], col=line_cols[base_pop_idxs[1]])
    }
    else
    {
        for (pop_idx in 1:n_other_pops)
        {
            # first row of yvals is popA in popA, so ignore that one
            grid::grid.lines(spdf_vals$xvals,spdf_vals$yvals[pop_idx+1,], default.units="native",
                             gp=grid::gpar(col=line_cols[base_pop_idxs[pop_idx+1]], cex=1.3,
                                     lwd=line_widths[base_pop_idxs[pop_idx+1]]))

            grid::grid.polygon(c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                               c(min(spdf_vals$yvals[pop_idx+1,]), spdf_vals$yvals[pop_idx+1,],
                                 min(spdf_vals$yvals[pop_idx+1,])),default.units="native",
                               gp=grid::gpar(col=NA, fill=line_cols[base_pop_idxs[pop_idx+1]],alpha=0.3))
        }
        # legend(ISP_legend_xy[1],ISP_legend_xy[2], display_names[base_pop_idxs], cex=1.3,
        #        lty=c(2,rep(1,times=n_other_pops)),
        #        lwd=line_widths[base_pop_idxs], col=line_cols[base_pop_idxs])
    }

    grid::popViewport()  ## pop the spdf plot viewport i.e. navigate back up to the level of the overall viewport

    ### Part 4: the x-axis label -----------------------------------------------
    vp5 <- grid::viewport(layout.pos.row = show_title+1+1, layout.pos.col = 2)
    grid::pushViewport(vp5)

    grid::grid.text(xlab, y = grid::unit(0.25, "cm"))

    grid::popViewport()  ## pop the spdf plot viewport i.e. navigate back up to the level of the overall viewport

    ### Part 5: the numeric measure values -------------------------------------
    if (show_statistics_on_plot)
    {
        for (pop_idx in 1:n_other_pops)
        {
            grid::grid.text(paste("Incumbent Selection Prob for", display_names[base_pop_idxs[pop_idx+1]],"=",
                                  round(spdf_vals$probvals[pop_idx+1],6)),
                            y=grid::unit(1,"lines"),
                            vp=grid::viewport(layout.pos.row = (show_title+1+1+pop_idx),
                                              layout.pos.col = 2),
                            gp=grid::gpar(fontsize=12, fontface="bold"))

        }
    }

    grid::popViewport()  ## pop the spdf plot viewport i.e. navigate back up to the level of the overall viewport
}

plot_spdf_single_refpop_grid_rotate <- function(spdf_vals, display_names,
                                                includepopnames=NULL,
                                                only_plot_baseline_pop=FALSE,
                                                show_quantiles=TRUE,
                                                heightMultiplier=9,
                                                line_cols=NULL, line_widths=NULL,
                                                ISP_xlim=NULL, ISP_ylim=NULL,
                                                lgp_axis=FALSE,
                                                prob_density_axis=TRUE)
{
    ### Setup plot params ------------------------------------------------------
    if (is.null(spdf_vals$xvals) || is.null(spdf_vals$yvals)) stop("Missing inputs. Need xvals and yvals within spdf_vals object.")

    base_pop_idxs <- spdf_vals$base_pop_idxs
    n_other_pops <- spdf_vals$n_other_pops

    if (is.null(display_names)) display_names <- c(spdf_vals$refpopA, spdf_vals$otherRefPops)

    if (is.null(ISP_ylim))
    {
        yvals_vec <- as.vector(spdf_vals$yvals)
        ylim_actual <- range(yvals_vec[which(!is.nan(yvals_vec) & !is.infinite(yvals_vec))])
        ylim_actual[1] <- 0
    }
    else ylim_actual <- ISP_ylim
    if (is.null(ISP_xlim)) xlim_actual <- range(spdf_vals$xvals[is.finite(spdf_vals$xvals)])
    else xlim_actual <- ISP_xlim

    if (is.null(line_widths)) line_widths <- rep(4,times=(1+n_other_pops))
    if (is.null(line_cols)) line_cols <- grDevices::rainbow(1+n_other_pops, s=0.5, start=0.625, end=0.42)

    if (lgp_axis) ylab <- paste('Log-probabilities for multilocus genotypes')
    else ylab <- ""

    if (prob_density_axis) xlab <- paste('Probability density within',display_names[base_pop_idxs[1]],'dist')
    else xlab <- ""

    ### Setup overall plot -----------------------------------------------------

    if (lgp_axis)
    {
        element_width_units <- c("null","cm")
        element_widths <- c(1,2)
    } else {
        element_width_units <- "null"
        element_widths <- 1
    }

    if (prob_density_axis)
    {
        element_height_units <- c("null","cm")
        element_heights <- c(1,2)
    } else {
        element_height_units <- "null"
        element_heights <- 1
    }

    lay <- grid::grid.layout(nrow=length(element_heights), ncol=length(element_widths),
                             heights=grid::unit(element_heights,element_height_units),
                             widths=grid::unit(element_widths, element_width_units))

    vp1 <- grid::viewport(width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), layout = lay)

    grid::grid.newpage()
    grid::pushViewport(vp1)

    ### Part 1: the main plot --------------------------------------------------
    ### Note that the scale for the x-axis needs to be reversed so that the base
    ### of the plot (i.e. the 0 point on the prob. density axis) is on the RIGHT
    vp2 <- grid::viewport(layout.pos.row = 1, layout.pos.col = 1,
                          xscale = c(ylim_actual[2],ylim_actual[1]), yscale = xlim_actual, gp=grid::gpar(cex=1))
    grid::pushViewport(vp2)

    ## Note that what would normally be the x-axis (ie. the LGP axis) is being
    ## drawn as the y-axis, and what would normally be the y-axis is being drawn
    ## as the x-axis
    grid::grid.xaxis(main=TRUE)
    if (lgp_axis) grid::grid.yaxis(main=FALSE)

    grid::grid.lines(spdf_vals$yvals[1,], spdf_vals$xvals, default.units="native",
                     gp=grid::gpar(col=line_cols[base_pop_idxs[1]], lty=2, cex=1.3,
                             lwd=line_widths[base_pop_idxs[1]]))

    grid::grid.polygon(c(min(spdf_vals$yvals[1,]), spdf_vals$yvals[1,], min(spdf_vals$yvals[1,])),
                       c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                       default.units="native",
                       gp=grid::gpar(col=NA, fill=line_cols[base_pop_idxs[1]], alpha=0.3))

    ## Plot the individuals from includepopnames, if required
    if (!is.null(includepopnames))
    {
        include_logprobs <- spdf_vals$logprob[which(spdf_vals$logprob$pop %in% includepopnames),spdf_vals$refpopA]
        grid::grid.points(rep(0,length(include_logprobs)), include_logprobs,
                          pch=18, gp=grid::gpar(cex=2))
    }

    if (show_quantiles)
    {
        ## Add column headings to the quantile table
        quantiles_vec <- attributes(spdf_vals$quantiles)$quantiles_vec
        for(idx in 1:length(quantiles_vec))
        {
            grid::grid.function(function(z) list(x=z, y=spdf_vals$quantiles[spdf_vals$refpopA,idx]),
                                range="x", gp=grid::gpar(col="black", lty=2, lwd=2))
            grid::grid.text(paste0(round(100*quantiles_vec[idx],0),"%"),
                            x=grid::unit(0.5, "lines"),
                            y=grid::unit(spdf_vals$quantiles[spdf_vals$refpopA,],"native"),
                            gp=grid::gpar(col="black"))
        }
    }

    if (!only_plot_baseline_pop)
    {
        for (pop_idx in 1:n_other_pops)
        {
            # first row of yvals is popA in popA, so ignore that one
            grid::grid.lines(spdf_vals$yvals[pop_idx+1,], spdf_vals$xvals,default.units="native",
                             gp=grid::gpar(col=line_cols[base_pop_idxs[pop_idx+1]], cex=1.3,
                                     lwd=line_widths[base_pop_idxs[pop_idx+1]]))

            grid::grid.polygon(c(min(spdf_vals$yvals[pop_idx+1,]), spdf_vals$yvals[pop_idx+1,],
                                 min(spdf_vals$yvals[pop_idx+1,])),
                               c(min(spdf_vals$xvals), spdf_vals$xvals, max(spdf_vals$xvals)),
                               default.units="native",
                               gp=grid::gpar(col=NA, fill=line_cols[base_pop_idxs[pop_idx+1]],alpha=0.3))
        }
    }

    grid::popViewport()  ## pop the spdf plot viewport i.e. navigate back up to the level of the overall viewport

    ### Part 2: the LGP label and the y-axis -----------------------------------
    # Note that this is not shown by default, nor is the y-axis shown by default,
    # because this plot is most likely to be used in the manual construction of
    # a combined Geneplot-SPDF plot, in which case the y-axis of the GenePlot
    # will double as the y-axis of the rotated SPDF plot
    if (lgp_axis)
    {
        vp3 <- grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
        grid::pushViewport(vp3)

        grid::grid.text(ylab, x = grid::unit(1.5, "cm"), rot = 90)

        grid::popViewport()  ## pop the title viewport i.e. navigate back up to the level of the overall viewport
    }

    ### Part 3: the probability density label and the x-axis -------------------
    # Note that this is shown by default, because even if the plot is used as
    # part of a manually-constructed combined plot, this section will still need
    # its own separate probability density axis (see the comment about LGP axis above)
    if (prob_density_axis) {
        vp4 <- grid::viewport(layout.pos.row = 2, layout.pos.col = 1)
        grid::pushViewport(vp4)

        grid::grid.text(xlab, y = grid::unit(0.5, "cm"))

        grid::popViewport()  ## pop the spdf plot viewport i.e. navigate back up to the level of the overall viewport
    }
}
