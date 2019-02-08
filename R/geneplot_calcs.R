### This is the code containing all the functions required to run GenePlot (apart
### from the saddlepoint functions, which are lengthy and are stored in a separate
### file).

#' Run GenePlot calculations without plotting.
#'
#' Run GenePlot calculations to obtain Log-Genotype-Probabilities for all
#' individuals from populations in \code{refpopnames} or \code{includepopnames},
#' with respect to each of the populations in \code{refpopnames}.
#'
#' All individuals in \code{dat} whose population label in the 'pop' column of
#' \code{dat} matches one of the populations in \code{refpopnames} or
#' \code{includepopnames} will be included.
#'
#' NOTE that if a population is not in \code{includepopnames}/\code{refpopnames},
#' then any alleles private to that population will NOT be included in the
#' prior / posterior.  Thus the posterior for a given refpop will change slightly
#' depending on which populations are in \code{includepopnames}/\code{refpopnames}.
#'
#' Leave-one-out will be used when calculating the log-genotype probability for
#' an individual with respect to their own reference population, if specified in
#' the inputs (default is NON leave-one-out).
#' Default is NON leave-one-out but WE STRONGLY RECOMMEND USING LEAVE-ONE-OUT,
#' ESPECIALLY FOR SMALL SAMPLES (<30).
#'
#' @param dat The data, in a data frame, with two columns labelled as 'id' and
#'    'pop', and with two additional columns per locus. Missing data at any
#'    locus should be marked as '0' for each allele.
#'    The locus columns must be labelled in the format Loc1.a1, Loc1.a2,
#'    Loc2.a1, Loc2.a2, etc.
#'    Missing data must be for BOTH alleles at any locus. Missing data for
#'    ONE allele at any locus will produce an error.
#'    See \code{\link{read_genepop_format}} for details of how to import Genepop
#'    format data files into the appropriate format.
#'
#' @param refpopnames Character vector of reference population names, that must
#'    match the values in the 'pop' column of \code{dat}.
#'
#' @param locnames Character vector, names of the loci, which must match the
#'      column names in the data so e.g. if dat has columns
#'      id, pop, EV1.a1, EV1.a2, EV14.a1, EV14.a2, etc.
#'      then you could use 'locnames = c("EV1","EV14") etc.
#'      The locnames do not need to be in any particular order but all of them
#'      must be in \code{dat}.
#'
#' @param includepopnames Character vector (default NULL) of population names to
#'      be included in the calculations. All individuals with 'pop' value in
#'      \code{includepopnames} will have their Log-Genotype-Probabilities
#'      calculated with respect to all populations in \code{refpopnames}. This
#'      input parameter should be used to specify the population labels of any
#'      new/additional individuals that you want to compare to the reference pops.
#'      For example, if the reference pops are Pop1 and Pop2, and you have some
#'      new individuals which you have labelled as PopNew, then use
#'      \code{includepopnames=c("PopNew")} to compare those individuals to Pop1
#'      and Pop2.
#'      You can specify the populations in any order, provided that they are all
#'      in \code{dat}.
#'      If NULL (default) then only the individuals from \code{refpopnames} will
#'      have their Log-Genotype-Probabilities calculated.
#'
#' @param prior (default="Rannala") String, either "Rannala" or "Baudouin",
#'      giving the choice of prior parameter for the Dirichlet priors for the
#'      allele frequency estimates. Both options define parameter values that
#'      depend on the number of alleles at each locus, k.
#'      "Baudouin" gives slightly more weight to
#'      rare alleles than "Rannala" does, or less weight to the data, so
#'      Baudouin may be more suitable for small reference samples, but there is
#'      no major difference between them. For more details, see McMillan and
#'      Fewster (2017), Biometrics.
#'      Additional options are "Half" or "Quarter" which specify parameters 1/2
#'      or 1/4, respectively. These options have priors whose parameters do not
#'      depend on the number of alleles at each locus, and so may be more suitable
#'      for microsatellite data with varying numbers of alleles at each locus.
#'
#' @param saddlepoint Boolean (default TRUE), indicates whether or not to use the
#'      saddlepoint method for imputing missing data/leave-one-out results.
#'      For more details, see McMillan and Fewster (2017), Biometrics.
#'
#' @param leave_one_out Boolean (default FALSE), indicates whether or not to
#'      calculate leave-one-out results for any individual from the reference
#'      pops. If TRUE, any individual from a reference population will have their
#'      Log-Genotype-Probability with respect to their own reference population
#'      after temporarily removing the individual's genotype from the sample
#'      data for that reference population. The individual's Log-Genotype-Probabilities
#'      with respect to all populations they are not a member of will be calculated
#'      as normal.
#'      We STRONGLY RECOMMEND using leave-one-out=TRUE for any small reference
#'      samples (<30).
#'
#' @param logten (default TRUE) Boolean, indicates whether to use base 10 for the
#'      logarithms, or base e (i.e. natural logarithms). logten=TRUE is default
#'      because it's easier to recalculate the original non-log numbers in your
#'      head when looking at the plots. Use FALSE for natural logarithms.
#'
#' @param min_loci (default 6) is the minimum number of loci that an individual
#'      must have (within the set of loci defined in \code{locnames}) for the
#'      individual to have its Log-Genotype-Probabilities calculated. Any
#'      individuals with fewer loci than this (i.e. with too many missing loci)
#'      will not have its Log-Genotype-Probabilities plotted.
#'      Thus if the calculations are based on 9 loci i.e. locnames is length 9,
#'      and min_loci=6, then every individual included in the results has data
#'      for at least 6 of these 9 loci.
#'
#' @param quantiles (default c(0.01,1.00)) Vector of probabilities, specifying the
#'      quantiles of the posterior distribution to be calculated.
#'      Default plots the 1\% and 100\% quantiles of the Log-Genotype-Probability
#'      distributions for each of the reference populations.
#'      For example, only 1\% of all possible genotypes that could arise from the
#'      given population will have Log-Genotype-Probabilities below the 1\% quantile,
#'      and 99\% of all possible genotypes arising from that population will have
#'      Log-Genotype-Probabilities above the 1\% quantile.
#'      The 100\% quantile is the maximum possible Log-Genotype-Probability that
#'      any genotype can have with respect to this population.
#'      Quantile values will be provided as attributes to the output object of
#'      calc_logprob (see the Value section.)
#'      If no quantiles are wanted, supply quantiles=NULL.
#'
#' @param Ndraw (default 100000) is only used if saddlepoint=FALSE. Defines the
#'      number of draws that will be taken from the distribution of
#'      log-posterior genotype probabilities for each reference population.
#'      These draws, i.e. simulated genotypes from the posterior distributions
#'      of the reference populations, are used when imputing the
#'      log-genotype-probabilities for individuals with
#'      missing data, or when calculating quantiles of the distribution. For
#'      more details, see McMillan and Fewster (2017), Biometrics.
#'
#' @return The structure of the output from \code{calc_logprob} and/or \code{geneplot}
#' is a data frame, with one row per individual.
#'
#' The first two columns are "id" and "pop", as in the input data.
#'
#' The next column (col3) is "status" which is "complete" or "impute" depending on
#' whether the individual had data for all loci, or had some loci missing.
#'
#' The next column (col4) is "nloci" which is how many loci the individual has data for
#'
#' The next columns are the final/imputed log-genotype probabilities for the
#' individual with respect to each of the reference populations.
#' They are named in the form "Pop1", "Pop2" etc. corresponding to the names in
#' the refpopnames input.
#'
#' Then the final columns are the "raw" log-genotype probabilities for the same pops.
#' These are named in the form "Pop1.raw", Pop2.raw", etc. again corresponding
#' to the names in refpopnames.
#'
#' For individuals with full data at all loci, i.e. no missing data, these two
#' sets of columns will be the same, and give the individual's log-genotype
#' probabilities with respect to each of the reference populations.
#'
#' For individuals with missing data at some loci then the raw values are the
#' log-genotype probabilities calculated based on the loci that *are* present in
#' the data, and the final/imputed columns, at the start of the results data frame,
#' are the imputed log-genotype probabilities for the full set of loci i.e.
#' the final LGPs for the missing-data individuals are comparable to the final LGPs
#' for the complete-data individuals.
#'
#' ---- Additional attributes of the results object ----------------------------
#'
#' At the end of \code{calc_logprob} the details of the algorithm used to calculate
#' the results are attached as attributes to the results object.
#' If your call to \code{calc_logprob} or \code{geneplot} is e.g.
#'
#'      \code{Pop1_vs_Pop2_results <- calc.logprob.func(dat, c("Pop1","Pop2"), locnames=whaleLocnames)}
#'
#'  then you would find out the attributes using \code{attributes(Pop1_vs_Pop2_results)$saddlepoint} etc.
#'
#' Other attributes attached to the results object are:
#'      \code{attributes(results)$min.loci} -- the minimum number of loci to require
#'          for any individual to be assigned, so any individual with fewer loci
#'          will be excluded from analysis
#'
#'      \code{attributes(results)$n.too.few} -- the number of individuals that have
#'          been excluded from the analysis because they had too few loci
#'
#'      \code{attributes(results)$percent.missing} -- the percentage of individuals
#'          that have been excluded, out of all those in the samples listed in
#'          allpopnames
#'
#'      \code{attributes(results)$qmat} -- the values of the plotted quantiles
#'          for the populations, with the \% labels of the quantiles as the column names
#'          e.g. if quantiles=c(0.05,0.99) was the input to chart.func then
#'          qmat will be of the form
#'
#'          \tabular{lrr}{
#'          \tab 5\% \tab 99\% \cr
#'          Pop1 \tab xx \tab xx\cr
#'          Pop2 \tab xx \tab xx}
#'
#'      \code{attributes(results)$allele_freqs} -- the posterior estimates of the allele
#'          frequencies for the populations, as a list, where each element of the
#'          list corresponds to one locus (and the list elements are named with
#'          the loci names), and at a single locus the allele
#'          frequencies are given as a matrix with the allele type names as the
#'          columns and the reference populations as the rows
#'          e.g. one locus example
#'
#'          \code{$TR3G2}
#'          \tabular{lrrrrrr}{
#'              \tab        150 \tab 158   \tab 168    \tab 172    \tab 176    \tab 180 \cr
#'              Pop1 \tab 0.125 \tab 0.125 \tab 12.125 \tab 26.125 \tab 22.125 \tab 12.125 \cr
#'              Pop2 \tab 1.125 \tab 2.125 \tab 13.125 \tab 29.125 \tab 21.125 \tab 10.125}
#'
#'          These are allele COUNT estimates, NOT PROPORTION estimates, so
#'          they do not need to add up to 1.
#'
#'      \code{attributes(results)$allpopnames} -- a vector of refpopnames, followed by include.pops names
#'          i.e. allpopnames <- c(refpopnames, include.pops)
#'
#'      \code{attributes(results)$refpopnames} -- vector of reference population names
#'
#'      \code{attributes(results)$include.pops} -- vector of included pop names for assignment
#'
#'      \code{attributes(results)$saddlepoint} -- TRUE/FALSE for whether saddlepoint was used
#'
#'      \code{attributes(results)$leave.one.out} -- TRUE/FALSE for whether leave.one.out was used
#'
#'      \code{attributes(results)$logten} -- TRUE/FALSE for whether log_10 was used (TRUE) or
#'          log_e was used (FALSE)
#'
#'      \code{attributes(results)$prior} -- "Rannala"/"Baudouin", for whether Rannala
#'          and Mountain or Baudouin and Lebrun prior was used (see McMillan & Fewster, 2017 Biometrics)
#'
#' @references McMillan, L. and Fewster, R. "Visualizations for genetic assignment
#'  analyses using the saddlepoint approximation method" (2017) \emph{Biometrics}.
#'
#'  Rannala, B., and Mountain, J. L. (1997). Detecting immigration by using multilocus
#'  genotypes. \emph{Proceedings of the National Academy of Sciences} \strong{94}, 9197--9201.
#'  Piry, S., Alapetite, A., Cornuet, J.-M., Paetkau, D., Baudouin, L., and
#'
#'  Estoup, A. (2004). GENECLASS2: A software for genetic assignment and
#'  first-generation migrant detection. \emph{Journal of Heredity} \strong{95}, 536--539.
#'
#' @author Log-Genotype-Probability calculations based on the method of Rannala
#' and Mountain (1997) as implemented in GeneClass2, updated to allow for individuals
#' with missing data and to enable accurate calculations of quantiles of the
#' Log-Genotype-Probability distributions of the reference populations.
#' See McMillan and Fewster (2017) for details.
#' Coded by Rachel Fewster and Louise McMillan.
#'
#' @importFrom stats dbeta dnorm pnorm rgamma rbeta rbinom rmultinom ecdf prcomp quantile uniroot
#'
#' @export
calc_logprob <- function(dat, refpopnames, locnames,
                         includepopnames=NULL,
                         prior="Rannala",
                         saddlepoint=T,
                         leave_one_out=F,
                         logten=T,
                         min_loci=6,
                         quantiles=c(0.01, 1.00),
                         Ndraw=100000)
{
    ## Validate inputs ---------------------------------------------------------
    validate_geneplot_inputs(dat,refpopnames,includepopnames,locnames,
                             prior,saddlepoint,leave_one_out,logten,min_loci,quantiles)

    ## Check quantiles are valid
    if (!is.null(quantiles) & (!is.vector(quantiles) || !is.numeric(quantiles) || any(quantiles > 1) || any(quantiles < 0))) stop("quantiles must be a vector of values between 0 and 1.")

    ## Make sure that dat$pop is stored as a character vector ------------------
    dat$pop <- as.character(dat$pop)
    allpopnames <- unique(c(refpopnames, includepopnames))

    ## Reduce the data frame "dat" to include only the populations in includepopnames:
    dat <- dat[!is.na(match(dat$pop, allpopnames)),]

    ## If you want to reorder refpopnames to the same order as includepopnames, use the line below; by default,
    ## keep the order specified in refpopnames because this sets the axes:
    ##      refpopnames <- includepopnames[sort(match(refpopnames, includepopnames))]

    ## Make popx2 --------------------------------------------------------------
    ## popx2 is useful so we can concatenate the vector for the first allele followed by
    ## the vector for the second allele, keeping the population information intact:
    popx2 <- rep(dat$pop, 2)

    ## Calculate the posterior allele frequencies for all populations at each locus
    ## This function returns an array with rows matching the refpopnames, and with entries being the
    ## vector (nu_1, ..., nu_k) for this row. For example the output could be (for a locus with 8 alleles):
    ##           150 154  156  164   172  178  182   184
    ## Pop1      1    1   19     6    12    1   18    10
    ## Pop2      1    1   16     1    42    1    1     1
    ## meaning that the posterior for allele frequencies for population Pop2 is (p1, ..., p8)~Dirichlet(1, 1, 16, ..., 1).
    posterior_nu_list <- lapply(locnames, calc_posterior_locus, dat, popx2, prior, refpopnames)
    names(posterior_nu_list) <- locnames

    ## Create for each locus a vector of nusums, such that the first element is the value nusum for the
    ## first refpop, the second for the second refpop, etc:
    posterior_nusumvec_list <- lapply(posterior_nu_list, rowSums)

    ## Split into complete and missing data individuals -------------------------------
    ## Now find the assignment probabilities, i.e. the log posterior genotype
    ## probabilities for each individual in dat.
    ## Split the data into complete-data individuals, dat_complete, and incomplete-data individuals, dat_missing.
    ## Just check all the locname.a1 columns for 0's, because the locname.a2 column will be the same:
    any.missing <- apply(dat[paste(locnames, ".a1", sep="")], 1, function(x) any(x==0))
    dat_complete <- dat[!any.missing,]
    dat_missing <- dat[any.missing,]
    n_complete <- nrow(dat_complete)
    n_missing <- nrow(dat_missing)

    ## Generate LGP distribution -----------------------------------------------
    ## In three situations, we need to simulate from the distribution of log-posterior genotype probabilities
    ## for each reference population.  These situations are (a) when there are any missing data, in which case
    ## the simulations enable us to place the missing-data individuals on the plot so they are represented with the same
    ## p-values as they would be if only their complete-data loci were shown; and (b) if the vector "quantiles"
    ## is supplied so that crosslines representing the corresponding quantiles are wanted on the twopop or
    ## 3spin plot; and (c) if doing leave-one-out, where we need to impute the log-posterior genotype probabilities
    ## calculated after leaving the individual out on the same scale as the full-population LGPs.

    ## Initially set up both the all_loci_SCDF_qsearch_params and all_loci_sim_logprob objects
    ## so that they can be passed in to functions (otherwise we end up with multiple versions of
    ## function calls, one for each option)
    all_loci_SCDF_qsearch_params = NULL
    all_loci_sim_logprob = NULL

    if((n_missing > 0 | !is.null(quantiles) | leave_one_out) & saddlepoint){

            # Calculate the SCDF for all loci
            all_loci_SCDF_qsearch_params <- calc.qsearch.params(posterior_nu_list, refpopnames, logten=logten,
                                                                leave_one_out=leave_one_out)
    } else if (!saddlepoint) {
        ## If saddlepoint has not been used, will need to have access to
        ## all_loci_sim_logprob as an attribute of logprob_results, in order to plot
        ## the bars in the multi-pop bar plot.
        ## So even if is.null(quantiles) and n_missing==0 and leave_one_out=FALSE
        ## we still need to calculate all_loci_sim_logprob when saddlepoint=FALSE.

        ## The logs will be base 10 if logten=T, and natural logs otherwise.
        sim.logprob.bylocus <- lapply(posterior_nu_list, sim_logprob_locus, Nsim=Ndraw, logten=logten)

        ## For both situations (n_missing>0 or !is.null(quantiles))), we need the all-loci aggregation of these
        ## probabilities: i.e. the distribution of all-locus genotype probabilities under the posterior distributions
        ## for each reference population.
        all_loci_sim_logprob <- aggregate_logprob(sim.logprob.bylocus, which.loci.names=locnames)
    }

    ## Quantiles ---------------------------------------------------------------
    ## Create the matrix of quantile information if quantiles is non-NULL, for plotting via plot_logprob_twopop or threespin.plot:
    ## quantile.mat below has the following format:
    ##              1\%   99\%
    ## Pop1   -37.3 -25.7
    ## Pop2 -29.8 -11.9
    ## meaning that an individual with a log-genotype probability of less than -37.3 is below the 1% percentile for
    ## individuals genuinely simulated from the posterior genotype distribution for population Pop1.

    if(!is.null(quantiles)){
        ## If only one quantile is required, the quantile matrix loses dimensionality and loses the nice names.
        ## Cheap hack to stop this happening is to ask for the same quantile twice:
        if(length(quantiles)==1) quantiles.vec <- rep(quantiles, 2)
        else quantiles.vec <- quantiles
        ## Make sure quantiles are in order:
        quantiles.vec <- sort(quantiles.vec)

        ###### Here, change to using quantiles from all-loci SCDF
        if (saddlepoint)
        {
            if (leave_one_out) qsearch.params <- calc.qsearch.params(posterior_nu_list, refpopnames, logten=logten,
                                                                     leave_one_out=leave_one_out)
            else qsearch.params <- all_loci_SCDF_qsearch_params
            ## Make sure to use leave-one-out for calc.quantiles.uniroot as well as
            ## for calc.qsearch.params, otherwise the K functions will not work
            ## correctly -- need to be using the twoPops/leave-one-out mode for
            ## all K functions when calculating leave-one-out quantiles
            quantile.mat = calc.quantiles.uniroot(qsearch.params,
                                                  posterior_nu_list, refpopnames,
                                                  quantiles.vec, logten=logten,
                                                  leave_one_out = leave_one_out)

            ## Add column headings to the quantile table
            colnames(quantile.mat) <- paste0(round(100*quantiles.vec,0),"%")
        }
        else
        {
            quantile.mat <- t(apply(all_loci_sim_logprob, 1, quantile, probs=quantiles.vec))
        }

        if(length(quantiles)==1) quantile.mat <- quantile.mat[,1, drop=F]
    }
    else quantile.mat <- NULL

    ## Deal with complete-data individuals --------------------------------------------

    if(n_complete > 0){
        dat_complete$status <- rep("complete", n_complete)
        dat_complete$nloc <- rep(length(locnames), n_complete)

        logprob.complete <- t(sapply(1:n_complete, calc_logprob_complete, dat_complete, refpopnames, locnames,
                                     posterior_nu_list, posterior_nusumvec_list,
                                     leave_one_out, saddlepoint, logten, Ndraw,
                                     all_loci_sim_logprob, all_loci_SCDF_qsearch_params))

        ## Put logprob together with whatever columns out of id, pop, subpop, and nloc are present,
        ## into logprob_results, a data frame:
        logprob_results_complete <- data.frame(dat_complete[,!is.na(match(names(dat_complete),
                                                                      c("id", "pop", "subpop", "status", "nloc")))],
                                               logprob.complete)
        names.complete <- names(logprob_results_complete)
    }
    else{
        logprob_results_complete <- NULL
        if (!is.na(match("subpop",names(dat)))) basicCols <- c("id", "pop", "subpop", "status", "nloc")
        else basicCols <- c("id", "pop", "status", "nloc")
        names.complete <- c(basicCols, refpopnames, paste(refpopnames, "raw", sep="."))
    }

    ## Deal with missing-data individuals ---------------------------------------------

    if(n_missing > 0){
        dat_missing$status <- rep("impute", n_missing)
        logprob.missing <- t(sapply(1:n_missing, calc_logprob_missing, dat_missing, refpopnames, locnames,
                                    posterior_nu_list, posterior_nusumvec_list,
                                    leave_one_out, saddlepoint, logten, Ndraw, sim.logprob.bylocus,
                                    all_loci_sim_logprob, all_loci_SCDF_qsearch_params))
        logprob_results_missing <- data.frame(dat_missing[,!is.na(match(names(dat_missing),
                                                                    c("id", "pop", "subpop", "status")))], logprob.missing)
        ## Reorder the columns of logprob_results_missing if necessary, in order to match those of logprob_results_complete:
        logprob_results_missing <- logprob_results_missing[ , match(names.complete, names(logprob_results_missing))]
        ## If they still don't match, I don't know why: it would need investigating.
        if(any(names(logprob_results_missing)!=names(logprob_results_complete))) stop("Naming of the complete and missing arrays is different. They shouldn't be.")

        ## Individuals with too few loci ----------------------------------------------
        ## Remove any individuals from logprob_results_missing that don't meet min_loci
        n.too.few <- length((1:n_missing)[logprob_results_missing$nloc < min_loci])
        if(n.too.few > 0){
            cat("Removing ", n.too.few, "individuals with number of loci < ", min_loci, "\n\n")
            logprob_results_missing <- logprob_results_missing[logprob_results_missing$nloc >= min_loci,]
        }

    }
    else logprob_results_missing <- NULL

    ## Combine complete-data results and missing-data results ------------------
    ## The two arrays should have been set up so that they have the same names and the function should
    ## already have dumped if they don't.
    logprob_results <- rbind(logprob_results_complete, logprob_results_missing)

    ## Add attributes to the results object ------------------------------------
    attributes(logprob_results)$min_loci <- min_loci
    if (n_missing > 0) attributes(logprob_results)$n.too.few <- n.too.few
    else attributes(logprob_results)$n.too.few <- 0

    attributes(logprob_results)$percent_missing <- count_missing_loci(logprob_results, dat)

    attributes(logprob_results)$qmat <- quantile.mat
    attributes(logprob_results)$allele_freqs <- posterior_nu_list

    attributes(logprob_results)$refpopnames <- refpopnames
    attributes(logprob_results)$includepopnames <- includepopnames

    attributes(logprob_results)$saddlepoint <- saddlepoint
    attributes(logprob_results)$leave_one_out <- leave_one_out
    attributes(logprob_results)$logten <- logten
    attributes(logprob_results)$prior <- prior

    attributes(logprob_results)$allpopnames <- allpopnames

    ## Attach the full-loci simulated distribution, which will be NULL except when
    ## saddlepoint=FALSE and quantiles are not NULL
    attributes(logprob_results)$all_loci_sim_logprob <- all_loci_sim_logprob

    logprob_results
}

validate_geneplot_inputs <- function(dat,refpopnames,includepopnames,locnames,
                                     prior,saddlepoint,leave_one_out,logten,min_loci,quantiles)
{
    if (is.null(dat)) stop("dat cannot be null.")
    if (is.null(refpopnames)) stop("refpopnames cannot be null.")
    if (is.null(locnames)) stop("locnames cannot be null.")

    if (!is.data.frame(dat)) stop("dat must be a data frame.")
    if (!("pop" %in% names(dat))) stop("dat does not have a column named 'pop' containing the population/sample labels.")
    if (is.factor(dat$pop)) stop("The 'pop' column in dat is a factor. Please convert to character and rerun.")
    # if (!is.character(dat$pop)) stop("The 'pop' column in dat must be a character vector of the population/sample labels.")

    ## Check refpopnames and includepopnames are valid
    if (!is.vector(refpopnames) || !is.character(refpopnames)) stop("refpopnames must be a character vector of population/sample names.")
    if (length(refpopnames) == 0) stop("refpopnames is an empty vector.")
    if (length(refpopnames) < 2) stop("You must specify at least two populations in refpopnames.")
    if (!all(refpopnames %in% unique(dat$pop))) stop("At least one refpopname is not found in the dataset. Please check the inputs.")

    ## Check if refpopnames have spaces, because spaces get filled in as "." in
    ## calc_logprob but then that doesn't match with plot_logprob still using
    ## the original refpopnames with spaces
    if (any(grepl(" ",refpopnames))) stop("refpopnames cannot have spaces in them. Please edit refpopnames and the pop field in the dataset before rerunning.")

    ## Check that all the locnames valid and are in the data columns
    if (!is.vector(locnames) || !is.character(locnames)) stop("locnames must be a character vector of locus names.")
    if (length(locnames) == 0) stop("locnames is an empty vector.")
    datColNames <- as.vector(outer(c(".a1",".a2"),locnames, function(allele,loc) paste0(loc,allele)))
    if (!all(datColNames %in% colnames(dat))) stop("Some of the locnames are not present in the dataset. Please check that dat has two columns for each locus, named <<locusname>>.a1 and <<locusname>>.a2.")

    if (!is.null(includepopnames))
    {
        if (!is.vector(includepopnames) || !is.character(includepopnames)) stop("includepopnames must be a character vector of sample names to be compared to the reference populations.")
        if (length(includepopnames) == 0) stop("includepopnames is an empty vector.")
        if (!all(includepopnames %in% unique(dat$pop))) stop("At least one of includepopnames is not found in the dataset. Please check the inputs.")
    }

    ## Check if any individuals are missing all their data
    allMissing <- sapply(dat$id,function(id){
        all(dat[dat$id == id,datColNames] == 0)
    })
    if (any(allMissing)) stop("At least one individual is missing data at all loci. Please remove these individuals before rerunning analysis.")

    ## Check optional inputs
    if (is.null(prior)) stop("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'.")
    if (!is.character(prior) || !is.vector(prior) || length(prior) != 1 || !(prior %in% c("Rannala","Baudouin","Half","Quarter"))) stop("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'.")

    if (is.null(saddlepoint)) stop("'saddlepoint' cannot be null. Please set to TRUE or FALSE.")
    if (!is.logical(saddlepoint) || !is.vector(saddlepoint) || length(saddlepoint)!=1) stop("'saddlepoint' must be TRUE or FALSE.")

    if (is.null(leave_one_out)) stop("'leave_one_out' cannot be null. Please set to TRUE or FALSE.")
    if (!is.logical(leave_one_out) || !is.vector(leave_one_out) || length(leave_one_out)!=1) stop("'leave_one_out' must be TRUE or FALSE.")

    if (is.null(logten)) stop("'logten' cannot be null. Please set to TRUE or FALSE.")
    if (!is.logical(logten) || !is.vector(logten) || length(logten)!=1) stop("'logten' must be TRUE or FALSE.")

    if (is.null(min_loci)) stop("'min_loci' cannot be null. Please specify a positive integer, smaller than the length of locnames.")
    if (!is.numeric(min_loci) || !is.vector(min_loci) ||
        length(min_loci)!=1 || min_loci%%1 != 0 ||
        min_loci < 1 || min_loci >= length(locnames)) stop("'min_loci' must be a positive integer, smaller than the length of locnames.")

    if (!is.null(quantiles) & (!is.vector(quantiles) || !is.numeric(quantiles) ||
                               !all(quantiles >= 0) || !all(quantiles <= 1))) stop("'quantiles' must be null, or a vector of values between 0 and 1.")
}

## calc_posterior_locus takes a single locus, and returns the parameters of the
## Dirichlet posterior for allele frequencies for each refpop at this locus.
## Let (p1, ..., pk) be the allele frequencies for a population at one locus with
## k possible alleles.
## The prior is (p1, ..., pk) ~ Dirichlet(alpha/k, ..., alpha/k) where
## alpha = 1 if prior="Rannala" or
## alpha = k if prior = "Baudouin" or
## alpha = k/2 if prior="Half" or
## alpha = k/4 if prior="Quarter".
## The posterior is (p1, ..., pk) ~ Dirichlet(nu_1, ..., nu_k) where nu_1, ..., nu_k
## are determined by the data.
## This function returns an array with rows matching the refpopnames, and with
## entries being the vector (nu_1, ..., nu_k) for this row. For example the
## output could be (for a locus with 8 alleles):
##         150  154 156  164 172  178 182  184
## Pop1    1    1   19   6   12   1   18   10
## Pop2    1    1   16   1   42   1   1    1
## meaning that the posterior for allele frequencies for population Pop2 is
## (p1, ..., p8)~Dirichlet(1, 1, 16, ..., 1).
calc_posterior_locus <- function(loc, dat, popx2, prior, refpopnames){

    ## First find all alleles present in dat (where dat includes only the
    ## reference pops specified in the call to calc_logprob) at this locus:
    a1.name <- paste0(loc, ".a1")
    a2.name <- paste0(loc, ".a2")
    dat_alleles <- c(dat[, a1.name], dat[, a2.name])

    ## allele.tab is a table of all possible alleles in the specified reference
    ## populations, and their frequency in those populations, e.g.:
    ##         150  154 156    164  172   178 182    184
    ## Pop1    0    0   18     5    11    0   17     9
    ## Pop2    0    0   15     0    41    0    0     0

    ## Using dnn=NULL below prevents extra names "popx2" and "dat_alleles" being
    ## inherited by other objects later in the coding.
    allele.tab <- table(popx2, dat_alleles, dnn=NULL)
    ## Remove missing data (code = 0):
    allele.tab <- allele.tab[ , colnames(allele.tab)!="0", drop=F]

    ## Determine choice for prior: the prior distribution of allele frequencies
    ## at each locus will be
    ## prior(p1, ..., pk) ~ Dirichlet(alpha/k, ..., alpha/k)
    ## for a locus with k possible alleles, and user chooses the value of alpha.
    k <- ncol(allele.tab)
    if(prior=="Rannala") alpha <- 1
    else if(prior=="Baudouin") alpha <- k

    else if(prior=="Half") alpha <- k/2
    else if(prior=="Quarter") alpha <- k/4

    ## Find the posterior nuvec for each population in refpopnames.
    ## It's as simple as adding alpha/k to the existing allele.tab, because the
    ## posterior is
    ## (p1, ..., pk) ~ Dirichlet ( alpha/k + n_1, ...., alpha/k + n_k ),
    ## and the values n_1, ..., n_k are exactly what are found in allele.tab.
    posterior_nu_tab <- alpha/k + allele.tab[refpopnames, , drop=F]

    ## Using the example above, if prior="Baudouin"
    ## (so prior ~ Dirichlet(1, 1, ..., 1))
    ## then posterior_nu_tab would just have 1 added to every table entry:
    ## posterior ~ Dirichlet(nuvec)
    ##         150  154 156    164  172   178 182   184
    ## Pop1    1    1   19     6    12    1   18    10
    ## Pop2    1    1   16     1    42    1    1     1
    ## To extract the posterior nuvec for a particular population, e.g. such
    ## that the posterior for Pop2 is Dirichlet(nuvec), just use posterior_nu_tab["Pop2",].

    posterior_nu_tab
}

## Functions to calculate the log genotype probability for each an individual
## with respect to each population.
## The posterior genotype probability is
## prod_loci { int_pvec   P(locus genotype | pvec)  f_posterior(pvec)  d pvec }
## and the locus-wise integral is given by the Dirichlet compound multinomial
## distribution:
## suppose the posterior ~ Dirichlet(nuvec), and nusum = sum(nuvec).
## Then:
## (a) if the individual has two distinct alleles s and t, we get:
## posterior genotype prob for this locus
##   = int_pvec P( genotype
##      = (s,t) | pvec) dDirichlet(pvec | nuvec) dpvec
##      = 2 * nu[s] * nu[t] / ( nusum * (nusum+1) )
## (b) if the individual has one repeated allele, t, we get:
## posterior genotype prob for this locus
##   = int_pvec P( genotype
##      = (t,t) | pvec) dDirichlet(pvec | nuvec) dpvec
##      = nu[t] * (nu[t] + 1) / ( nusum * (nusum+1) )
calc_logprob_complete <- function(row, dat_complete, refpopnames, locnames,
                                  posterior_nu_list, posterior_nusumvec_list,
                                  leave_one_out, saddlepoint, logten, Ndraw,
                                  all_loci_sim_logprob = NULL,
                                  all_loci_SCDF_qsearch_params = NULL){

    ## calc_logprob_complete takes a row of the data frame "dat_complete",
    ## and finds the posterior genotype probability for that
    ## individual at every locus, then combines by summing the logs.
    ## log_indiv_probvec is the output log-posterior genotype probability for this individual.
    indiv_dat <- dat_complete[row,]

    ## For leave-one-out, create a separate copy of posterior_nu_list to store the left-one-out
    ## allele frequencies for this individual
    if (leave_one_out && indiv_dat$pop %in% refpopnames) LOO.posterior_nu_list <- posterior_nu_list

    indiv_allele_freqs <- vector()
    indiv_allele_freqs.LOO <- vector()

    log_indiv_probvec <- 0
    for(loc in locnames){

        loc.table <- posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")
        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        if (length(col.a1) == 0 || length(col.a2) == 0) browser()

        ## For leave-one-out, correct the posterior allele frequencies by subtracting 1
        ## for each of the individual's alleles -- but only for the individual's own labelled population
        ## Do not do this for individuals in includepopnames
        if (leave_one_out && indiv_dat$pop %in% refpopnames)
        {
            loc.table[indiv_dat$pop, indiv_a1] = loc.table[indiv_dat$pop, indiv_a1] - 1
            loc.table[indiv_dat$pop, indiv_a2] = loc.table[indiv_dat$pop, indiv_a2] - 1
            nusumvec[indiv_dat$pop] = nusumvec[indiv_dat$pop] - 2

            ## Also change the estimated allele frequencies at this locus by removing
            ## the rat's own alleles
            ## Must subtract 1 from LOO.posterior_nu_list so that for homozygotes,
            ## take same allele away twice (if subtracting each time directly from
            ## from posterior_nu_list, end up only taking allele away once)
            LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a1] = LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a1] - 1
            LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a2] = LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a2] - 1

            indiv_allele_freqs.LOO <- c(indiv_allele_freqs.LOO,loc.table[,col.a1],loc.table[,col.a2])
        }

        indiv_allele_freqs <- c(indiv_allele_freqs,loc.table[,col.a1]+1,loc.table[,col.a2]+1)

        if(col.a1==col.a2)
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        else
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote

        if(logten) log_indiv_probvec <- log_indiv_probvec + log10(prob.loc)
        else log_indiv_probvec <- log_indiv_probvec + log(prob.loc)
    } # end for-locus loop

    ## The log_indiv_probvec returned is of the following form:
    ##       Pop1       Pop2       DBM
    ##      -13.9      -22.3        -17.4

    if (leave_one_out && indiv_dat$pop %in% refpopnames) ## If doing leave-one-out, impute individual LGP for its own population
    {
        # if (saddlepoint)
        # {
        #     post.nu.pop <- lapply(LOO.posterior_nu_list, function(x) x[individual.dat$pop,]) ## Note: have to use "lapply" to preserve correct format!
        #     post.info <- calc.multi.locus.probs.func(post.nu.pop)
        #     multi.K.params <- make.K.params(post.info$dist)
        #     mean.pop <- mu(multi.K.params)
        #
        #     indiv_p <- Fhat(log_indiv_probvec[indiv_dat$pop], post.info, mean.pop, logten=logten)
        #
        #     if (is.na(indiv_p))
        #     {
        #         print("indiv_p is NA")
        #         test <- Fhat(log_indiv_probvec[,indiv_dat$pop], post.info, mean.pop, logten=logten)
        #     }
        #
        #     if (is.infinite(indiv_p))
        #     {
        #         print("indiv_p is infinite")
        #         test <- Fhat(log_indiv_probvec[,indiv_dat$pop], post.info, mean.pop, logten=logten)
        #     }
        # }
        # else
        if (!saddlepoint)
        {
            if (is.null(all_loci_sim_logprob)) stop("all_loci_sim_logprob cannot be null when using non-saddlepoint Ndraw method.")

            ## First need to know the probability in indiv_agg of getting an answer <= the one in
            ## log_indiv_present,
            ## e.g. length((1:Ndraw)[indiv_agg["Pop1",] <= log_indiv_present["Pop1"]])/Ndraw
            indiv_p <- length((1:Ndraw)[all_loci_sim_logprob[indiv_dat$pop,] <= log_indiv_probvec[indiv_dat$pop]])/Ndraw
        }

        ###### Here, use Qhat for the full SCDF
        # if (saddlepoint)
        # {
        #     if (is.null(all_loci_SCDF_qsearch_params)) stop("all_loci_SCDF_qsearch_params cannot be null when using saddlepoint.")
        #
        #     #print("Starting to calculate the quantiles for a single population")
        #     q <- calc.quantiles.uniroot(all_loci_SCDF_qsearch_params, posterior_nu_list, indiv_dat$pop, indiv_p, logten=logten)
        #     #print("Calculated the quantiles for a single population")
        # }
        # else
        if (!saddlepoint)
        {
            ## The result is the p-value for each reference pop:
            ## e.g. indiv_p = 0.501 for rpop = Pop1 means that 50.1% of individuals from the genuine Pop1 posterior
            ## have a log-genotype probability less than or equal to the one that this individual has.
            ## Now find the corresponding quantile in all_loci_sim_logprob: this will be the imputed value for
            ## this individual.
            q <- quantile(all_loci_sim_logprob[indiv_dat$pop,], probs=indiv_p)
        }

        leave_one_out.probvec <- log_indiv_probvec
        # leave_one_out.probvec[indiv_dat$pop] <- q
        if (!saddlepoint) leave_one_out.probvec[indiv_dat$pop] <- q

        ## Add the raw (non-imputed) columns to the imputed columns to return both:
        names(log_indiv_probvec) <- paste(names(log_indiv_probvec), "raw", sep=".")
        log_indiv_complete <- c(leave_one_out.probvec, log_indiv_probvec)
    }
    else
    {
        ## For complete-data individuals when not using leave-one-out, the "raw" results are the same as the logprob.complete results:
        logprob.raw <- log_indiv_probvec
        names(logprob.raw) <- paste(names(log_indiv_probvec), "raw", sep=".")
        log_indiv_complete <- c(log_indiv_probvec, logprob.raw)
    }
    log_indiv_complete
} # end calc_logprob_complete

calc_logprob_missing <- function(row, dat_missing, refpopnames, locnames,
                                 posterior_nu_list, posterior_nusumvec_list,
                                 leave_one_out, saddlepoint, logten, Ndraw,
                                 sim.logprob.bylocus=NULL,
                                 all_loci_sim_logprob = NULL,
                                 all_loci_SCDF_qsearch_params = NULL){
    ## calc_logprob_missing takes a row of the data frame "dat_missing".
    ## It finds the posterior genotype probability for that individual at the
    ## loci that are present, and combines them by summing the logs.
    ## It then finds the p-value of that individual's genotype probability among
    ## those loci that are present, such that if the p-value is indiv_p for
    ## population A, then indiv_p is the probability that an individual generated
    ## from the known posterior genotype distribution for A has genotype
    ## probability less than or equal to that of this target individual.
    ## Once we know indiv_p, we find the quantile in the all-locus genotype
    ## probability distribution corresponding to indiv_p.  This quantile is the
    ## logprob output for this individual for refpop A.
    ## A similar process is applied to all the other reference populations.
    ##
    ## log_indiv_probvec is the output "imputed" log-posterior genotype probability
    ## for this individual.
    indiv_dat <- dat_missing[row,]

    ## For leave-one-out, create a separate copy of posterior_nu_list to store
    ## the left-one-out allele frequencies for this individual
    if (leave_one_out && indiv_dat$pop %in% refpopnames) LOO.posterior_nu_list <- posterior_nu_list

    ## Find which loci are present for this individual:
    indiv_loc.a1 <- indiv_dat[paste(locnames, ".a1", sep="")]
    locnames.missing <- locnames [which( indiv_loc.a1 ==0)]
    locnames.present <- locnames [which( indiv_loc.a1 !=0)]

    if (length(locnames.present) == 0) browser()

    ## Find the aggregate log genotype probabilities for the loci present, so we
    ## can look up this individual's values:
    if (!saddlepoint)
    {
        if (is.null(sim.logprob.bylocus)) stop("sim.logprob.bylocus cannot be null when using non-saddlepoint Ndraw method.")
        indiv_agg <- aggregate_logprob(sim.logprob.bylocus, which.loci.names=locnames.present)
    }

    log_indiv_present <- 0
    for(loc in locnames.present){
        loc.table <- posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")

        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        ## For leave-one-out, correct the posterior allele frequencies by
        ## subtracting 1 for each of the individual's alleles -- but only for
        ## the individual's own labelled population
        if (leave_one_out && indiv_dat$pop %in% refpopnames)
        {
            loc.table[indiv_dat$pop, indiv_a1] = loc.table[indiv_dat$pop, indiv_a1] - 1
            loc.table[indiv_dat$pop, indiv_a2] = loc.table[indiv_dat$pop, indiv_a2] - 1
            nusumvec[indiv_dat$pop] = nusumvec[indiv_dat$pop] - 2

            ## Also change the estimated allele frequencies at this locus by removing
            ## the rat's own alleles
            ## Must subtract 1 from LOO.posterior_nu_list so that for homozygotes,
            ## take same allele away twice (if subtracting each time directly from
            ## from posterior_nu_list, end up only taking allele away once)
            LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a1] = LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a1] - 1
            LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a2] = LOO.posterior_nu_list[[loc]][indiv_dat$pop, indiv_a2] - 1
        }

        if(col.a1==col.a2)
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        else
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote

        if(logten) log_indiv_present <- log_indiv_present + log10(prob.loc)
        else log_indiv_present <- log_indiv_present + log(prob.loc)
    } # end for-locus loop for loci with data intact

    ## We now have something like this for log_indiv_present (in the case of three
    ## refpops):
    ##      Pop1       Pop2       Pop3
    ## -8.099003 -10.186118  -8.159208
    ##
    ## We also have something like this for indiv_agg if using saddlepoint=FALSE:
    ##
    ## Pop1 -7.36 -9.00 -8.38 -8.89  .... and so on for Ndraw=10000 columns.
    ## Pop2 -4.97 -3.71 -4.81 -7.47
    ## Pop3 -6.36 -9.25 -8.45 -7.48
    ##
    ## Or for saddlepoint=TRUE we have post.info, which defines the distributions
    ## of log-genotype-probabilities for the complete set of all loci, and for
    ## the reduced set of loci present in this individual.

    impute_logprob_refpop <- function(rpop){
        ## For each refpop, find the "imputed" log-genotype probability for the individual in this refpop.

        if (saddlepoint)
        {
            ###### Here, calculate the SCDF for this set of loci
            ###### For leave-one-out, if this individual is from this refpop,
            ###### use the Leave-One-Out allele frequencies for the
            ###### reduced set of loci
            if (leave_one_out & indiv_dat$pop==rpop) post.nu.pop <- lapply(LOO.posterior_nu_list[locnames.present], function(x) x[rpop,]) ## Note: must use "lapply" to preserve correct format
            else post.nu.pop <- lapply(posterior_nu_list[locnames.present], function(x) x[rpop,]) ## Note: must use "lapply" to preserve correct format
            if (length(post.nu.pop) == 0) browser()
            # post.info <- calc.multi.locus.probs.func(post.nu.pop)
            post.info <- calc.multi.locus.probs.func(post.nu.pop, leave_one_out=leave_one_out)
            multi.K.params <- make.K.params(post.info$dist)
            # mean.pop <- mu(multi.K.params)
            # indiv_p <- Fhat(log_indiv_present[rpop], post.info, mean.pop=mean.pop, logten=logten)
            mean.pop <- mu(multi.K.params, twoPops=leave_one_out)
            indiv_p <- Fhat(log_indiv_present[rpop], post.info, mean.pop=mean.pop, logten=logten, twoPops=leave_one_out)

            if (is.na(indiv_p))
            {
                print("indiv_p is NA")
                test <- Fhat(log_indiv_present[rpop], post.info, mean.pop=mean.pop, logten=logten, twoPops=leave_one_out)
            }

            if (is.infinite(indiv_p))
            {
                print("indiv_p is infinite")
                test <- Fhat(log_indiv_present[rpop], post.info, mean.pop=mean.pop, logten=logten, twoPops=leave_one_out)
            }
        }
        else
        {
            ## First need to know the probability in indiv_agg of getting an answer <= the one in
            ## log_indiv_present,
            ## e.g. length((1:Ndraw)[indiv_agg["Pop1",] <= log_indiv_present["Pop1"]])/Ndraw
            indiv_p <- length((1:Ndraw)[indiv_agg[rpop,] <= log_indiv_present[rpop]])/Ndraw
        }

        ###### Here, use Qhat for the full SCDF
        if (saddlepoint)
        {
            if (is.null(all_loci_SCDF_qsearch_params)) stop("all_loci_SCDF_qsearch_params cannot be null when using saddlepoint.")

            #print("Starting to calculate the quantiles for a single population")
            q <- calc.quantiles.uniroot(all_loci_SCDF_qsearch_params, posterior_nu_list,
                                        rpop, indiv_p, logten=logten, leave_one_out=leave_one_out)
            #print("Calculated the quantiles for a single population")
        }
        else
        {
            if (is.null(all_loci_sim_logprob)) stop("all_loci_sim_logprob cannot be null when using non-saddlepoint Ndraw method.")

            ## The result is the p-value for each reference pop:
            ## e.g. indiv_p = 0.501 for rpop = Pop1 means that 50.1% of individuals from the genuine Pop1 posterior
            ## have a log-genotype probability less than or equal to the one that this individual has.
            ## Now find the corresponding quantile in all_loci_sim_logprob: this will be the imputed value for
            ## this individual.
            q <- quantile(all_loci_sim_logprob[rpop,], probs=indiv_p)
        }

        q
    }

    log_indiv_probvec <- sapply(refpopnames, impute_logprob_refpop)

    ## For missing-data individuals, tack the number of complete loci on to the beginning so it
    ## gels with the equivalent for the complete-data individuals data frame:
    log_indiv_probvec <- c(length(locnames.present), log_indiv_probvec)
    names(log_indiv_probvec) <- c("nloc", refpopnames)

    ## The log_indiv_probvec returned is of the following form:
    ## nloc       Pop1       Pop2       DBM
    ##   8         -13.9      -22.3        -17.4

    ## Add the raw (non-imputed) columns to the imputed columns to return both:
    names(log_indiv_present) <- paste(names(log_indiv_present), "raw", sep=".")
    c(log_indiv_probvec, log_indiv_present)
} # end calc_logprob_missing


find_nloc <- function(dat, locnames){
    ## find_nloc 10/6/10
    ## Takes input data frame and adds a column nloc specifying how many loci are successful for each individual.
    ## It also ensures that dat$pop is a character vector.
    ## EXAMPLE:
    ## gbi.dat <- find_nloc(gbi.dat)
    ## saveres("dat")
    ##
    ## OR:
    ## find_nloc(gbi.dat, locnames=loc.names.fixed[-c(1, 3, 8)])$nloc
    dat$pop <- as.character(dat$pop)
    nrats <- nrow(dat)
    ## nloc is a record for each individual giving the number of successful loci:
    nloc <- rep(0, nrats)
    loc.col.names <- paste(rep(locnames, each=2), c("a1", "a2"), sep=".")
    for(loc.col in loc.col.names){
        dat[,loc.col][is.na(dat[,loc.col])] <- 0
        ## Add 0.5 to the number of present loci for each individual which DOESN'T have missing
        ## record at this loc.allele.  The 0.5 should be replicated by the other allele for the same locus.
        nloc[dat[,loc.col]!=0] <- nloc[dat[,loc.col]!=0] + 0.5
    }
    dat$nloc <- nloc
    dat
}

count_missing_loci <- function(logprob_results, dat){
    nloci <- (ncol(dat) - 2)/2

    percent_missing <- 100*(1-sum(logprob_results$nloc)/(nloci*nrow(logprob_results)))
}

sim_logprob_locus <- function(nutab, Nsim, logten){
    ## sim_logprob_locus 19/5/11
    ## Returns the log-probabilities of animals genuinely simulated from the
    ## posterior distribution of allele frequencies specified in nutab.
    ## Operates in log_10 if logten=T, otherwise in log_e.
    ## nutab is a table of nu-vectors which are the posterior parameters for
    ## the different populations, one row per population, e.g.
    ##        150  154 156    164  172   178  182   184
    ## Pop1   1    1   19     6    12    1    18    10
    ## Pop2   1    1   16     1    42    1    1     1
    ## meaning that the posterior distribution for allele frequencies is:
    ## (p1, ..., p8) ~ Dirichlet(1, 1, 19, 6, 12, 1, 18, 10) in population Pop1;
    ## (p1, ..., p8) ~ Dirichlet(1, 1, 16, 1, 42, 1, 1, 1) in population Pop2.
    ##
    ## A single individual generated from one of these populations has the Dirichlet
    ## compound multinomial distribution for its genotype probability at this
    ## locus (assuming Hardy-Weinberg equilibrium): so if nu is the vector of
    ## Dirichlet parameters for a single population, i.e. the posterior
    ## distribution parameters, and nusum=sum(nu), then the genotype probability is:
    ## 2 * nu[s] * nu[t] / (nusum * (nusum + 1))  with probability 2 * nu[s] * nu[t] / (nusum * (nusum + 1)) ;
    ## nu[t] * (nu[t]+1) / (nusum * (nusum + 1))  with probability nu[t] * (nu[t]+1) / (nusum * (nusum + 1)).
    ##
    ## (This is equivalent to saying, let X be a discrete random variable, and
    ## let f(X) be the probability function of X: it's a function of X and so it
    ## is itself a random variable.  Then  (X=x)  =>  (f(X) = f(x)) and this
    ## happens with probability f(x).  So we add f(x) to the probability of
    ## getting (f(X) = f(x)).  Here, X is the genotype, e.g. X = (s, t) or X = (t, t);
    ## and f(X) is the genotype probability which is what we're most interested in.)
    ##
    ## This function generates Nsim realisations of the genotype probabilities:
    ## i.e. it generates Nsim probabilities between 0 and 1 that would correspond
    ## to the posterior genotype probability of an individual *genuinely* drawn from
    ## this population (as estimated by the posterior distribution) at this locus.
    ##  For example, if the results are:
    ## Pop1 : (0.2, 0.1, 0.05, 0.05, 0.3, 0.2, 0.05, ...) then we learn that
    ## individuals genuinely drawn from the posterior Pop1 population would have
    ## genotype probabilities 0.05 (commonly), 0.1, 0.2, 0.3; but maybe that an
    ## individual with genotype probability 0.002 would be vanishingly unlikely
    ## from this population.  On the other hand, if there are many alleles then
    ## we might get Pop1 : (0.002, 0.001, 0.03, 0.001, 0.002, 0.007, ...) making
    ## 0.002 perfectly respectable for an individual drawn from this population.
    ##
    ## The output from this function is a matrix with rows corresponding the the
    ## populations in nutab, and Nsim columns corresponding to the Nsim draws
    ## from each population.  The results are presented as LOG genotype probabilities:
    ## e.g.  Pop1:(log(0.02),  log(0.03),  log(0.05),  log(0.05), log(0.3), log(0.2),  log(0.05), ....)
    ## e.g.  Pop2:(log(0.30),  log(0.13),  log(0.24),  log(0.15), log(0.1), log(0.13), log(0.1),  ....)

    sim_logprob_onepop <- function(nupop){
        ## Do the simulations for a single population.
        ## To avoid problems with sample() in the case that there is only one
        ## allele, return rep(0, Nsim) if this happens, because all animals
        ## sampled from this population must have genotype probability 1 if
        ## there is only one allele, so the log genotype probability is 0.
        if(length(nupop)==1) return(rep(0, Nsim))

        ## If there is more than one allele, continue.
        ## First create the vector of all possible genotype probabilities.
        ## Because all probabilities have the same denominator (nusum * (nusum+1)),
        ## calculate the common denominator to divide by at the end:
        nusum <- sum(nupop)
        common.denom <- nusum* (nusum + 1)

        ## To find all homozygotes, the probabilities are
        ## nupop * (nupop+1) / (nusum * (nusum+1)).
        ## Create the vector of numerators of these probabilities:
        homozyg.top <- nupop * (nupop + 1)

        ## For heterozygotes, we need all pairs of alleles (call them s and t),
        ## then the probabilities are 2 nu[s] nu[t] / common.denom.
        ## Using outer() gives the values nu[s] * nu[t] for all possible pairings
        ## of s and t.  Use the lower triangle of these only, i.e. the pairs where
        ## s is not equal to t, and multiply by 2. Dscard the diagonal homozygotes.
        ## The result in heterozyg.top is a vector.
        outer.nu <- outer(nupop, nupop)
        heterozyg.top <- 2 * outer.nu[lower.tri(outer.nu, diag=F)]

        ## The vector of all genotype probabilities is geno.probvec:
        geno.probvec <- c(homozyg.top, heterozyg.top) / common.denom

        ## Now simulate Nsim log-genotype probabilities with the probabilities
        ## given in geno.probvec:
        if(logten)
            log_geno_prob_sim <- sample(log10(geno.probvec), size=Nsim, replace=T, prob=geno.probvec)
        else
            log_geno_prob_sim <- sample(log(geno.probvec), size=Nsim, replace=T, prob=geno.probvec)
        return(log_geno_prob_sim)
    }

    ## Run the simulation for each row i.e. each population in nutab
    t(apply(nutab, 1, sim_logprob_onepop))
}

aggregate_logprob <- function(sim.logprob, which.loci.names){
    ## aggregate_logprob
    ## Given a sim.logprob object corresponding to the output from
    ## lapply(posterior_nu_list, sim_logprob_locus, Nsim),
    ## this function aggregates across the loci specified by which.loci.names,
    ## and returns a matrix with rows corresponding to populations, and Nsim
    ## columns such that the Nsim columns are the SUMMED results across the
    ## specified loci of the log-genotype probabilities.
    ## It doesn't matter whether the sim.logprob is in log_10 or log_e, because
    ## it will all be in the same base.
    ##
    ## EXAMPLE: suppose which.loci.names=c("D10Rat20", "D11Mgh5"),
    ## and the input sim.logprob has Nsim=5 and looks like this:
    ##   $D10Rat20
    ##   Pop1 -0.62 -0.62 -2.39 -0.62 -0.62
    ##   Pop2 -3.49 -2.95 -3.55 -2.90 -4.80
    ##
    ##   $D11Mgh5
    ##   Pop1 -1.06 -0.823 -1.06 -4.52 -0.823
    ##   Pop2 -1.84 -3.999 -2.08 -1.84 -2.076
    ##
    ##   $D15Rat77
    ##   Pop1 -2.54 -1.49 -1.84 -6.73 -1.84
    ##   Pop2 -1.51 -6.01 -1.90 -1.51 -1.90
    ##
    ## Then the output from this function is:
    ##
    ##   Pop1 -1.680 -1.443 -3.450 -5.140 -1.443
    ##   Pop2 -5.330 -6.949 -5.630 -4.740 -6.876
    ##
    ## where the results for Pop1 are added across the two specified loci only,
    ## and likewise for Pop2:
    ## e.g. -0.62 - 1.06 = -1.680 for the first element; etc.

    subset.list <- sim.logprob[which.loci.names]
    subset.sum <- subset.list[[1]]
    n.in.list <- length(which.loci.names)

    ## If there is more than one locus, add the others to the first one:
    if(n.in.list > 1)
        for(locn in 2:n.in.list) subset.sum <- subset.sum + subset.list[[locn]]

    subset.sum
}

rDirichlet <- function(ndraw, alphvec){
    ## See Wikipedia on the Dirichlet distribution for confirmation:
    gamdraw <- rgamma(ndraw, shape=alphvec, rate=1)
    gamdraw / sum(gamdraw)
}
