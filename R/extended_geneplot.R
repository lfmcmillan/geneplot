#' Plot GenePlot with two additional side plots showing directional analyses of
#' the connectivity of the two reference populations.
#'
#' The side plots show the genetic distribution of each reference population and
#' the other population compared with it. The overlap area within each of the
#' side plots corresponds to the Overlap Area measures obtained by
#' \link{calc_spdfs}.
#'
#' Suppose populations A and B are plotted on the GenePlot with the x-axis
#' showing fit to B (i.e. Log-Genotype Probabilities with respect to B) and the
#' y-axis showing fit to A (i.e. Log-Genotype Probabilities with respect to A).
#'
#' Then the auxiliary plot at the bottom of the GenePlot will show the
#' saddlepoint distribution plots for baseline B. The solid curve shows the
#' genetic distribution of baseline B with itself: it is the distribution of
#' Log-Genotype Probabilities for all genotypes that could arise from B. The
#' dashed curve matching the colour of the other reference population, A, shows
#' the comparison distribution of A into B. For all the genotypes that could
#' arise from B, this shows how often those would arise in A. If B and A are
#' very different genetically then a lot of the genotypes that would have a very
#' good fit to B may only occur rarely in A, and the genotypes that occur
#' commonly in A may have a poor fit to B, and this would be shown by a low
#' amount of overlap between these two curves. On the other hand, if B and A are
#' genetically similar then genotypes that commonly occur in A also have a good
#' fit to B, and so there will be a high amount of overlap between the two
#' curves.
#'
#' The auxiliary plot on the left of the GenePlot will show similar curves, but
#' for baseline A, so the solid curve is the genetic distribution of A with itself
#' and the dashed curve is the genetic distribution of B into A. A low overlap
#' here means that genotypes that commonly occur in B would have a poor fit to A.
#' A high overlap here means that genotypes that commonly occur in B would have
#' a good fit to A as well as B.
#'
#' The user can also specify \code{includepopnames}, which are additional
#' populations or groups of individuals in the dataset to be plotted on the
#' central GenePlot, to show their fit to reference populations A and B.
#'
#' This function runs the GenePlot and saddlepoint distribution calculations
#' and produces the plot. If you want to split these up into separate code steps,
#' then run \link{calc_geneplot_spdfs} to perform the calculations,
#' \link{prepare_plot_params} to set up the plotting details and
#' \link{plot_geneplot_spdfs} to produce the extended plot.
#'
#'
#' @returns Output is the output of \link{calc_geneplot_spdfs}
#'
#' @seealso \link{calc_spdfs}, \link{calc_geneplot_spdfs}, \link{plot_geneplot_spdfs}, \link{prepare_plot_params}
#'
#' @references McMillan, L. F., "Concepts of statistical analysis,
#'   visualization, and communication in population genetics" (2019). Doctoral
#'   thesis, https://researchspace.auckland.ac.nz/handle/2292/47358.
#'
#' @references McMillan, L. and Fewster, R. "Visualizations for genetic assignment analyses
#'  using the saddlepoint approximation method" (2017) \emph{Biometrics}.
#'
#' @references Rannala, B., and Mountain, J. L. (1997). Detecting immigration by using
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
#' @inherit calc_geneplot_spdfs
#' @inherit plot_geneplot_spdfs
#' @inherit prepare_plot_params
#'
#' @export
extended_geneplot <- function(dat,
                              refpopnames,
                              locnames,
                              includepopnames=NULL,
                              quantiles_vec,
                              colvec=NULL, shapevec=NULL, line_widths=NULL,
                              orderpop=NULL, axispop=NULL,
                              display_names=NULL, xyrange=NULL,
                              ylim_input=NULL, mark_impute=T,
                              geneplot_multiplier=3,
                              show_overlap_areas=F, show_legend=T,
                              show_legend_below=T, legend_width=NULL,
                              show_title=T, title_text=NULL,
                              grayscale_quantiles=F,
                              show_include_ids=F,
                              prior="Rannala",
                              leave_one_out=F,
                              logten=T,
                              saddlepoint=T,
                              rel_tol=.Machine$double.eps^0.25,
                              abs_tol=rel_tol, npts=1000) {

    if (!is.vector(refpopnames) || length(refpopnames) != 2) stop("refpopnames should be a vector of length 2 giving the 2 refpops. This function only works for 2 refpops.")

    if (!is.null(axispop)) {
        if (length(axispop) != 2 || any(!(axispop %in% refpopnames))) stop("Population names in axispop do not match refpopnames.")
        if (any(!(c("x","y") %in% names(axispop)))) stop("axispop must be of the form of the form c(x='Pop1', y='Pop2')")
    }
    if (!is.null(orderpop)) {
        if (length(orderpop) != 2 || any(!(orderpop %in% refpopnames))) stop("Population names in orderpop do not match refpopnames.")
    }

    calc_results <- calc_geneplot_spdfs(dat=dat,refpopnames=refpopnames,locnames=locnames,
                                        includepopnames=includepopnames,prior=prior,logten=logten,
                                        saddlepoint=saddlepoint,leave_one_out=leave_one_out,
                                        quantiles_vec=quantiles_vec,rel_tol=rel_tol,abs_tol=abs_tol)

    plot_params <- prepare_plot_params(logprob=calc_results$logprob,
                                       logten=logten,xyrange=xyrange,
                                       colvec=colvec,shapevec=shapevec,
                                       orderpop=orderpop, axispop=axispop,
                                       display_names=display_names,
                                       ylim_input=ylim_input,
                                       geneplot_multiplier=geneplot_multiplier,
                                       mark_impute=mark_impute,
                                       grayscale_quantiles=grayscale_quantiles)

    plot_geneplot_spdfs(calc_results$logprob,calc_results$spdf_vals1,
                        calc_results$spdf_vals2,plot_params,
                        show_overlap_areas=show_overlap_areas,
                        show_quantiles=!is.null(quantiles_vec),
                        show_legend=show_legend, show_legend_below=show_legend_below,
                        legend_width=legend_width,
                        show_title=show_title, title_text=title_text,
                        show_include_ids=show_include_ids)

    calc_results
}

#' Lower-level function for preparing the plot details step of
#' \link{extended_geneplot}. Prepare the colour and shape for each population,
#' the order in which to plot the populations and the axis labels.
#'
#' Use after \link{calc_geneplot_spdfs} and before \link{plot_geneplot_spdfs}.
#'
#' @param logprob Log-Genotype Probabilities for individual genotypes with
#'   respect to each of the reference populations. Required for GenePlot,
#'   available as an output from \link{calc_geneplot_spdfs}.
#'
#' @param refpopnames Character vector of reference population names, that must
#'   match two values in the 'pop' column of \code{dat}. The SPDF methods
#'   currently only work for a pair of baseline populations, so
#'   \code{refpopnames} must be length 2.
#'
#' @param colvec (Optional) Vector of colours for the populations, starting with
#'   the reference populations in the order \code{refpopnames} and followed by
#'   any included populations in the order \code{includepopnames}. The same
#'   colours are used in the GenePlot and the side plots. Colours can be
#'   specified using rgb objects, hexadecimal codes, or any of the R colour
#'   names (see \url{http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf} for
#'   a PDF of R colours).
#'
#' @param shapevec (Optional) Vector of shapes for the populations plotted on
#'   the central GenePlot. These are named shapes from the following list:
#'   "Circle", "Square", "Diamond", "TriangleUp", "TriangleDown", "OpenSquare",
#'   "OpenCircle", "OpenTriangleUp", "Plus", "Cross", "OpenDiamond",
#'   "OpenTriangleDown", "Asterisk" which correspond to the following pch values
#'   for R plots: 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 8. Do not use the
#'   numbers, use the words, which will be automatically converted within
#'   plot_logprob into the appropriate codes.
#'
#' @param line_widths (Optional) Vector of line widths for all the populations
#'   for the genetic distribution side plots. Default widths are 4. Length must
#'   be the number of reference populations plus the number of additional
#'   included populations.
#'
#' @param cexpts (Optional) Point size for points in the central GenePlot
#'   (default is 1).
#'
#' @param orderpop (Optional) Vector of names of the reference populations and
#'   include populations, to indicate the order their points should be plotted
#'   within the GenePlot. The first will be plotted first and so will appear to
#'   be "beneath" all the other populations, and the last will be plotted last
#'   and so will appear to be "above" all the the other populations. Use if you
#'   have a particular population whose individuals you are interested in and
#'   need to see on top of the rest: put this population last in
#'   \code{orderpop}. Default is NULL, in which case populations are plotted in
#'   order of size, so the population with the largest number of points is
#'   plotted at the bottom, and the population with the smallest number of
#'   individuals/points is plotted over the top, so as not to be obscured.
#'
#' @param axispop (Optional) Vector of reference populations, indicating which
#'   to plot as the baseline on the x-axis and which on the y-axis. The default
#'   is that the first entry in \code{refpopnames} is plotted on the x-axis.
#'
#' @param display_names (default \code{refpopnames}) Use this to supply
#'   alternative display names for the populations. The refpopnames, as columns
#'   in the dataset, cannot have spaces, for example, whereas the display names
#'   can have spaces.
#'
#' @param axis_labels (Optional) Vector of axis labels. Defaults are constructed
#'   depending on whether the \code{logten} argument is TRUE or FALSE, and
#'   whether \code{short_axis_labels} is TRUE or FALSE.
#'
#' @param short_axis_labels (Optional) Default is FALSE, If TRUE, a shortened
#'   version of the default axis labels is used e.g. Log-Genotype Probability is
#'   shortened to LGP. Use this if you want to display a small plot without
#'   space for the full axis labels.
#'
#' @param xyrange (Optional) Numerical limits for the GenePlot axes, which also
#'   form the x-axes of the two side plots. The plot is always symmetrical, so
#'   the same limits will be used for the x and the y axes. This is to ensure
#'   that comparisons between the two reference populations are fair. The values
#'   are a range of Log-Genotype Probability values so the maximum is 0. Default
#'   is slightly wider than the range of the calculated Log-Genotype
#'   Probabilities for all individuals in the plot.
#'
#' @param ylim_input (Optional) Numerical limits for the y-axes of the two side
#'   plots. Has no effect on the central GenePlot. Run with the defaults first
#'   to see what the range of values is for the given populations.
#'
#' @param geneplot_multiplier (default 3) Ratio of width of GenePlot to
#'   sideplots. Default is 3, i.e. the GenePlot height will be three times the
#'   height of the bottom plot and the GenePlot width will be three times the
#'   width of the left-hand plot.
#'
#' @param mark_impute (default FALSE) Boolean, indicates whether to mark
#'   individuals with missing data using asterisks on the GenePlot.
#'
#' @param grayscale_quantiles (default FALSE) Used for plots with 2 reference
#'   pops. FALSE (default) plots the quantile lines using colvec colours TRUE
#'   plots the quantile lines in gray (as the default colours can be quite pale,
#'   the grayscale quantile lines can be easier to see than the default coloured
#'   ones).
#'
#' @param logten (default TRUE) Boolean, indicates whether to use base 10 for
#'   the logarithms, or base e (i.e. natural logarithms). logten=TRUE is default
#'   because it's easier to recalculate the original non-log numbers in your
#'   head when looking at the plots. Use FALSE for natural logarithms.
#'
#' @export
prepare_plot_params <- function(logprob, refpopnames, logten=TRUE,
                                colvec=NA, shapevec=NA,
                                line_widths=NULL, cexpts=1,
                                orderpop=NULL, axispop=NULL,
                                display_names=NULL, axis_labels=NULL,
                                short_axis_labels=FALSE, xyrange=NULL,
                                ylim_input=NULL, geneplot_multiplier=3,
                                mark_impute=TRUE, grayscale_quantiles=FALSE) {

    refpopnames <- attributes(logprob)$refpopnames
    allpopnames <- attributes(logprob)$allpopnames
    npop <- length(allpopnames)
    quantile_mat <- attributes(logprob)$qmat

    if (is.null(display_names)) display_names <- allpopnames

    ## Shapes and colours: -----------------------------------------------------
    shapevec_default <- rep(21:25, npop)[1:npop]
    if (length(shapevec)>1)
    {
        ## If some shapes are supplied then shapevec has length npop, otherwise shapevec=NA (length 1)
        shapes_ref <- c("Circle", "Square", "Diamond", "TriangleUp", "TriangleDown", "OpenSquare", "OpenCircle",
                        "OpenTriangleUp", "Plus", "Cross", "OpenDiamond", "OpenTriangleDown", "Asterisk")
        numbers_ref <- c(21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 8)
        shapes_inds <- match(shapevec, shapes_ref)
        shapevec <- numbers_ref[shapes_inds]
        shapevec[is.na(shapevec)] <- shapevec_default[is.na(shapevec)]
    } else {
        shapevec <- shapevec_default
    }

    colvec_default <- grDevices::hcl.colors(npop,"Roma")
    ## If some colours are supplied then colvec has length npop, otherwise colvec=NA (length 1)
    if (length(colvec)>1) {
        colvec[is.na(colvec)] <- colvec_default[is.na(colvec)]
    } else {
        colvec <- colvec_default
    }

    ## In case xyrange is entered the wrong way round by the user:
    if(!is.null(xyrange)) xyrange <- sort(xyrange)

    ## Populations and loci: ---------------------------------------------------

    ## Order the populations before plotting
    if(is.null(orderpop)){
        ## If no plotting order is specified (orderpop not supplied), arrange
        ## allpopnames to be in decreasing order of sample size:
        allpopsizes <- numeric(length(npop))
        for(pp in 1:npop) allpopsizes[pp] <- length(logprob$pop[logprob$pop==allpopnames[pp]])
        order_inds <- order(allpopsizes, decreasing=T)
        orderpop <- allpopnames[order_inds]
    } else {
        ## If orderpop is supplied, order_inds is gained from matching orderpop
        ## straight to allpopnames:
        order_inds <- match(orderpop, allpopnames)
    }
    allpopnames <- allpopnames[order_inds]

    ## Reorder shapevec and colvec by the same ordering to keep the correct
    ## assignment of colours and shapes to populations:
    shapevec <- shapevec[order_inds]
    colvec <- colvec[order_inds]

    ## Sort out colours for plotting points. Any populations that have pch
    ## between 21 and 25 should have rimvec set to 1 for the rim:
    rimvec <- colvec
    rimvec[21 <= shapevec & shapevec <= 25] <- 1

    ## If axispop is not specified, and there are two reference populations,
    ## make axispop the same as refpopnames:
    if(is.null(axispop) & length(refpopnames)==2){
        axispop <- refpopnames
        names(axispop) <- c("x", "y")
    }

    ## Now reorder the pops in the order of axispop:
    orig_refpopnames <- refpopnames
    refpopnames <- axispop[c("x", "y")]
    display_names <- display_names[match(orig_refpopnames,refpopnames)]

    ## Set up for plotting:
    if(is.null(xyrange)) xyrange <- range(c(logprob[, refpopnames[1]], logprob[, refpopnames[2]])) +
        c(-0.05, 0.05) * diff(range(c(logprob[, refpopnames[1]], logprob[, refpopnames[2]])))

    if (is.null(axis_labels) || !is.vector(axis_labels) || !is.character(axis_labels) || length(axis_labels) != 2) {
        if (logten && short_axis_labels) label <- "LGP10 for population"
        else if (logten && !short_axis_labels) label <- "Log10 genotype probability for population"
        else if (!logten && short_axis_labels) label <- "LGP for population"
        else if (!logten && !short_axis_labels) label <- "Log genotype probability for population"

        if (!is.null(display_names) && length(display_names) == length(refpopnames)) {
            x_lgp_label <- paste(label, display_names[1])
            y_lgp_label <- paste(label, display_names[2])
        } else {
            x_lgp_label <- paste(label, refpopnames[1])
            y_lgp_label <- paste(label, refpopnames[2])
        }
    } else {
        x_lgp_label <- axis_labels[1]
        y_lgp_label <- axis_labels[2]
    }

    split_dat <- split(logprob, logprob$pop)

    ## Prepare inputs for SPDF plots -------------------------------------------

    if (is.null(line_widths)) line_widths <- rep(4,times=2)

    ## Don't try to mark imputed LGPs if there aren't any
    if (length(which(logprob$status == "impute")) == 0) mark_impute <- FALSE

    list(allpopnames=allpopnames, logten=logten, split_dat=split_dat,
         colvec=colvec, shapevec=shapevec, rimvec=rimvec,
         orderpop=orderpop, axispop=axispop,
         line_widths=line_widths, cexpts=cexpts,
         xyrange=xyrange, ylim_input=ylim_input, display_names=display_names,
         x_lgp_label=x_lgp_label, y_lgp_label=y_lgp_label,
         geneplot_multiplier=geneplot_multiplier,
         mark_impute=mark_impute, grayscale_quantiles=grayscale_quantiles)
}

#' Lower-level function for carrying out the plotting step of
#' \link{extended_geneplot}.
#'
#' Use after running \link{calc_geneplot_spdfs} to
#' produce \code{logprob}, \code{spdf_vals1} and \code{spdf_vals2} and
#' \link{prepare_plot_params} to set up the plotting details.
#'
#' @param logprob Log-Genotype Probabilities for individual genotypes with
#' respect to each of the reference populations. Required for GenePlot,
#' available as an output from \link{calc_geneplot_spdfs}.
#'
#' @param spdf_vals1 Genetic Distributions with first reference population as
#' the baseline population. Required for the auxiliary side plots, available as
#' an output from \link{calc_geneplot_spdfs}.
#'
#' @param spdf_vals2 Genetic Distributions with second reference population as
#' the baseline population. Required for the auxiliary side plots, available as
#' an output from \link{calc_geneplot_spdfs}.
#'
#' @param plot_params Plotting options and settings, output from
#' \link{prepare_plot_params}. List includes \code{allpopnames, logten, split_dat,
#' colvec, shapevec, rimvec, orderpop, axispop, line_widths, cexpts,
#' mark_impute, xyrange, ylim_input, display_names, x_lgp_label, y_lgp_label,
#' grayscale_quantiles}.
#'
#' @param show_overlap_areas (default FALSE) If TRUE, print the numerical
#'   Overlap Area values beneath the extended plot.
#'
#' @param show_legend (default TRUE) If TRUE, plot the legend.
#'
#' @param show_legend_below (default TRUE) If TRUE, show the legend as a row
#'   of coloured shapes and labels below the extended plot. If FALSE, show the
#'   legend within the GenePlot. Default is TRUE because putting the legend
#'   inside the central GenePlot can hide points from individual genotypes.
#'
#' @param legend_width (default 4cm) Change the width of the displayed legend.
#'   Units are cm.
#'
#' @param show_title (default TRUE) Include a title for the plot.
#'
#' @param title_text (Optional) Provide alternative title text for the plot.
#'
#' @param only_plot_baseline_pop (default FALSE) Plot the baseline population in
#'   each side plot i.e. the corresponding reference population, and do not plot
#'   the comparison population.
#'
#' @param show_quantiles (default TRUE) If TRUE, plot quantiles of each baseline
#'   population in the side plots.
#'
#' @param show_include_ids (default FALSE) Ignored unless \code{includepopnames}
#'   is supplied. If TRUE, plots individual genotypes from the additional
#'   included pops as points but also displays the ID for each genotype. This
#'   can be useful if one or more genotype shows unusual patterns of fit to the
#'   reference populations and you want to identify those genotypes.
#'
#' @returns A list with the following components:
#'
#'     \code{allpopnames} Vector of populations included in the plot, starting
#'     with the reference populations followed by the included populations.
#'
#'     \code{logten} Whether the log scales in the plot are base 10 (default)
#'     or base e.
#'
#'     \code{split_dat} The Log-Genotype Probabilities for the supplied genotypes,
#'     split up into a list of data frames, with each data frame corresponding
#'     to the genotypes from a single named population and containing their
#'     Log-Genotype Probabilities with respect to both of the reference
#'     populations.
#'
#'     \code{colvec, shapevec, rimvec} Details of the colours and shapes used
#'     for the populations, and the outline colours for the GenePlot points
#'
#'     \code{orderpop, axispop} Order in which the populations were plotted in
#'     the central GenePlot, from back to front, and the x-y order of the two
#'     reference pops.
#'
#'     \code{line_widths} Line widths used in the side plot saddlepoint
#'     distribution curves.
#'
#'     \code{cexpts} Point sizes used in the central GenePlot.
#'
#'     \code{mark_impute, xyrange, ylim_input, display_names} Records the input
#'     options used to produce the plot.
#'
#'     \code{x_lgp_label, y_lgp_label} Full axis labels used in the plot.
#'
#'     \code{grayscale_quantiles} Whether the quantiles are plotted in grey or
#'     using the assigned population colours.
#'
#' @importFrom graphics polygon
#'
#' @export
plot_geneplot_spdfs <- function(logprob, spdf_vals1, spdf_vals2, plot_params,
                                show_overlap_areas=TRUE,
                                show_legend=TRUE, show_legend_below=TRUE,
                                show_title=TRUE, title_text=NULL,
                                show_quantiles=TRUE, only_plot_baseline_pop=FALSE,
                                legend_width=4,show_include_ids=FALSE)
{
    ## Take the refpops from axispop (which defaults to the original refpopnames anyway)
    refpopnames <- plot_params$axispop[c("x","y")]
    includepopnames <- attributes(logprob)$includepopnames
    quantile_mat <- attributes(logprob)$qmat

    if (spdf_vals1$refpopA != refpopnames[1] || spdf_vals2$refpopA != refpopnames[2]) {
        ## If the axispops are the other way around, need to swap the data to go with them
        if (spdf_vals1$refpopA == refpopnames[2] && spdf_vals2$refpopA == refpopnames[1]) {
            spdf_vals1_temp <- spdf_vals1
            spdf_vals1 <- spdf_vals2
            spdf_vals2 <- spdf_vals1
            rm(spdf_vals1_temp)
        }
        else stop("The two refpopnames from logprob do not match the baseline pops for spdf_vals1 and spdf_vals2.")
    }

    allpopnames <- plot_params$allpopnames
    npop <- length(allpopnames)
    logten <- plot_params$logten
    split_dat <- plot_params$split_dat
    colvec <- plot_params$colvec
    shapevec <- plot_params$shapevec
    line_widths <- plot_params$line_widths
    display_names <- plot_params$display_names
    x_lgp_label <- plot_params$x_lgp_label
    y_lgp_label <- plot_params$y_lgp_label
    xyrange <- plot_params$xyrange

    refpop_display_names <- display_names[allpopnames %in% refpopnames]
    includepop_display_names <- display_names[allpopnames %in% includepopnames]

    refpop_colvec <- colvec[allpopnames %in% refpopnames]
    includepop_colvec <- colvec[allpopnames %in% includepopnames]

    refpop_shapevec <- shapevec[allpopnames %in% refpopnames]
    includepop_shapevec <- shapevec[allpopnames %in% includepopnames]

    # Set up (4+?)x4 layout, of which first row will be the title.
    # In next two rows top left will be SPDF plot of 1st refpop, top right will
    # be GenePlot and bottom right will be SPDF plot of 2nd refpop.
    # fourth row will be x-axis label for bottom SPDF plot and GenePlot.
    # Then two optional rows for the overlap areas.
    nrefpop <- length(refpopnames)
    refpop_legend_height <- 1+ceiling(nrefpop/4)
    nincludepop <- length(includepopnames)
    show_include <- (nincludepop > 0)
    includepop_legend_height <- 1+ceiling(nincludepop/4)
    element_height_units <- c("cm","null","null","cm","cm","cm","cm","cm")
    element_heights <- c(show_title,plot_params$geneplot_multiplier,1,1.7,
                         show_legend_below, show_legend_below*show_include,
                         0.9*show_overlap_areas, 0.9*show_overlap_areas)
    element_height_units <- element_height_units[(element_heights != 0)]
    element_heights <- element_heights[(element_heights != 0)]

    lay <- grid::grid.layout(nrow=length(element_heights), ncol=4,
                             heights = grid::unit(element_heights,element_height_units),
                             widths = grid::unit(c(0.5,1,plot_params$geneplot_multiplier,2),c("cm","null","null","cm")))

    # Set up the overall viewport
    vp1 <- grid::viewport(layout = lay)

    grid::grid.newpage()
    grid::pushViewport(vp1)

    ### Part 1 : the title -----------------------------------------------------
    if (show_title){
        if (is.null(title_text)) title_text <- paste("GenePlot for ",display_names[1]," vs. ",display_names[2])
        grid::grid.text(title_text,vp=grid::viewport(layout.pos.row = 1, layout.pos.col = 2:4),
                        gp=grid::gpar(fontface="bold"))
    }

    ### Part 2 : the main GenePlot ---------------------------------------------
    vp3 <- grid::viewport(layout.pos.row = show_title+1, layout.pos.col = 3,
                          xscale = xyrange, yscale = xyrange, gp=grid::gpar(cex=1))
    grid::pushViewport(vp3)

    # Draw the y-axis on the right hand side of the plot rather than the left
    geneplotRect <- grid::grid.rect()
    grid::grid.xaxis(main=TRUE, name="xaxisGenePlot")

    # Now reduce the spacing between the tick marks and the labels
    grid::grid.force()
    grid::grid.edit("xaxisGenePlot::labels", y=grid::unit(-0.5,"cm"))

    grid::grid.yaxis(main=FALSE, name="yaxisGenePlot")

    plot_pop <- function(whichpop){
        ## Plot results for a single population.  If whichpop=1 then plot the
        ## first population in allpopnames, etc.
        pop <- allpopnames[whichpop]
        popdat <- split_dat[[pop]]
        p1 <- popdat[[refpopnames[1]]]
        p2 <- popdat[[refpopnames[2]]]
        points = grid::grid.points(x=p1, y=p2, pch=shapevec[whichpop],
                                   gp=grid::gpar(col=plot_params$rimvec[whichpop], fill=colvec[whichpop],
                                           cex=plot_params$cexpts))

        if (pop %in% includepopnames && show_include_ids==TRUE) {
            ## Have to specify units as "native" because whereas grid::grid.points has
            ## "native" as default value for default.units, grid::grid.text has "npc"
            ## as default value for default.units
            grid::grid.text(x=p1, y=p2, label=popdat$id, default.units = "native",
                            gp=grid::gpar(col="black", fontsize=12, fontface="bold"))
        }
    }
    sapply(1:npop, plot_pop)

    # Add the diagonal lines
    grid::grid.abline(0, 1, gp=grid::gpar(col="black", lwd=2))
    ## Next two lines if wanting to show the Geneclass scores = 1 lines:
    if (logten) {
        logDiff <- log10(10)
    } else {
        logDiff <- log(10)
    }
    grid::grid.function(function(z) list(x=z-logDiff,y=z), n=101,
                        gp=grid::gpar(col="darkgrey", lty=1, lwd=1))
    grid::grid.function(function(z) list(x=z+logDiff,y=z), n=101,
                        range=c(xyrange[1],xyrange[2]-logDiff),
                        gp=grid::gpar(col="darkgrey", lty=1, lwd=1))

    ## Add crosslines for quantiles if required:
    if(!is.null(quantile_mat)){
        if (plot_params$grayscale_quantiles)
        {
            col1 <- "black"
            col2 <- "black"
        }
        else
        {
            col1 <- colvec[1]
            col2 <- colvec[2]
        }

        for (idx in 1:ncol(quantile_mat))
        {
            ## ONLY show quantiles that are within the xyrange of the graph
            if (quantile_mat[refpopnames[1],idx] > xyrange[1] &&
                quantile_mat[refpopnames[1],idx] < xyrange[2]) {
                ## Note that the labels for the quantiles on the x-axis need to be further
                ## from the quantile lines than are the labels for the quantiles on the
                ## y-axis, because the labels are all displayed horizontally
                grid::grid.text(label=colnames(quantile_mat)[idx], gp=grid::gpar(col=col1, fontface="bold"),
                                x=grid::unit(quantile_mat[refpopnames[1],idx], "native"),
                                y=grid::unit(xyrange[1]+diff(xyrange)/20,"native"), rot=0, just="right")
                ## For vertical lines, have to use grid::grid.function rather than grid::grid.abline
                grid::grid.lines(x=rep(quantile_mat[refpopnames[1],idx],times=101),
                                 y=seq(xyrange[1],xyrange[2],length.out=101),
                                 default.units="native",gp=grid::gpar(col=col1, lty=2, lwd=1.5))
            }
            if (quantile_mat[refpopnames[2],idx] > xyrange[1] &&
                quantile_mat[refpopnames[2],idx] < xyrange[2]) {
                grid::grid.text(label=colnames(quantile_mat)[idx], gp=grid::gpar(col=col2, fontface="bold"),
                                y=grid::unit(quantile_mat[refpopnames[2],idx]-diff(xyrange)/40,"native"),
                                x=grid::unit(xyrange[1]+diff(xyrange)/20,"native") ,rot=0)
                grid::grid.abline(intercept=quantile_mat[refpopnames[2],idx],
                                  slope=0, gp=grid::gpar(col=col2, lty=2, lwd=1.5))
            }
        }
    }

    ## Flag the points gained by imputation if required:
    if(plot_params$mark_impute)
    {
        grid::grid.points(x=logprob[logprob$status=="impute", refpopnames[1]],
                          y=logprob[logprob$status=="impute", refpopnames[2]], pch="*")
    }

    ### Part 2a : the in-Geneplot legend ----
    if (show_legend & !show_legend_below) {
        ## Add the legend, putting it in top left to minimise chance of it hiding any
        ## of the GenePlot points
        if (is.null(legend_width)) legend_width <- 4
        ## Height of legend row in INCHES -- must be inches because when we fetch the
        ## height of the whole GenePlot using heightDetails, the output is in inches,
        ## so if we want to combine that with the legend row height, the legend row
        ## height must also be in inches
        legendRowHeight <- 0.335
        geneplotHeight <- grid::heightDetails(geneplotRect)
        vp=grid::viewport(x=grid::unit(0,"cm"),y=geneplotHeight-grid::unit(legendRowHeight/2*npop,"inches"), just="left",
                          width=grid::unit(legend_width,"cm"),height=grid::unit(legendRowHeight*npop,"inches"))
        grid::pushViewport(vp)
        grid::grid.rect()
        grid::grid.legend(labels=display_names, ncol=npop, nrow=npop, pch=shapevec, vgap=grid::unit(0.7,"lines"),
                          gp=grid::gpar(col=plot_params$rimvec, fill=colvec, cex=1))
        ## Add in a chunk of blank space to widen the left margin of the legend
        grid::grid.pack(grid::grid.grep("frame", grep=TRUE, strict=TRUE),
                        grid::rectGrob(width=grid::unit(0.5, "cm"), gp=grid::gpar(col=NA)),
                        side="left")
        ## leave the legend and geneplot viewports and navigate back up to the
        ## level of the overall viewport
        grid::upViewport(2)
    } else {
        grid::upViewport(1)
    }

    ### Part 3 : the left-hand rotated SPDF plot -------------------------------
    ### This plot has refpop 2 as the baseline pop

    base_pop_idxs <- spdf_vals2$base_pop_idxs

    if (is.null(plot_params$ylim_input))
    {
        yvals_vec <- as.vector(spdf_vals2$yvals)
        ylim_actual2 <- range(yvals_vec[which(!is.nan(yvals_vec) & !is.infinite(yvals_vec))])
        ylim_actual2[1] <- 0
    } else ylim_actual2 <- plot_params$ylim_input

    xlim_actual <- xyrange

    vp4 <- grid::viewport(layout.pos.row = show_title+1, layout.pos.col = 2,
                          xscale = c(ylim_actual2[2],ylim_actual2[1]), yscale = xlim_actual, gp=grid::gpar(cex=1))
    grid::pushViewport(vp4)

    grid::grid.rect()

    ## Note that what would normally be the x-axis (ie. the LGP axis) is being
    ## drawn as the y-axis, and what would normally be the y-axis is being drawn
    ## as the x-axis
    ## Plot the axis with only the highest-value default tickmark
    xaxis_ticks_default <- pretty(c(ylim_actual2[2],ylim_actual2[1]))
    xaxis_final_tick <- max(xaxis_ticks_default[xaxis_ticks_default < ylim_actual2[2]])

    xaxis_left_spdf <- grid::grid.xaxis(main=TRUE,name="xaxisLeftSPDF", at=xaxis_final_tick)

    grid::grid.force()
    grid::grid.edit("xaxisLeftSPDF::labels", y=grid::unit(-0.5,"cm"))

    # grid::grid.yaxis(main=FALSE)

    ## Crop the raw x- and y-vals of this SPDF plot so that only the points at
    ## x-values within xyrange are actually plotted
    xvals2 <- spdf_vals2$xvals
    included_pts <- which(xvals2 >= xyrange[1] & xvals2 <= xyrange[2])
    xvals2 <- xvals2[included_pts]
    yvals2 <- spdf_vals2$yvals[,included_pts]

    grid::grid.lines(yvals2[1,], xvals2, default.units="native",
                     gp=grid::gpar(col=colvec[base_pop_idxs[1]], lty=2, cex=1.3,
                             lwd=line_widths[base_pop_idxs[1]]))

    grid::grid.polygon(c(0, min(yvals2[1,]), yvals2[1,], min(yvals2[1,]), 0),
                       c(min(xvals2), min(xvals2), xvals2, max(xvals2), max(xvals2)),
                       default.units="native",
                       gp=grid::gpar(col=NA, fill=colvec[base_pop_idxs[1]], alpha=0.3))

    ## Plot the individuals from includepopnames, if required
    if (!is.null(includepopnames)) {
        include.logprobs <- logprob[which(logprob$pop %in% includepopnames),spdf_vals2$refpopA]
        grid::grid.points(rep(0,length(include.logprobs)), include.logprobs,
                          pch=3, gp=grid::gpar(cex=1.3))
    }

    if (show_quantiles) {
        for(idx in 1:ncol(quantile_mat))
        {
            if (quantile_mat[refpopnames[2],idx] > xyrange[1] &&
                quantile_mat[refpopnames[2],idx] < xyrange[2]) {
                grid::grid.lines(x=seq(ylim_actual2[2],ylim_actual2[1],length.out=101),
                                 y=rep(quantile_mat[refpopnames[2],idx],times=101),
                                 default.units="native",gp=grid::gpar(col=col1, lty=2, lwd=1.5))
            }
        }
    }

    if (!only_plot_baseline_pop)
    {
        # first row of yvals is popA in popA, so ignore that one
        grid::grid.lines(yvals2[2,], xvals2,default.units="native",
                         gp=grid::gpar(col=colvec[base_pop_idxs[2]], cex=1.3,
                                 lwd=line_widths[base_pop_idxs[2]]))

        grid::grid.polygon(c(0, min(yvals2[2,]), yvals2[2,], min(yvals2[2,]), 0),
                           c(min(xvals2), min(xvals2), xvals2, max(xvals2), max(xvals2)),
                           default.units="native",
                           gp=grid::gpar(col=NA, fill=colvec[base_pop_idxs[2]],alpha=0.3))
    }

    ## navigate back up to the level of the overall viewport
    grid::upViewport()

    ### Part 4 : the bottom SPDF plot ------------------------------------------
    ### This plot has refpop 1 as the baseline pop
    base_pop_idxs <- spdf_vals1$base_pop_idxs

    if (is.null(plot_params$ylim_input))
    {
        yvals_vec <- as.vector(spdf_vals1$yvals)
        ylim_actual1 <- range(yvals_vec[which(!is.nan(yvals_vec) & !is.infinite(yvals_vec))])
        ylim_actual1[1] <- 0
    } else ylim_actual1 <- plot_params$ylim_input

    xlim_actual <- xyrange

    vp4 <- grid::viewport(layout.pos.row = show_title+1+1, layout.pos.col = 3,
                          xscale = xlim_actual, yscale = ylim_actual1, gp=grid::gpar(cex=1))
    grid::pushViewport(vp4)

    grid::grid.rect(gp=grid::gpar(fill="transparent"))

    ## Plot the x-axis below the bottom SPDF plot
    grid::grid.xaxis(main=TRUE,name="xaxisBottomSPDF")
    grid::grid.force()
    grid::grid.edit("xaxisBottomSPDF::labels", y=grid::unit(-0.5,"cm"))

    ## Plot the y-axis for the bottom SPDF plot, using only the highest tickmark
    yaxis_ticks_default <- pretty(ylim_actual1)
    yaxis_final_tick <- max(xaxis_ticks_default[xaxis_ticks_default < ylim_actual1[2]])

    yaxis_left_spdf <- grid::grid.yaxis(main=TRUE,name="yaxisBottomSPDF", at=yaxis_final_tick)
    grid::grid.force()
    grid::grid.edit("yaxisBottomSPDF::labels", x=grid::unit(-0.4,"cm"))

    ## Crop the x- and y-vals of this SPDF plot so that only the points at x-values
    ## within xyrange are actually plotted
    xvals1 <- spdf_vals1$xvals
    included_pts <- which(xvals1 >= xyrange[1] & xvals1 <= xyrange[2])
    xvals1 <- xvals1[included_pts]
    yvals1 <- spdf_vals1$yvals[,included_pts]

    grid::grid.lines(xvals1, yvals1[1,], default.units="native",
                     gp=grid::gpar(col=colvec[base_pop_idxs[1]], lty=2, cex=1.3,
                             lwd=line_widths[base_pop_idxs[1]]))

    grid::grid.polygon(c(min(xvals1), min(xvals1), xvals1, max(xvals1), max(xvals1)),
                       c(0, min(yvals1[1,]), yvals1[1,], min(yvals1[1,]),0),
                       default.units="native",
                       gp=grid::gpar(col=NA, fill=colvec[base_pop_idxs[1]], alpha=0.3))

    ## Plot the individuals from includepopnames, if required
    if (!is.null(includepopnames))
    {
        include.logprobs <- logprob[which(logprob$pop %in% includepopnames),spdf_vals1$refpopA]
        grid::grid.points(include.logprobs, rep(0,length(include.logprobs)),
                          pch=3, gp=grid::gpar(cex=1.3))
    }

    if (show_quantiles)
    {
        for(idx in 1:ncol(quantile_mat))
        {
            if (quantile_mat[refpopnames[1],idx] > xyrange[1] &&
                quantile_mat[refpopnames[1],idx] < xyrange[2]) {
                grid::grid.lines(x=rep(quantile_mat[refpopnames[1],idx],times=101),
                                 y=seq(ylim_actual1[1],ylim_actual1[2],length.out=101),
                                 default.units="native",gp=grid::gpar(col=col1, lty=2, lwd=1.5))
            }
        }
    }

    if (!only_plot_baseline_pop) {
        # first row of yvals is popA in popA, so ignore that one
        grid::grid.lines(xvals1,yvals1[2,], default.units="native",
                         gp=grid::gpar(col=colvec[base_pop_idxs[2]], cex=1.3,
                                 lwd=line_widths[base_pop_idxs[2]]))

        grid::grid.polygon(c(min(xvals1), min(xvals1), xvals1, max(xvals1), max(xvals1)),
                           c(0,min(yvals1[2,]), yvals1[2,], min(yvals1[2,]),0),
                           default.units="native",
                           gp=grid::gpar(col=NA, fill=colvec[base_pop_idxs[2]],alpha=0.3))
    }

    ## navigate back up to the level of the overall viewport
    grid::upViewport()

    ### Part 5 : "y-axis" label for the rotated left-hand SPDF plot ------------
    ### i.e. what would normally be the y-axis but is now being displayed as the
    ### x-axis of that subplot
    spdf2_label <- paste('PDFs')
    # grid::grid.text(spdf2_label,x=grid::unit(0.5,"cm"), y=grid::unit(0.5,"cm"),
    #           vp=grid::viewport(layout.pos.row=show_title+1+1, layout.pos.col=2))
    grid::grid.text(spdf2_label,x=grid::unit(ylim_actual2[2]-diff(ylim_actual2)/10,"native"),
                    y=grid::unit(xyrange[1]+diff(xyrange)/30,"native"),just="left",
                    vp=grid::viewport(layout.pos.row=show_title+1,layout.pos.col=2,
                                      xscale=c(ylim_actual2[2],ylim_actual2[1]),yscale=xyrange))

    ### Part 6 : the y-axis label for the bottom SPDF plot, --------------------
    ### noting that since this appears in the space to the left, the xscale of
    ### its grid::viewport is the "ylim_actual2" scale from the left-hand SPDF plot
    spdf1_label <- paste('PDFs')
    grid::grid.text(spdf1_label, x=grid::unit(ylim_actual2[1]+diff(ylim_actual2)/6.5,"native"),
                    y=grid::unit(diff(ylim_actual1)/4,"native"),rot=90,
                    vp=grid::viewport(layout.pos.row=show_title+1+1, layout.pos.col=2,
                                      xscale=c(ylim_actual2[2],ylim_actual2[1]),yscale=ylim_actual1))

    ### Part 7 : the y-axis label to the right of the main GenePlot ------------
    grid::grid.text(y_lgp_label, x=grid::unit(1.5,"cm"), rot=90,
                    vp=grid::viewport(layout.pos.row=show_title+1, layout.pos.col=4))

    ### Part 8 : the x-axis label below the bottom SPDF plot -------------------
    grid::grid.text(x_lgp_label, y=grid::unit(0.5,"cm"),
                    vp=grid::viewport(layout.pos.row=show_title+1+1+1, layout.pos.col=3))

    ### Part 2b : the below-Geneplot legend ------------------------------------
    if (show_legend & show_legend_below) {

        grid::grid.legend(labels=refpop_display_names,
                          ncol=min(4, nrefpop), nrow=ceiling(nrefpop/4),
                          pch=refpop_shapevec, vgap=grid::unit(0.7,"lines"),
                          gp=grid::gpar(col=plot_params$rimvec, fill=refpop_colvec, cex=1, fontsize=12),
                          vp=grid::viewport(layout.pos.row = show_title+3+1, layout.pos.col = 3))

        if (show_include) {
            grid::grid.legend(labels=includepop_display_names,
                              ncol=min(4, nincludepop), nrow=ceiling(nincludepop/4),
                              pch=includepop_shapevec, vgap=grid::unit(0.7,"lines"),
                              gp=grid::gpar(col=plot_params$rimvec, fill=includepop_colvec, cex=1, fontsize=12),
                              vp=grid::viewport(layout.pos.row = show_title+3+2, layout.pos.col = 3))
        }
    }

    ### Part 9 : the overlap area values ---------------------------------------
    if (show_overlap_areas)
    {
        grid::grid.text(paste("Area of overlap of ", refpopnames[2], " for baseline pop ",
                              refpopnames[1]," = ", round(spdf_vals1$oavals[2],6)),
                        y=grid::unit(0.5,"cm"),
                        vp=grid::viewport(
                            layout.pos.row = show_title+3+show_legend*show_legend_below*(1+show_include)+1,
                            layout.pos.col = 3),
                        gp=grid::gpar(fontsize=12, fontface="bold"))
        grid::grid.text(paste("Area of overlap of ", refpopnames[1], " for baseline pop ",
                              refpopnames[2]," = ", round(spdf_vals2$oavals[2],6)),
                        y=grid::unit(0.5,"cm"),
                        vp=grid::viewport(
                            layout.pos.row = show_title+3+show_legend*show_legend_below*(1+show_include)+2,
                            layout.pos.col = 3),
                        gp=grid::gpar(fontsize=12, fontface="bold"))
    }
}

#' Lower-level function for running the calculation step of \link{extended_geneplot}.
#' Carries out the GenePlot calculations and the saddlepoint distribution
#' calculations.
#'
#' Use before \link{prepare_plot_params} and \link{plot_geneplot_spdfs}.
#'
#' @param dat The data, in a data frame, with two columns labelled as 'id' and
#'   'pop', and with two additional columns per locus. Missing data at any locus
#'   should be marked as '0' for each allele. The locus columns must be labelled
#'   in the format Loc1.a1, Loc1.a2, Loc2.a1, Loc2.a2, etc.
#'   Missing data must be for BOTH alleles at any locus.
#'   See \link{read_genepop_format} for details of how to import Genepop
#'   format data files into the appropriate format.
#'
#' @param refpopnames Character vector of reference population names, that must
#'   match two values in the 'pop' column of \code{dat}. The SPDF methods
#'   currently only work for a pair of baseline populations, so
#'   \code{refpopnames} must be length 2.
#'
#' @param locnames Character vector, names of the loci, which must match the
#'   column names in the data so e.g. if dat has columns id, pop, EV1.a1,
#'   EV1.a2, EV14.a1, EV14.a2, etc. then you could use 'locnames =
#'   c("EV1","EV14") etc. The locnames do not need to be in any particular order
#'   but all of them must be in \code{dat}.
#'
#' @param includepopnames Character vector (default NULL) of population names to
#'   be included in the calculations as comparison populations. The reference
#'   populations are automatically used as comparison populations for each
#'   other, but you can also add additional comparison populations using
#'   \code{includepopnames}. For example, if the reference pops are Pop1 and
#'   Pop2, and you have some new individuals which you have labelled as PopNew,
#'   then use \code{includepopnames=c("PopNew")} to compare those individuals to
#'   Pop1 and Pop2. You can specify the populations in any order, provided that
#'   they are all in \code{dat}.
#'
#' @param quantiles_vec Specify which quantiles to show on the plots, as a
#'   vector of numbers between 0 and 1. They do not have to be ordered. If
#'   NULL, quantiles will not be plotted.
#'
#' @param prior (default "Rannala") String, either "Rannala" or "Baudouin",
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
#' @param leave_one_out (default TRUE) Boolean, indicates whether or not to
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
#' @param saddlepoint (default TRUE) If TRUE, use saddlepoint approximation to
#'      impute Log-Genotype Probability for individual genotypes with missing
#'      data. If not, use an empirical approximation to impute the LGPs.
#'      Defaults to TRUE because the side plots in the extended GenePlot use the
#'      saddlepoint approximation process.
#'
#' @param rel_tol (default NULL) Specify the relative tolerance for the
#'      numerical integration function that is used to calculate the overlap area
#'      and also the normalization constants for the various distributions. The
#'      default value corresponds to the \code{integrate} default i.e.
#'      \code{.Machine$double.eps^0.25}.
#'
#' @param abs_tol (default NULL) Specify the absolute tolerance for the
#'      numerical integration function that is used to calculate the overlap area
#'      and also the normalization constants for the various distributions. The
#'      default value corresponds to the \code{integrate} default i.e.
#'      \code{.Machine$double.eps^0.25}.
#'
#' @param npts (default 1000) Number of values to use when calculating numerical
#'      integrals (for the overlap area measures, and for the normalization
#'      constants of the distributions). Increasing this value will increase the
#'      precision of the numerical integrals but will also increase the
#'      computational cost. Reducing this below 1000 may save some computation
#'      time if you are not too concerned with the precision of the results.
#'
#' @returns A list with the following components:
#'
#'    \code{logprob} GenePlot calculation results: Log-Genotype Probability
#'    values for all individuals with respect to all of the reference populations.
#'    If there are individuals with missing values, their raw LGPs are shown
#'    which are based on the loci that are present, and also the imputed LGPs
#'    for the full set of loci. This output is the same as the output from
#'    \link{calc_logprob}.
#'
#'    \code{spdf_vals1} Saddlepoint distribution approximations with the first
#'    reference population as the baseline. List contains \code{xvals},
#'    \code{yvals}, \code{wvals} and \code{zvals}, the raw distribution curves
#'    for plotting the distributions.
#'    This list also contains \code{oavals}, \code{probvals} and
#'    \code{diffvals}, which are the Overlap Area, Incumbent Selection
#'    Probability and Home Assignment Probability values for the given baseline
#'    population with all the other populations as comparisons. The sub-list
#'    also records the quantile values requested, the name of the given baseline
#'    pop for this sub-list as \code{refpopA}, the indices of the baseline pop
#'    and comparison pop in the reference pops, the name of the other reference
#'    population and the number of other reference populations (always equal to
#'    one).
#'
#'    \code{spdf_vals2} As for \code{spdf_vals1}, but with the second reference
#'    population as the baseline.
#'
#' @export
calc_geneplot_spdfs <- function(dat,refpopnames,locnames,includepopnames=NULL,
                                quantiles_vec,prior,logten=T,
                                saddlepoint=T,leave_one_out=F,
                                rel_tol=NULL,abs_tol=NULL,npts=1000) {

    if (is.null(rel_tol)) rel_tol <- .Machine$double.eps^0.25
    if (is.null(abs_tol)) abs_tol <- .Machine$double.eps^0.5

    if (!is.vector(refpopnames) || length(refpopnames) != 2) stop("refpopnames should be a vector of length 2 giving the 2 refpops. This function only works for 2 refpops.")

    logprob <- calc_logprob(dat=dat, refpopnames=refpopnames,
                            includepopnames=includepopnames,
                            locnames=locnames, prior=prior, logten=logten,
                            saddlepoint=saddlepoint,
                            leave_one_out=leave_one_out,
                            quantiles=quantiles_vec)
    result <- calc_spdfs(dat=dat, refpopnames=refpopnames, locnames=locnames,
                         prior=prior, logten=logten,
                         leave_one_out=leave_one_out,
                         plot_spdfs=FALSE, calc_vecs=TRUE,
                         rel_tol=rel_tol, abs_tol=abs_tol,
                         quantiles_vec=quantiles_vec,
                         npts=npts)

    spdf_vals1 <- result$spdf_results[[1]]
    spdf_vals2 <- result$spdf_results[[2]]

    list(logprob=logprob, spdf_vals1=spdf_vals1, spdf_vals2=spdf_vals2)
}
