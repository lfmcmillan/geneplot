#' Produce GenePlots of the results from \code{\link{calc_logprob}}.
#'
#' Produce GenePlots of the results obtained from running
#' \code{\link{calc_logprob}}, or replot the results obtained from running
#' \code{\link{geneplot}}.
#'
#' @param logprob_results A data frame containing the results of the GenePlot
#'     calculations, as obtained by running \code{\link{geneplot}} or
#'     \code{\link{calc_logprob}}. See those functions for details.
#'
#' @param plot_type (default NULL) Can be used to specify "twopop" or "manypop"
#'     plots. Defaults to "twopop" for 2 reference pops (i.e. 2 pops listed in
#'     \code{refpopnames} in the call to \code{\link{geneplot}} or
#'     \code{\link{calc_logprob}}) and "manypop" for >2 reference pops.
#'
#' @param plot_bars (default FALSE) Specify what type of plot to use for >2
#'     reference populations.
#'     FALSE (default) plots PCA of the outputs from\code{\link{calc_logprob}} i.e.
#'     runs PCA  the log-genotype-probabilities for all the reference pops and
#'     plots two of the PCs (by default, PC1 and PC2).
#'     TRUE plots multiple bar charts, one per reference pop, with all individuals
#'     as bars, coloured according to their original pop. For the bar plots,
#'     individuals that are in one of the reference pops are ordered according
#'     to their Log-Genotype-Probability with respect to their own pop, and that
#'     ordering is then used to display them in all the other bar plots as well,
#'     so that all the bar plots show the individuals in the same order.
#'
#' @param colvec (default=rep(RColorBrewer::brewer.pal(12,"Paired")[c(1:10,12)], npop)[1:npop])
#'     Vector of colours for plotting. The colours correspond to populations
#'     specified in the order of \code{c(refpopnames, includepopnames)}. Thus
#'     the first element of colvec corresponds to the first element of
#'     \code{refpopnames}; the last element of colvec corresponds to the last
#'     element of \code{includepopnames}.
#'     Colours can be specified using rgb objects, hexadecimal codes, or any of
#'     the R colour names (see \url{http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf}
#'     for a PDF of R colours).
#'
#' @param shapevec Vector of shapes for the plotting points.
#'     These are named shapes from the following list:
#'     "Circle", "Square", "Diamond", "TriangleUp", "TriangleDown", "OpenSquare",
#'     "OpenCircle", "OpenTriangleUp", "Plus", "Cross", "OpenDiamond",
#'     "OpenTriangleDown", "Asterisk"
#'     which correspond to the following pch values for R plots:
#'     21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 8.
#'     Do not use the numbers, use the words, which will be automatically
#'     converted within plot_logprob into the appropriate codes.
#'     The elements of shapevec correspond to the populations specified in the
#'     order of \code{c(refpopnames, includepopnames)}. Thus the first element
#'     of shapevec corresponds to the first element of \code{refpopnames}; the
#'     last element of shapevec corresponds to the last element of
#'     \code{includepopnames}.
#'     Defaults to the list above, looping through as many times as required for
#'     all the populations.
#'
#' @param mark_impute (default FALSE) Boolean, indicates whether to mark
#'     individuals with missing data using asterisks.
#'
#' @param txt (default "points") Defines whether to plot individuals as points
#'     on the GenePlot (\code{"points"}), or whether to display the name of
#'     their population (\code{"pop"}), subpopulation (\code{"subpop"}), or ID
#'     (\code{"id"}) as text. For \code{"subpop"}, must have 'subpop' as one of
#'     the columns in the \code{dat} input to calc_logprob. Then 'subpop' will
#'     automatically be included in the \code{logprob_results} object.
#'
#' @param use_legend (default TRUE) Plot the legend (or FALSE for don't plot the legend).
#'
#' @param legend_pos (default "bottomleft") Define where to plot the legend,
#'     uses the same position labels as in the \code{legend} function e.g. "topright".
#'
#' @param xyrange (default NULL) Specify the xyrange as a vector, will be the
#'     same range for both axes. Default is slightly wider than the range of the
#'     calculated Log-Genotype-Probabilities for all individuals in the plot.
#'
#' @param orderpop Specify the plotting order for the populations.
#'     E.g. if orderpop=c("Pop4", "Pop2"), then points for individuals from Pop4
#'     will be plotted first, then individuals from Pop2 will be plotted over the
#'     top of them, etc.
#'     Default is NULL, in which case populations are plotted in order of size,
#'     so the population with the largest number of points is plotted at the
#'     bottom, and the population with the smallest number of individuals/points
#'     is plotted over the top, so as not to be obscured.
#'
#' @param axispop is used when \code{length(refpopnames) == 2} i.e. when
#'     plot_type="twopop".
#'     It is of the form axispop=c(x="Pop1", y="Pop2"), meaning that the Pop1
#'     reference population will be plotted on the x axis, and the Pop2 reference
#'     population will be plotted on the y axis.
#'     Default is NULL, which plots the first population in \code{refpopnames}
#'     on the x axis and the second population in \code{refpopnames} on the y axis.
#'
#' @param axis_labels (default NULL) Used for plots with 2 reference pops.
#'     Character vector, 2 elements, can be used to specify more readable axis
#'     labels. Defaults to the 'pop' labels in \code{logprob_results}.
#'
#' @param short_axis_labels (default FALSE) Used for plots with 2 reference pops.
#'     FALSE (default) gives full-length axis labels of the form
#'     "Log10 genotype probability for population Pop1"
#'     TRUE gives short-form axis labels of the form "LGP10 for population Pop1"
#'
#' @param grayscale_quantiles (default FALSE) Used for plots with 2 reference pops.
#'     FALSE (default) plots the quantile lines using colvec colours
#'     TRUE plots the quantile lines in gray (as the default colours can be
#'     quite pale, the grayscale quantile lines can be easier to see than the
#'     default coloured ones).
#'
#' @param dim1 (default 1) Used for plots with more than 2 reference pops, when
#'     plot_bars=FALSE. Specifies which principal component should be plotted on
#'     x-axis.
#' @param dim2 (default 2) Used for plots with more than 2 reference pops, when
#'     plot_bars=FALSE. Specifies which principal component should be plotted on
#'     y-axis.
#'
#' @param layout_already_set (default=FALSE) Boolean, used for plots with more
#'     than 2 reference pops, when plot_bars=TRUE.
#'     Indicates whether the \code{layout} command, for arranging plots, has
#'     already been called by a higher-level function (or if the \code{par(mfrow)}
#'     command has been called earlier). TRUE prevents the layout command within
#'     \code{plot_logprob_barplot} from clashing with the higher-level command.
#'
#' @param cexpts (default 1.4) Specify the size of the points in the plot.
#'
#' @return Displays the plot.
#'
#' @references McMillan, L. and Fewster, R. "Visualizations for genetic
#'   assignment analyses using the saddlepoint approximation method" (2017)
#'   \emph{Biometrics}.
#'
#'   Rannala, B., and Mountain, J. L. (1997). Detecting
#'   immigration by using multilocus genotypes. \emph{Proceedings of the
#'   National Academy of Sciences} \strong{94}, 9197--9201.
#'
#'   Piry, S., Alapetite, A., Cornuet, J.-M., Paetkau, D., Baudouin, L., and
#'   Estoup, A. (2004). GENECLASS2: A software for genetic assignment and
#'   first-generation migrant detection. \emph{Journal of Heredity} \strong{95},
#'   536--539.
#'
#' @author Log-Genotype-Probability calculations based on the method of Rannala
#' and Mountain (1997) as implemented in GeneClass2, updated to allow for individuals
#' with missing data and to enable accurate calculations of quantiles of the
#' Log-Genotype-Probability distributions of the reference populations.
#' See McMillan and Fewster (2017) for details.
#'
#' @importFrom graphics abline axis barplot layout legend mtext par plot points text
#'
#' @export
plot_logprob <- function(logprob_results,
                         plot_type=switch(as.character(length(refpopnames)), "2"="twopop", "manypop"),
                         plot_bars=F,
                         colvec=NA,
                         shapevec=NA,
                         mark_impute=F,
                         txt="points",
                         use_legend=T,
                         legend_pos="bottomleft",
                         xyrange=NULL,
                         orderpop=NULL,
                         axispop=NULL,
                         axis_labels=NULL,
                         short_axis_labels=F,
                         grayscale_quantiles=F,
                         dim1=1, dim2=2,
                         layout_already_set=F,
                         cexpts=1.4)
{
    includepopnames <- attributes(logprob_results)$includepopnames
    refpopnames <- attributes(logprob_results)$refpopnames
    if (is.null(refpopnames)) refpopnames <- rownames(attributes(logprob_results)$qmat)
    if (is.null(refpopnames)) refpopnames <- rownames(attributes(logprob_results)$allele_freqs[[1]])
    if (is.null(refpopnames)) stop("Error! No refpops defined as attributes of logprob_results or in qmat or allele freqs attributes!")
    allpopnames <- unique(c(refpopnames, includepopnames))
    npop <- length(allpopnames)
    posterior_nu_list <- attributes(logprob_results)$allele_freqs
    saddlepoint <- attributes(logprob_results)$saddlepoint
    logten <- attributes(logprob_results)$logten
    quantile_mat <- attributes(logprob_results)$qmat

    if (is.null(plot_type))
    {
        plot_type <- switch(as.character(length(refpopnames)), "2"="twopop", "manypop")
    }

    ## Error checking
    if (plot_type=="twopop")
    {
        if(!is.null(axispop) &&
           (!is.character(axispop) || !is.vector(axispop) ||
            names(axispop) != c("x","y"))) stop("axispop must be specified in the form c(x='Pop1',y='Pop2')")
    }

    ## Shapes and colours ------------------------------------------------------
    shapevec_default <- rep(21:25, npop)[1:npop]
    if(length(shapevec)>1){
        ## If some shapes are supplied then shapevec has length npop, otherwise shapevec=NA (length 1)
        shapes_ref <- c("Circle", "Square", "Diamond", "TriangleUp",
                        "TriangleDown", "OpenSquare", "OpenCircle",
                        "OpenTriangleUp", "Plus", "Cross", "OpenDiamond",
                        "OpenTriangleDown", "Asterisk")
        numbers_ref <- c(21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 8)
        shapes_inds <- match(shapevec, shapes_ref)
        shapevec <- numbers_ref[shapes_inds]
        shapevec[is.na(shapevec)] <- shapevec_default[is.na(shapevec)]
    }
    else shapevec <- shapevec_default

    ## Note: using the RColorBrewer brewer.pal() function NOT the grDevices
    ## palette() function because the latter only has 10 colours in the "Paired"
    ## palette so loops in sync with the 5 default shapes, whereas the original
    ## RColorBrewer "Paired" palette has 12 colours
    colvec_default <- rep(RColorBrewer::brewer.pal(12,"Paired")[c(1:10,12)], npop)[1:npop]
    ## If some colours are supplied then colvec has length npop, otherwise colvec=NA (length 1)
    if(length(colvec)>1) colvec[is.na(colvec)] <- colvec_default[is.na(colvec)]
    else colvec <- colvec_default

    ## In case xyrange is entered the wrong way round by the user:
    if(!is.null(xyrange)) xyrange <- sort(xyrange)

    ## Populations and loci ----------------------------------------------------

    ## Order the populations before plotting
    if(is.null(orderpop)){
        ## If no plotting order is specified (orderpop not supplied), arrange allpopnames to be
        ## in decreasing order of sample size:
        allpop_sizes <- numeric(length(npop))
        for(pp in 1:npop) allpop_sizes[pp] <- length(logprob_results$pop[logprob_results$pop==allpopnames[pp]])
        order_inds <- order(allpop_sizes, decreasing=T)
    }
    ## If orderpop is supplied, order_inds is gained from matching orderpop straight to allpopnames:
    else order_inds <- match(orderpop, allpopnames)
    allpopnames <- allpopnames[order_inds]

    ## Reorder shapevec and colvec by the same ordering to keep the correct assignment of
    ## colours and shapes to populations:
    shapevec <- shapevec[order_inds]
    colvec <- colvec[order_inds]

    ## If axispop is not specified, and there are two reference populations, make axispop the same as refpopnames:
    if(is.null(axispop) & length(refpopnames)==2){
        axispop <- refpopnames
        names(axispop) <- c("x", "y")
    }

    ## Call specific plot functions --------------------------------------------
    if (plot_bars)
    {
        if (saddlepoint)
        {
            plot_logprob_barplot(logprob=logprob_results, refpopnames=refpopnames,
                                 allpopnames=allpopnames, colvec=colvec, logten=logten,
                                 quantile_mat=quantile_mat, all_loci_sim_logprob=NULL,
                                 posterior_nu_list=posterior_nu_list,
                                 axis_labels=axis_labels, use_legend=use_legend,
                                 layout_already_set=layout_already_set)
        }
        else
        {
            all_loci_sim_logprob <- attributes(logprob_results)$all_loci_sim_logprob

            plot_logprob_barplot(logprob=logprob_results, refpopnames=refpopnames,
                                 allpopnames=allpopnames, colvec=colvec,
                                 logten=logten, quantile_mat=quantile_mat,
                                 all_loci_sim_logprob=all_loci_sim_logprob,
                                 posterior_nu_list=NULL, axis_labels=axis_labels,
                                 use_legend=use_legend, layout_already_set=layout_already_set)
        }
    }
    else
    {
        if(plot_type=="twopop") plot_logprob_twopop(logprob=logprob_results,
                                                    axispop=axispop,
                                                    allpopnames=allpopnames, txt=txt,
                                                    shapevec=shapevec, colvec=colvec,
                                                    legend_pos=legend_pos,
                                                    mark_impute=mark_impute,
                                                    logten=logten,
                                                    quantile_mat=quantile_mat,
                                                    xyrange=xyrange, cexpts=cexpts,
                                                    axis_labels=axis_labels,
                                                    short_axis_labels=short_axis_labels,
                                                    use_legend=use_legend,
                                                    grayscale_quantiles=grayscale_quantiles)
        else plot_logprob_manypop(logprob=logprob_results, refpopnames=refpopnames,
                                  allpopnames=allpopnames, txt=txt,
                                  shapevec=shapevec, colvec=colvec,
                                  legend_pos=legend_pos, mark_impute=mark_impute,
                                  dim1=dim1, dim2=dim2, xyrange=xyrange,
                                  cexpts=cexpts, use_legend=use_legend)
    }
}

plot_logprob_twopop <- function(logprob, axispop, allpopnames, txt, shapevec, colvec,
                                legend_pos, mark_impute, logten, quantile_mat, xyrange,
                                cexpts, use_legend=T, axis_labels=NULL,
                                short_axis_labels=F, grayscale_quantiles=F){

    ## If quantile_mat=NULL, then no significance crosslines are plotted.
    ## Otherwise, quantile_mat is a matrix of values, e.g.
    ##          1%   99%
    ## Pop1   -37.3 -25.7
    ## Pop2   -29.8 -11.9
    ## giving the quantiles of the corresponding distributions of genotype
    ## probabilities.  The base of the logs will already be log_10 if logten=T,
    ## otherwise it will be log_e.

    ## Extract refpopnames from axispop, in the required order.  For example,
    ## if axispop=c(y="Pop1", x="Pop2")
    ## then refpopnames = c("Pop1", "Pop2")
    refpopnames <- axispop[c("x", "y")]
    if(length(refpopnames)!=2) stop("Need two reference populations for plot_logprob_twopop.")

    ## Set up for plotting:
    if(is.null(xyrange)) xyrange <- range(c(logprob[, refpopnames[1]], logprob[, refpopnames[2]])) +
            c(-0.05, 0.05) * diff(range(c(logprob[, refpopnames[1]], logprob[, refpopnames[2]])))

    if (is.null(axis_labels) || !is.vector(axis_labels) ||
        !is.character(axis_labels) || length(axis_labels) != 2)
    {
        if (logten && short_axis_labels) label <- "LGP10 for population"
        else if (logten && !short_axis_labels) label <- "Log10 genotype probability for population"
        else if (!logten && short_axis_labels) label <- "LGP for population"
        else if (!logten && !short_axis_labels) label <- "Log genotype probability for population"

        xlabel <- paste(label, refpopnames[1])
        ylabel <- paste(label, refpopnames[2])
    }
    else
    {
        xlabel <- axis_labels[1]
        ylabel <- axis_labels[2]
    }

    plot(1, xlim=xyrange, ylim=xyrange, type="n", xlab=xlabel, ylab=ylabel,
         cex.lab=1, cex.main=1, cex.axis=1)

    split_dat <- split(logprob, logprob$pop)

    ## Sort out colours if plotting points.  Any populations that have pch between
    ## 21 and 25 should have rimvec set to 1 for the rim:
    rimvec <- colvec
    if(txt=="points" | txt=="id") rimvec[21 <= shapevec & shapevec <= 25] <- "black"

    plot_single_pop <- function(whichpop){
        ## Plot results for a single population.  If whichpop=1 then plot the
        ## first population in allpopnames, etc.
        pop <- allpopnames[whichpop]
        popdat <- split_dat[[pop]]
        p1 <- popdat[[refpopnames[1]]]
        p2 <- popdat[[refpopnames[2]]]
        switch(txt,
               pop = text(p1, p2, pop, col=colvec[whichpop]),
               subpop = text(p1, p2, popdat$subpop, col=colvec[whichpop]),
               id = text(p1, p2, popdat$id, col=colvec[whichpop]),
               points = points(p1, p2, col=rimvec[whichpop],
                               pch=shapevec[whichpop],
                               bg=colvec[whichpop], cex=cexpts),
               stop("txt should be pop, subpop, id, or points"))

    }

    npop <- length(allpopnames)
    sapply(1:npop, plot_single_pop)

    abline(0, 1, col="black", lwd=2)
    ## Next two lines if wanting to show the Geneclass scores = 1.0 lines:
    abline(-1, 1, col="darkgrey", lty=1, lwd=1)
    abline(1, 1, col="darkgrey", lty=1, lwd=1)

    ## Add crosslines for quantiles if required --------------------------------
    if(!is.null(quantile_mat)){

        if (grayscale_quantiles)
        {
            axis(side=3, at=quantile_mat[refpopnames[1],], labels=colnames(quantile_mat),
                 col.axis="black", tick=F, cex.axis=1, line=-0.8)
            abline(v=quantile_mat[refpopnames[1],], col="black", lty=2, lwd=1.5)
            axis(side=4, at=quantile_mat[refpopnames[2],], labels=colnames(quantile_mat),
                 col.axis="black", tick=F, cex.axis=1, line=-0.8, las=1)
            abline(h=quantile_mat[refpopnames[2],], col="black", lty=2, lwd=1.5)
        }
        else
        {
            col1 <- colvec[which(allpopnames==refpopnames[1])]
            col2 <- colvec[which(allpopnames==refpopnames[2])]
            axis(side=3, at=quantile_mat[refpopnames[1],], labels=colnames(quantile_mat),
                 col.axis=col1, tick=F, cex.axis=1, line=-0.8)
            abline(v=quantile_mat[refpopnames[1],], col=col1, lty=2, lwd=1.5)
            axis(side=4, at=quantile_mat[refpopnames[2],], labels=colnames(quantile_mat),
                 col.axis=col2, tick=F, cex.axis=1, line=-0.8, las=1)
            abline(h=quantile_mat[refpopnames[2],], col=col2, lty=2, lwd=1.5)
        }
    }

    ## Add legend if needed ----------------------------------------------------
    if((txt=="points" | txt=="id") && use_legend){
        legend(legend_pos, pch=shapevec[1:npop], col=rimvec[1:npop], pt.bg=colvec[1:npop],
               legend=allpopnames, cex=1, bg="white")
    }

    ## Flag the points gained by imputation if required:
    if(mark_impute)
        points(logprob[logprob$status=="impute", refpopnames[1]],
               logprob[logprob$status=="impute", refpopnames[2]], pch="*")
}

plot_logprob_manypop <- function(logprob, refpopnames, allpopnames, txt,
                                 shapevec, colvec, legend_pos, mark_impute,
                                 dim1, dim2, xyrange, cexpts, use_legend=T){

    ## PCA ---------------------------------------------------------------------
    ## Do the PCA on the matrix of imputed log-likelihoods, and set up plot:

    ## catline("Colours", colvec)
    pca.mat <- prcomp(logprob[refpopnames])
    if(is.null(xyrange)) xyrange <- range(pca.mat$x[ , c(dim1, dim2)]) +
                                    c(-0.2, 0.05)*diff(range(pca.mat$x[ , c(dim1, dim2)]))
    plot(1, xlim=xyrange, ylim=xyrange, type="n", cex.lab=1.3,
         cex.main=1.4, cex.axis=1.2, xlab="", ylab="")

    ## Join the PCA matrix onto the columns of logprob.  The new columns are
    ## called PC1, PC2, etc.
    logprob <- cbind(logprob, pca.mat$x)
    PC1_name <- paste("PC", dim1, sep="")
    PC2_name <- paste("PC", dim2, sep="")

    split_dat <- split(logprob, logprob$pop)

    ## Sort out colours if plotting points.  Any populations that have pch
    ## between 21 and 25 should have rimvec set to 1 for the rim:
    rimvec <- colvec
    if(txt=="points" | txt=="id") rimvec[21 <= shapevec & shapevec <= 25] <- "black"

    plot_single_pop <- function(whichpop){
        ## Plot results for a single population.  If whichpop=1 then plot the
        ## first population in allpopnames, etc.
        pop <- allpopnames[whichpop]
        popdat <- split_dat[[pop]]
        p1 <- popdat[[PC1_name]]
        p2 <- popdat[[PC2_name]]
        switch(txt,
               pop = text(p1, p2, pop, col=colvec[whichpop]),
               subpop = text(p1, p2, popdat$subpop, col=colvec[whichpop]),
               id = text(p1, p2, popdat$id, col=colvec[whichpop]),
               points = points(p1, p2, col=rimvec[whichpop], pch=shapevec[whichpop],
                               bg=colvec[whichpop], cex=cexpts),
               stop("txt should be pop, subpop, id, or points"))
    }

    npop <- length(allpopnames)
    sapply(1:npop, plot_single_pop)

    if((txt=="points" | txt=="id") && use_legend) {
        legend(legend_pos, pch=shapevec[1:npop], col=rimvec[1:npop],
               pt.bg=colvec[1:npop], legend=allpopnames, cex=1.2)
    }

    ## Find the percentage variances explained by each of the axes.
    ## This technique (squaring the sdevs and taking overall proportion)
    ## has been confirmed to be correct by Brian McArdle, for the PCA.

    perc.var <- pca.mat$sdev^2/sum(pca.mat$sdev^2)
    perc.var1 <- format(perc.var[dim1], digits=2)
    perc.var2 <- format(perc.var[dim2], digits=2)
    perc.var12 <- format(perc.var[dim1] + perc.var[dim2], digits=2)
    #         title(sub=paste("%var: ", PC1_name, "=", perc.var1, "; ", PC2_name, "=", perc.var2, "; Tot=",
    #               perc.var12, sep=""), cex.sub=1.5, line=2.5)
    mtext(paste(PC1_name, "%var",perc.var1), side=1, line=2, cex=1.3)
    mtext(paste(PC2_name, "%var",perc.var2), side=2, line=2, cex=1.3)

    ## Flag the points gained by imputation if required:
    if(mark_impute)
        points(logprob[logprob$status=="impute", PC1_name],
               logprob[logprob$status=="impute", PC2_name], pch="*")
}

plot_logprob_barplot <- function(logprob, refpopnames, allpopnames, colvec, logten,
                                 quantile_mat, posterior_nu_list=NULL,
                                 all_loci_sim_logprob=NULL, axis_labels=NULL,
                                 use_legend=T, layout_already_set=F){
    ## bars.plot 29/04/2015 (based on earlier version from Masters project Rfunc_bars.R)
    ## Plot series of bar charts, one for each reference population, showing the
    ## log-probabilities for every single individual (ref pop individuals and
    ## other individuals), with different coloured bars for different
    ## populations
    ## Logs in base 10 if logten=T, or base e if logten=F
    ## If quantile_mat=NULL, then no significance crosslines are plotted.
    ## Otherwise, quantile_mat is a matrix of values, e.g.
    ##              1%   99%
    ## Pop1   -37.3 -25.7
    ## Pop2 -29.8 -11.9
    ## giving the quantiles of the corresponding distributions of genotype
    ## probabilities.
    ## If the quantiles are included, then they will be plotted as horizontal
    ## lines on top of the bar charts for each population

    if (!is.null(axis_labels) && (!is.vector(axis_labels) || !is.character(axis_labels)
                                  || length(axis_labels) != length(refpopnames)))
    {
        stop("axis_labels must be a character vector, the same length as refpopnames")
    }

    split_dat <- split(logprob, logprob$pop)

    bgvec <- colvec

    ## Set up layout of plots and legend
    if (!layout_already_set)
    {
        if (use_legend) layout(matrix(1:(length(refpopnames)+1), ncol=1),
                               heights=c(rep(2,times=length(refpopnames)),1))
        else layout(matrix(1:(length(refpopnames)), ncol=1))
    }

    ## Reduce the amount of whitespace between the plots
    if (!layout_already_set) par(mar=c(1.5,3,2.5,0))

    ## Get sorting order for individuals in reference populations
    ## Sort them in order of log-probs for their own population
    indiv_order = list()

    for (i in 1:length(refpopnames))
    {
        popdat = split_dat[[refpopnames[i]]]
        probs = popdat[[refpopnames[i]]]
        indiv_order[[i]] = order(probs, decreasing=TRUE)
    }

    all_heights = matrix(0, nrow = nrow(logprob), ncol = length(refpopnames))
    all_colours = matrix(0, nrow = nrow(logprob), ncol = length(refpopnames))

    for (i in 1:length(refpopnames))
    {
        heights = vector(mode="numeric")
        colours = vector(mode="numeric")

        ## Set up bars for all individuals, coloured by population
        for (j in 1:length(allpopnames))
        {
            pop <- allpopnames[j]
            popdat <- split_dat[[pop]]
            this_refpop_probs <- popdat[[refpopnames[i]]]

            ## If this is one of the ref pops, order the individuals by their log
            ## probabilities with respect to their own population
            ## i.e. order "Taik" individuals by their log probs for the "Taik" population
            ref_pop_idx = match(pop, refpopnames)
            if (!is.na(ref_pop_idx)) this_refpop_probs <- this_refpop_probs[indiv_order[[ref_pop_idx]]]

            ## If the simulated individual data or the saddlepoint stuff is available,
            ## convert to using CDF values instead
            if (!is.null(all_loci_sim_logprob))
            {
                ECDF <- ecdf(all_loci_sim_logprob[refpopnames[i],])

                this_refpop_cdf_values = ECDF(this_refpop_probs)

                ## Any heights that are too small should be replaced with zeros
                heights[which(heights < 1E-3)] <- 0

                heights <- c(heights, this_refpop_cdf_values)

            }
            else if (!is.null(posterior_nu_list))
            {
                ## Make sure the saddlepoint function is being calculated for the refpop,
                ## not the sample population being looked at!
                ## Extract the posterior allele frequencies for this pop, keeping
                ## them as a LIST
                post_nu_pop <- lapply(posterior_nu_list, function(x) x[refpopnames[i],])
                post_info <- calc.multi.locus.probs.func(post_nu_pop)
                multi_K_params <- make.K.params(post_info$dist)
                mean_pop <- mu(multi_K_params)
                this_refpop_cdf_values <- sapply(this_refpop_probs, Fhat,
                                                  post_info, mean.pop=mean_pop, logten=logten)

                ## Any heights that are too small should be replaced with zeros
                heights[which(heights < 1E-3)] <- 0

                heights <- c(heights, this_refpop_cdf_values)
            }
            else
            {
                heights <- c(heights, this_refpop_probs)
            }

            colours <- c(colours, rep(colvec[j],times=length(this_refpop_probs)))
        }

        # Use 100*heights to get readable percent values
        all_heights[,i] <- heights*100
        all_colours[,i] <- colours
    }

    min.all.pops <- min(all_heights, na.rm = TRUE)

    for (i in 1:length(refpopnames))
    {
        heights <- all_heights[,i]
        colours <- all_colours[,i]

        ## Set up the individual plot title
        if (is.null(axis_labels)) title <- paste0("Percentiles for population ",
                                                  refpopnames[i], ", by individual")
        else title <- axis_labels[i]

        if (!is.null(all_loci_sim_logprob))
        {
            barplot(height = heights, col = colours, ylim=c(0,100))
            mtext(title,side=1, line=1, cex=0.7)
        }
        else
        {
            barplot(height = heights, col = colours, ylim=c(min.all.pops,100))
            mtext(title, side=1, line=1, cex=0.7)
        }

        ## Add quantile lines to the plot
        if (!is.null(quantile_mat))
        {
            if (nrow(logprob) <= 30) horiz_label_pos <- 0.5
            else horiz_label_pos <- 2
            quantile_labels <- colnames(quantile_mat)
            quantile_positions <- as.numeric(unlist(strsplit(quantile_labels,"%")))
            for (q in 1:length(quantile_labels)) {
                abline(h=quantile_positions[q], lty=2, lwd=2, col="black")
                if (quantile_positions[q] < 20) text(horiz_label_pos, quantile_positions[q]+4, quantile_labels[q], pos=2, font=2)
                else text(horiz_label_pos, quantile_positions[q]-3, quantile_labels[q], pos=2, font=2)
            }
        }
    }

    if (use_legend)
    {
        ## Increase the amount of whitespace around the legend plot
        if (!layout_already_set) par(mar=c(2,2,3,2))

        ## Add a legend plot
        if (length(refpopnames) <= 3) barHeight <- 1
        else barHeight <- 2
        barplot(height = rep(barHeight,times=length(allpopnames)), col = colvec,
                names.arg = allpopnames, horiz = FALSE, axes=FALSE, ylim=c(0,1.2))
    }

    if (!layout_already_set) {
        par(mfrow=c(1,1))
        par(mar=c(5,4,4,2)+0.1)
    }
}
