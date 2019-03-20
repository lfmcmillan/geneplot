#' Run GenePlot and produce results and plots for the selected populations.
#'
#' All individuals in \code{dat} whose population matches one of \code{refpopnames}
#' will be included. Leave-one-out will be used when calculating the log-genotype
#' probability for individual with respect to their own reference population,
#' if specified in the inputs (default is NON leave-one-out).
#'
#' Default is NON leave-one-out but WE STRONGLY RECOMMEND USING LEAVE-ONE-OUT,
#' ESPECIALLY FOR SMALL SAMPLES (<30).
#'
#' \code{includepopnames} specifies which additional populations are plotted.
#' It defaults to NULL, in which case only \code{refpopnames} will be plotted.
#' \code{includepopnames} can specify the populations in any order.
#' If \code{refpopnames} are missing from \code{includepopnames}, they will be added.
#'
#' NOTE that if a population is not in \code{includepopnames}/\code{refpopnames},
#' then any alleles private to that population will NOT be included in the
#' prior / posterior.  Thus the posterior for a given refpop will change slightly
#' depending on which populations are in \code{includepopnames}/\code{refpopnames}.
#'
#' @param plotit (default=TRUE) is whether to produce a plot, or whether to just
#' do the calculations and spit out the table of log-genotype probabilities.
#' FALSE just runs calc_logprob but not plot_logprob, so is equivalent to a call
#' to calc_logprob.
#'
#' @inherit calc_logprob
#' @inheritParams plot_logprob
#'
#' @examples
#' ## Example dataset created directly within R (usually you would read in data from a file instead):
#' ratLocnames <- c("D10Rat20","D11Mgh5","D15Rat77","D16Rat81","D18Rat96",
#'                  "D19Mit2","D20Rat46","D2Rat234","D5Rat83","D7Rat13")
#' ratData <- rbind(
#' c("Ki001","Kai",96,128,246,280,234,250,155,165,226,232,219,231,149,149,101,127,174,176,164,182),
#' c("Ki002","Kai",122,126,246,276,238,238,155,165,226,232,223,231,187,187,107,121,174,174,164,164),
#' c("Ki003","Kai",122,122,276,280,234,234,157,165,244,244,231,231,187,187,107,107,174,174,164,182),
#' c("Ki004","Kai",130,130,276,280,238,238,157,165,0,0,223,231,187,187,101,111,168,176,184,184),
#' c("Ki009","Kai",122,122,276,276,234,236,165,165,240,244,229,231,187,187,89,101,174,176,164,164),
#' c("Ki010","Kai",122,122,278,280,236,236,155,165,236,244,219,231,185,187,101,101,168,174,164,164),
#' c("Ki011","Kai",120,128,280,282,236,238,155,165,226,236,223,231,149,149,99,101,174,174,164,164),
#' c("Bi01","Brok",96,126,280,280,236,250,165,165,232,246,231,231,185,187,89,89,170,176,154,164),
#' c("Bi02","Brok",96,126,280,280,250,262,155,155,232,232,231,233,149,185,127,127,174,174,164,166),
#' c("Bi03","Brok",96,126,280,280,258,262,165,165,232,232,231,231,185,187,89,127,174,174,164,164),
#' c("Bi04","Brok",96,126,280,280,238,262,155,155,232,232,231,233,149,185,127,127,174,174,164,164),
#' c("Bi05","Brok",96,122,280,280,250,258,155,155,226,244,231,231,187,187,107,127,174,176,164,164),
#' c("Bi06","Brok",96,96,280,280,238,262,155,155,232,232,231,231,187,187,123,127,174,174,164,164),
#' c("Bi11","Brok",96,96,278,280,234,250,165,165,226,240,231,231,149,187,89,99,170,170,154,164),
#' c("Bi12","Brok",96,96,276,280,234,250,165,165,240,240,231,231,187,187,89,99,170,174,154,164),
#' c("Bi13","Brok",96,126,276,276,246,250,165,165,226,244,231,231,149,187,99,99,174,174,164,164),
#' c("Bi14","Brok",96,126,276,276,262,262,155,165,226,244,231,231,149,187,89,107,170,174,154,164),
#' c("Ki092","Main",122,126,280,282,234,238,165,165,236,240,231,231,149,187,95,95,0,0,164,164),
#' c("Ki093","Main",122,126,282,282,238,238,165,165,236,240,231,231,149,187,95,107,166,174,164,182),
#' c("Ki094","Main",122,126,280,282,238,238,165,165,226,240,231,231,173,187,95,127,174,176,154,182),
#' c("Ki095","Main",120,126,280,280,234,236,155,165,244,246,231,231,161,187,123,127,174,174,154,154),
#' c("Ki097","Main",122,126,280,280,236,236,163,165,236,242,219,231,149,161,107,115,166,174,164,166),
#' c("Ki098","Main",96,122,276,280,236,238,155,165,242,244,233,233,149,187,99,107,174,174,164,164),
#' c("Ki100","Main",122,122,280,280,234,234,155,165,236,236,219,235,0,0,107,107,174,176,164,164),
#' c("Ki101","Main",122,126,276,280,234,238,155,155,236,244,229,231,0,0,101,101,0,0,164,182),
#' c("Ki102","Main",122,126,0,0,0,0,155,163,0,0,229,231,0,0,107,107,0,0,0,0),
#' c("Ki103","Main",122,122,280,280,234,236,163,165,0,0,231,233,0,0,99,107,0,0,164,184),
#' c("Ki104","Main",96,126,276,280,236,238,157,165,230,246,231,231,149,187,107,107,0,0,164,164),
#' c("Ki105","Main",122,126,276,280,238,250,157,165,226,244,217,231,0,0,111,121,174,174,164,164),
#' c("R01","Erad10",128,128,280,288,234,244,155,165,242,244,231,231,149,149,107,107,174,174,164,166),
#' c("R02","Erad10",128,130,276,288,238,244,155,155,228,244,223,231,149,149,101,111,174,174,164,166),
#' c("R03","Erad10",128,130,276,288,238,244,155,155,244,244,223,231,149,187,107,111,174,176,164,166))
#' ratData <- as.data.frame(ratData, stringsAsFactors=FALSE)
#' names(ratData) <- c("id","pop","D10Rat20.a1","D10Rat20.a2","D11Mgh5.a1","D11Mgh5.a2",
#'                     "D15Rat77.a1","D15Rat77.a2","D16Rat81.a1","D16Rat81.a2",
#'                     "D18Rat96.a1","D18Rat96.a2","D19Mit2.a1","D19Mit2.a2",
#'                     "D20Rat46.a1","D20Rat46.a2","D2Rat234.a1","D2Rat234.a2",
#'                     "D5Rat83.a1","D5Rat83.a2","D7Rat13.a1","D7Rat13.a2")
#'
#' ## Run GenePlot for 2 reference populations:
#' geneplot(dat=ratData,refpopnames=c("Kai","Main"),locnames=ratLocnames,
#'          prior="Baudouin", leave_one_out=TRUE,
#'          colvec=c("darkorchid4","steelblue"), shapevec=c("Circle","Square"),
#'          axis_labels=c("Log10 genotype probability for Kaikoura Island",
#'                        "Log10 Genotype probability for Mainland"))
#'
#' ## Run GenePlot for 2 reference populations, and include an extra group of
#' ## individuals who are going to be compared to the 2 reference populations to
#' ## see which reference population they are most similar to (note that we
#' ## specify an additional colour and shape for the new individuals, but the
#' ## axis labels stay the same because they correspond to the two reference
#' ## populations):
#' results <- geneplot(dat=ratData,refpopnames=c("Kai","Main"),locnames=ratLocnames,
#'          includepopnames=c("Erad10"), prior="Baudouin", leave_one_out=TRUE,
#'          colvec=c("darkorchid4","steelblue","chartreuse4"),
#'          shapevec=c("Circle","Square","TriangleUp"),
#'          axis_labels=c("Log10 genotype probability for Kaikoura Island",
#'                        "Log10 Genotype probability for Mainland"))
#'
#' ## Barplot:
#' plot_logprob(results, plot_bars=TRUE,
#'          colvec=c("darkorchid4","steelblue","chartreuse4"))
#'
#' ## Rename Kai and Main as a single population and compare that population to Brok:
#' ratData2 <- ratData
#' ratData2$pop[which(ratData2$pop %in% c("Kai","Main"))] <- "KaiMain"
#' geneplot(dat=ratData2,refpopnames=c("KaiMain","Brok"),locnames=ratLocnames,
#'          prior="Rannala", leave_one_out=TRUE,
#'          colvec=c("forestgreen","darkgoldenrod"),
#'          shapevec=c("TriangleDown","Diamond"),
#'          axis_labels=c("Log10 genotype probability for Kaikoura and Mainland",
#'                        "Log10 Genotype probability for Broken Islands"))
#'
#' \donttest{
#' ## Example code for reading in a Genepop-format file
#'   genepopDat <- read_genepop_format("/home/data/genepop_format_example.gen",digits_per_allele=3)
#'   ## Extract the loci names (that were read in from the top of the file):
#'   locnames <- genepopDat$locnames
#'   ## Separate out the data:
#'   dat <- genepopDat$popData
#'   ## You could then run GenePlot on the data and locnames.
#'   ## Note that by default, data read in from Genepop format will have populations
#'   ## called Pop1, Pop2 etc. unless the individuals in that pop have non-unique
#'   ## names, in which case they will be given the ID of the first individual
#'   ## in that pop as their pop name and will be given auto-generated unique IDs.
#'   geneplot(dat, refpopnames=c("Mahu","Taik"), include.pops=c("Flat"), locnames=genepopDat$locnames)
#' }
#'
#' @references McMillan, L. and Fewster, R. "Visualizations for genetic assignment
#'  analyses using the saddlepoint approximation method" (2017) \emph{Biometrics}.
#'  Rannala, B., and Mountain, J. L. (1997). Detecting immigration by using multilocus
#'  genotypes. \emph{Proceedings of the National Academy of Sciences} \strong{94}, 9197--9201.
#'  Piry, S., Alapetite, A., Cornuet, J.-M., Paetkau, D., Baudouin, L., and
#'  Estoup, A. (2004). GENECLASS2: A software for genetic assignment and
#'  first-generation migrant detection. \emph{Journal of Heredity} \strong{95}, 536--539.
#'
#' @author Log-Genotype-Probability calculations based on the method of Rannala
#' and Mountain (1997) as implemented in GeneClass2, updated to allow for individuals
#' with missing data and to enable accurate calculations of quantiles of the
#' Log-Genotype-Probability distributions of the reference populations.
#' See McMillan and Fewster (2017) for details.
#'
#' @export
geneplot <- function(dat, refpopnames, locnames,
                     includepopnames=NULL,
                     prior="Rannala",
                     saddlepoint=T,
                     leave_one_out=F,
                     logten=T,
                     min_loci=6,
                     quantiles=c(0.01, 1.00),
                     Ndraw=100000,
                     plotit=T,
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
                     cexpts=1.4){

    logprob_results <- calc_logprob(dat=dat, refpopnames=refpopnames, locnames=locnames,
                                    includepopnames=includepopnames, prior=prior,
                                    logten=logten, leave_one_out=leave_one_out, saddlepoint=saddlepoint,
                                    min_loci=min_loci, quantiles=quantiles, Ndraw=Ndraw)

    if(plotit){

        plot_logprob(logprob_results, orderpop=orderpop, axispop=axispop, plot_type=plot_type,
                     shapevec=shapevec, colvec=colvec, xyrange=xyrange, cexpts=cexpts,
                     txt=txt, dim1=dim1, dim2=dim2, mark_impute=mark_impute,
                     use_legend=use_legend, legend_pos=legend_pos, plot_bars=plot_bars,
                     axis_labels=axis_labels, short_axis_labels=short_axis_labels,
                     grayscale_quantiles=grayscale_quantiles,
                     layout_already_set=layout_already_set)
    }

    logprob_results
}
