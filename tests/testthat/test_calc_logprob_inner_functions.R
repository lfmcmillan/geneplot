context("Testing inner GenePlot calculations")

test_that("calc_posterior_locus runs the correct calculations.", {

    dat <- data.frame(id=paste0("ID",1:6),pop=rep("Pop1",times=6),
                      "Loc1.a1"=c(1,1,1,1,1,2),"Loc1.a2"=c(1,2,2,2,3,3),
                      stringsAsFactors = FALSE)
    popx2 <- rep("Pop1",times=12)
    refpopnames <- "Pop1"

    expect_equivalent(calc_posterior_locus(loc="Loc1", dat=dat, popx2=popx2, prior="Rannala", refpopnames=refpopnames),
                      as.table(matrix(c(6+1/3,4+1/3,2+1/3),nrow=1)))

    dat <- data.frame(id=paste0("ID",1:8),pop=c(rep("Pop1",times=6),rep("Pop2",times=2)),
                      "Loc1.a1"=c(1,1,1,1,1,2,1,1),"Loc1.a2"=c(1,2,2,2,3,3,1,2),
                      stringsAsFactors = FALSE)
    popx2 <- c(rep("Pop1",times=6),rep("Pop2",times=2),rep("Pop1",times=6),rep("Pop2",times=2))
    refpopnames <- c("Pop1","Pop2")

    expect_equivalent(calc_posterior_locus(loc="Loc1", dat=dat, popx2=popx2, prior="Rannala", refpopnames=refpopnames),
                      as.table(matrix(c(6+1/3,4+1/3,2+1/3,3+1/3,1+1/3,1/3),nrow=2,byrow=TRUE)))

    expect_equivalent(calc_posterior_locus(loc="Loc1", dat=dat, popx2=popx2, prior="Baudouin", refpopnames=refpopnames),
                      as.table(matrix(c(6+1,4+1,2+1,3+1,1+1,1),nrow=2,byrow=TRUE)))

    expect_equivalent(calc_posterior_locus(loc="Loc1", dat=dat, popx2=popx2, prior="Half", refpopnames=refpopnames),
                      as.table(matrix(c(6+1/2,4+1/2,2+1/2,3+1/2,1+1/2,1/2),nrow=2,byrow=TRUE)))

    expect_equivalent(calc_posterior_locus(loc="Loc1", dat=dat, popx2=popx2, prior="Quarter", refpopnames=refpopnames),
                      as.table(matrix(c(6+1/4,4+1/4,2+1/4,3+1/4,1+1/4,1/4),nrow=2,byrow=TRUE)))
})


test_that("calc_logprob_complete runs the correct calculations.", {

    dat_complete <- data.frame(id=paste0("ID",1:8),pop=c(rep("Pop1",times=6),rep("Pop2",times=2)),
                               "Loc1.a1"=c(1,1,1,1,1,2,1,1),"Loc1.a2"=c(1,2,2,2,3,3,1,2),
                               stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    locnames <- "Loc1"
    posterior_nu_list <- list("Loc1"=as.table(matrix(c(6+1/3,4+1/3,2+1/3,3+1/3,1+1/3,1/3),nrow=2,byrow=TRUE,
                                                     dimnames=list(c("Pop1","Pop2"),c("1","2","3")))))
    posterior_nusumvec_list <- lapply(posterior_nu_list,rowSums)

    ## Leave-one-out = FALSE
    Pop1Row1 <- log10((6+1/3)*(6+1/3+1)/((12+1)*(12+1+1)))
    Pop2Row1 <- log10((3+1/3)*(3+1/3+1)/((4+1)*(4+1+1)))

    expect_equivalent(calc_logprob_complete(row=1, dat_complete=dat_complete,
                                            refpopnames=refpopnames,locnames=locnames,
                                            posterior_nu_list=posterior_nu_list,
                                            posterior_nusumvec_list = posterior_nusumvec_list,
                                            leave_one_out=FALSE, logten=TRUE),
                      c(Pop1Row1,Pop2Row1,Pop1Row1,Pop2Row1))

    Pop1Row2 <- log10(2*(6+1/3)*(4+1/3)/((12+1)*(12+1+1)))
    Pop2Row2 <- log10(2*(3+1/3)*(1+1/3)/((4+1)*(4+1+1)))

    expect_equivalent(calc_logprob_complete(row=2, dat_complete=dat_complete,
                                            refpopnames=refpopnames,locnames=locnames,
                                            posterior_nu_list=posterior_nu_list,
                                            posterior_nusumvec_list = posterior_nusumvec_list,
                                            leave_one_out=FALSE, logten=TRUE),
                      c(Pop1Row2,Pop2Row2,Pop1Row2,Pop2Row2))

    ## Leave-one-out = TRUE -- note that Row1 and Row2 individuals are both from
    ## Pop1, so will have LOO LGPs for Pop1, and normal LGPs for Pop2
    ## Leave-one-out homozygote:
    Pop1Row1_LOO <- log10((4+1/3)*(4+1/3+1)/((10+1)*(10+1+1)))
    Pop2Row1_LOO <- log10((3+1/3)*(3+1/3+1)/((4+1)*(4+1+1)))

    expect_equivalent(calc_logprob_complete(row=1, dat_complete=dat_complete,
                                            refpopnames=refpopnames,locnames=locnames,
                                            posterior_nu_list=posterior_nu_list,
                                            posterior_nusumvec_list = posterior_nusumvec_list,
                                            leave_one_out=TRUE, logten=TRUE,
                                            saddlepoint=TRUE),
                      c(Pop1Row1_LOO,Pop2Row1_LOO,Pop1Row1_LOO,Pop2Row1_LOO))

    ## Leave-one-out heterozygote:
    Pop1Row2_LOO <- log10(2*(5+1/3)*(3+1/3)/((10+1)*(10+1+1)))
    Pop2Row2_LOO <- log10(2*(3+1/3)*(1+1/3)/((4+1)*(4+1+1)))

    expect_equivalent(calc_logprob_complete(row=2, dat_complete=dat_complete,
                                            refpopnames=refpopnames,locnames=locnames,
                                            posterior_nu_list=posterior_nu_list,
                                            posterior_nusumvec_list = posterior_nusumvec_list,
                                            leave_one_out=TRUE, logten=TRUE,
                                            saddlepoint=TRUE),
                      c(Pop1Row2_LOO,Pop2Row2_LOO,Pop1Row2_LOO,Pop2Row2_LOO))

})

test_that("calc_logprob_missing runs the correct calculations.", {

    dat_missing <- data.frame(id=paste0("ID",1:8),pop=c(rep("Pop1",times=6),rep("Pop2",times=2)),
                              "Loc1.a1"=c(1,1,1,1,1,2,1,1),"Loc1.a2"=c(1,2,2,2,3,3,1,2),
                              "Loc2.a1"=c(0,0,1,1,1,1,2,2),"Loc2.a2"=c(0,0,1,1,1,1,2,2),
                              stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    locnames <- "Loc1"
    posterior_nu_list <- list("Loc1"=as.table(matrix(c(6+1/3,4+1/3,2+1/3,3+1/3,1+1/3,1/3),nrow=2,byrow=TRUE,
                                                     dimnames=list(c("Pop1","Pop2"),c("1","2","3")))),
                              "Loc2"=as.table(matrix(c(8+1/2,1/2,4+1/2,1/2),nrow=2,byrow=TRUE,
                                                     dimnames=list(c("Pop1","Pop2"),c("1","2")))))
    posterior_nusumvec_list <- lapply(posterior_nu_list,rowSums)

    all_loci_SCDF_qsearch_params <- calc_qsearch_params_usingRcpp(posterior_nu_list,c("Pop1","Pop2"),logten=TRUE,leave_one_out=FALSE)

    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop1$Fh_min > 0 && all_loci_SCDF_qsearch_params$Pop1$Fh_min < 1)
    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop1$Fh_max > 0 && all_loci_SCDF_qsearch_params$Pop1$Fh_max < 1)
    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop1$Fh_max > all_loci_SCDF_qsearch_params$Pop1$Fh_min)
    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop2$Fh_min > 0 && all_loci_SCDF_qsearch_params$Pop2$Fh_min < 1)
    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop2$Fh_max > 0 && all_loci_SCDF_qsearch_params$Pop2$Fh_max < 1)
    expect_equal(TRUE,all_loci_SCDF_qsearch_params$Pop2$Fh_max > all_loci_SCDF_qsearch_params$Pop2$Fh_min)

    ## Leave-one-out = FALSE
    Pop1Row1_raw <- log10((6+1/3)*(6+1/3+1)/((12+1)*(12+1+1)))
    Pop2Row1_raw <- log10((3+1/3)*(3+1/3+1)/((4+1)*(4+1+1)))

    ## Test the details of the Pop1 Loc1 distribution, which for the individual
    ## in row1 is the Pop1 distribution for the loci they have (i.e. Loc1 only)
    Pop1_post_nu_pop <- lapply(posterior_nu_list["Loc1"], function(x) x["Pop1",,drop=FALSE]) ## Note: must use "lapply" and "drop=FALS to preserve correct format
    Pop1_post_info <- rcpp_calc_multi_locus_dist(Pop1_post_nu_pop, leave_one_out=FALSE)
    Pop1_Loc1_dist <- Pop1_post_info$dist[[1]]
    Pop1_mean_pop <- rcpp_calc_mu(Pop1_post_info$dist)

    ## Test the mean of the distribution
    expect_equivalent(sum(Pop1_Loc1_dist$values*(Pop1_Loc1_dist$probs)),Pop1_mean_pop)

    Pop1_indiv_p <- rcpp_calc_Fhat(Pop1Row1_raw, Pop1_post_info$dist,
                                   Pop1_post_info$min, Pop1_post_info$max,
                                   Pop1_mean_pop, logten=TRUE)

    Pop1_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, posterior_nu_list,
                                       "Pop1", Pop1_indiv_p, logten=TRUE, leave_one_out=FALSE)

    Pop2_post_nu_pop <- lapply(posterior_nu_list["Loc1"], function(x) x["Pop2",,drop=FALSE]) ## Note: must use "lapply" to preserve correct format
    Pop2_post_info <- rcpp_calc_multi_locus_dist(Pop2_post_nu_pop, leave_one_out=FALSE)
    Pop2_mean_pop <- rcpp_calc_mu(Pop2_post_info$dist)
    Pop2_indiv_p <- rcpp_calc_Fhat(Pop2Row1_raw, Pop2_post_info$dist,
                                   Pop2_post_info$min, Pop2_post_info$max,
                                   Pop2_mean_pop, logten=TRUE)

    Pop2_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, posterior_nu_list,
                                       "Pop2", Pop2_indiv_p, logten=TRUE, leave_one_out=FALSE)

    expect_equivalent(calc_logprob_missing(row=1, dat_missing=dat_missing,
                                           refpopnames=refpopnames,locnames=locnames,
                                           posterior_nu_list=posterior_nu_list,
                                           posterior_nusumvec_list = posterior_nusumvec_list,
                                           leave_one_out=FALSE, logten=TRUE,
                                           saddlepoint=TRUE,
                                           all_loci_SCDF_qsearch_params = all_loci_SCDF_qsearch_params),
                      c(1,Pop1_q,Pop2_q,Pop1Row1_raw,Pop2Row1_raw))

    ## Leave-one-out = TRUE -- note that Row1 and Row2 individuals are both from
    ## Pop1, so will have LOO LGPs for Pop1, and normal LGPs for Pop2
    ## And individual 1 and 2 only have data for Loc1, so the LOO part only
    ## applies to Loc1
    ## Leave-one-out homozygote:
    Pop1Row1_LOO <- log10((4+1/3)*(4+1/3+1)/((10+1)*(10+1+1)))

    expect_equivalent(calc_logprob_missing(row=1, dat_missing=dat_missing,
                                           refpopnames=refpopnames,locnames=locnames,
                                           posterior_nu_list=posterior_nu_list,
                                           posterior_nusumvec_list = posterior_nusumvec_list,
                                           leave_one_out=TRUE, logten=TRUE,
                                           saddlepoint=TRUE,
                                           all_loci_SCDF_qsearch_params = all_loci_SCDF_qsearch_params)['Pop1.raw'],
                      Pop1Row1_LOO)

    ## Leave-one-out heterozygote:
    Pop1Row2_LOO <- log10(2*(5+1/3)*(3+1/3)/((10+1)*(10+1+1)))

    expect_equivalent(calc_logprob_missing(row=2, dat_missing=dat_missing,
                                           refpopnames=refpopnames,locnames=locnames,
                                           posterior_nu_list=posterior_nu_list,
                                           posterior_nusumvec_list = posterior_nusumvec_list,
                                           leave_one_out=TRUE, logten=TRUE,
                                           saddlepoint=TRUE,
                                           all_loci_SCDF_qsearch_params = all_loci_SCDF_qsearch_params)['Pop1.raw'],
                      Pop1Row2_LOO)

})

test_that("rcpp_calc_dist runs the correct calculations.", {

    nu1 <- matrix(2, nrow=1)
    dist1 <- rcpp_calc_dist(nu1,FALSE)
    expect_equivalent(dist1$values[1],0)
    expect_equivalent(dist1$probs[1],1)

    nu2 <- matrix(c(1.5,1.5), nrow=1)
    dist2 <- rcpp_calc_dist(nu2,FALSE)
    expect_equivalent(dist2$probs,c(1.5*2.5/3/4,1.5*2.5/3/4,2*1.5*1.5/3/4))
    expect_equivalent(log(1.5*2.5/3/4),rcpp_calc_min_beta(nu2))
    expect_equivalent(log(2*1.5*1.5/3/4),rcpp_calc_max_beta(nu2))

    nu3 <- matrix(c(9.5,1.5), nrow=1)
    dist3 <- rcpp_calc_dist(nu3,FALSE)
    expect_equivalent(dist3$probs,c(1.5*2.5/11/12,2*1.5*9.5/11/12,9.5*10.5/11/12))
    expect_equivalent(log(1.5*2.5/11/12),rcpp_calc_min_beta(nu3))
    expect_equivalent(log(9.5*10.5/11/12),rcpp_calc_max_beta(nu3))

    nu4 <- matrix(c(10+1/3,1/3,1/3), nrow=1)
    dist4 <- rcpp_calc_dist(nu4,FALSE)
    expect_equivalent(dist4$probs,c(2*1/3*1/3/11/12,
                                    1/3*(1/3+1)/11/12,
                                    1/3*(1/3+1)/11/12,
                                    2*1/3*(10+1/3)/11/12,
                                    2*1/3*(10+1/3)/11/12,
                                    (10+1/3)*(11+1/3)/11/12))
    expect_equivalent(log(2*1/3*1/3/11/12),rcpp_calc_min_beta(nu4))
    expect_equivalent(log((10+1/3)*(11+1/3)/11/12),rcpp_calc_max_beta(nu4))

})
