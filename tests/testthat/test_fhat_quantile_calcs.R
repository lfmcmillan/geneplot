context("Fhat and quantile saddlepoint tests.")

test_that("Fhat and quantile calculations in saddlepoint code match up - dataset 1", {

    dat_missing1 <- data.frame(id=paste0("ID",1:8),pop=c(rep("Pop1",times=6),rep("Pop2",times=2)),
                              "Loc1.a1"=c(1,1,1,1,1,2,1,1),"Loc1.a2"=c(1,2,2,2,3,3,1,2),
                              "Loc2.a1"=c(0,0,1,1,1,1,2,2),"Loc2.a2"=c(0,0,1,1,1,1,2,2),
                              "Loc3.a1"=c(1,1,2,2,1,2,1,1),"Loc3.a2"=c(1,2,2,2,3,3,1,2),
                              "Loc4.a1"=c(1,1,2,1,3,2,1,1),"Loc4.a2"=c(1,2,2,2,3,3,1,2),
                              "Loc5.a1"=c(1,1,1,2,3,2,1,1),"Loc5.a2"=c(1,2,2,2,3,3,1,2),
                              "Loc6.a1"=c(1,2,1,1,1,2,1,1),"Loc6.a2"=c(1,2,2,2,3,3,1,2),
                              stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    data1_posterior_nu_list <- lapply(paste0("Loc",1:6),calc_posterior_locus,dat_missing1,
                                      c(dat_missing1$pop,dat_missing1$pop),
                                      prior="Rannala",refpopnames=refpopnames)
    names(data1_posterior_nu_list) <- paste0("Loc",1:6)
    posterior_nusumvec_list <- lapply(data1_posterior_nu_list,rowSums)

    ## Calculate the LGP for individual ID1 who has missing data at Loc2
    locnames_present <- c("Loc1","Loc3","Loc4","Loc5","Loc6")
    indiv_dat <- dat_missing1[1,]
    indiv_lgp <- 0
    for(loc in locnames_present){
        loc.table <- data1_posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")

        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        if(col.a1==col.a2) {
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        } else {
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote
        }

        indiv_lgp <- indiv_lgp + log10(prob.loc)
    } # end for-locus loop for loci with data intact
    data1_Pop1_lgp <- indiv_lgp['Pop1']

    ## Calculate the Fhat value for individual ID1 based on the loci they have
    data1_Pop1_present_dist_info <- rcpp_calc_multi_locus_dist(data1_posterior_nu_list[locnames_present], leave_one_out=FALSE)
    data1_Pop1_present_mean <- rcpp_calc_mu(data1_Pop1_present_dist_info$dist)

    data1_Pop1_indiv_p <- rcpp_calc_Fhat(data1_Pop1_lgp,
                                         dist=data1_Pop1_present_dist_info$dist,
                                         min_dist=data1_Pop1_present_dist_info$min,
                                         max_dist=data1_Pop1_present_dist_info$max,
                                         mean_dist=data1_Pop1_present_mean, logten=TRUE)

    data1_x <- data1_Pop1_lgp/log10(exp(1))
    data1_x
    data1_sh <- rcpp_calc_shat(data1_x, data1_Pop1_present_dist_info$dist)
    data1_sh
    data1_wh <- rcpp_calc_what(data1_x, data1_sh, data1_Pop1_present_dist_info$dist, use_x=TRUE)
    data1_wh

    data1_x*data1_sh
    rcpp_calc_multi_locus_K(data1_Pop1_present_dist_info$dist, data1_sh)
    rcpp_calc_multi_locus_K1(data1_Pop1_present_dist_info$dist, data1_sh)

    ## Now calculate the full-loci distribution
    all_loci_SCDF_qsearch_params <- calc_qsearch_params_usingRcpp(data1_posterior_nu_list,c("Pop1","Pop2"),
                                                                  logten=TRUE,leave_one_out=FALSE)

    ## Calculate the quantile for individual ID1 in the full-loci distribution
    Pop1_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, data1_posterior_nu_list,
                                       "Pop1", data1_Pop1_indiv_p, logten=TRUE, leave_one_out=FALSE)

    Pop1_all_loci_nu <- lapply(data1_posterior_nu_list, function(x) x["Pop1",,drop=FALSE])
    Pop1_all_loci_info <- rcpp_calc_multi_locus_dist(Pop1_all_loci_nu,leave_one_out=FALSE)
    Pop1_all_loci_mean <- rcpp_calc_mu(Pop1_all_loci_info$dist)

    ## Reverse the quantile calculation to calculate Fhat for individual ID1 in the
    ## full-loci distribution, and check it's roughly the same as the loci-present
    ## Fhat value
    expect_equivalent(data1_Pop1_indiv_p,
                      rcpp_calc_Fhat(Pop1_q, Pop1_all_loci_info$dist, Pop1_all_loci_info$min,
                                     Pop1_all_loci_info$max, Pop1_all_loci_mean, logten=TRUE),
                      tolerance=1E-4)
})


test_that("Fhat and quantile calculations in saddlepoint code match up - dataset 2", {

    dat_missing2 <- data.frame(id=paste0("ID",1:12),pop=c(rep("Pop1",times=8),rep("Pop2",times=4)),
                              "Loc1.a1"=c(1,2,3,4,4,5,6,3,5,2,1,1),"Loc1.a2"=c(2,6,3,2,2,2,3,3,1,2,3,1),
                              "Loc2.a1"=c(0,0,1,1,2,1,2,2,1,1,1,1),"Loc2.a2"=c(0,0,1,1,1,1,2,2,1,1,1,1),
                              "Loc3.a1"=c(1,1,2,2,1,2,1,3,1,5,3,2),"Loc3.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc4.a1"=c(1,1,2,1,3,2,1,1,2,2,2,2),"Loc4.a2"=c(1,2,2,2,3,3,1,2,4,4,1,1),
                              "Loc5.a1"=c(1,1,1,2,3,2,1,1,2,1,1,1),"Loc5.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc6.a1"=c(1,2,1,1,1,2,1,2,2,2,2,2),"Loc6.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    data2_posterior_nu_list <- lapply(paste0("Loc",1:6),calc_posterior_locus,dat_missing2,
                                      c(dat_missing2$pop,dat_missing2$pop),
                                      prior="Rannala",refpopnames=refpopnames)
    names(data2_posterior_nu_list) <- paste0("Loc",1:6)
    posterior_nusumvec_list <- lapply(data2_posterior_nu_list,rowSums)

    ## Calculate the LGP for individual ID1 who has missing data at Loc2
    locnames_present <- c("Loc1","Loc3","Loc4","Loc5","Loc6")
    indiv_dat <- dat_missing2[1,]
    indiv_lgp <- 0
    for(loc in locnames_present){
        loc.table <- data2_posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")

        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        if(col.a1==col.a2) {
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        } else {
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote
        }

        indiv_lgp <- indiv_lgp + log10(prob.loc)
    } # end for-locus loop for loci with data intact
    data2_Pop1_lgp <- indiv_lgp['Pop1']

    ## Calculate the Fhat value for individual ID1 based on the loci they have
    data2_Pop1_present_dist_info <- rcpp_calc_multi_locus_dist(data2_posterior_nu_list[locnames_present], leave_one_out=FALSE)
    data2_Pop1_present_mean <- rcpp_calc_mu(data2_Pop1_present_dist_info$dist)

    data2_Pop1_indiv_p <- rcpp_calc_Fhat(data2_Pop1_lgp,
                                         dist=data2_Pop1_present_dist_info$dist,
                                         data2_Pop1_present_dist_info$min,
                                         data2_Pop1_present_dist_info$max,
                                         mean_dist=data2_Pop1_present_mean, logten=TRUE)

    data2_x <- data2_Pop1_lgp/log10(exp(1))
    data2_x
    data2_sh <- rcpp_calc_shat(data2_x, data2_Pop1_present_dist_info$dist)
    data2_sh
    data2_wh <- rcpp_calc_what(data2_x, data2_sh, data2_Pop1_present_dist_info$dist, use_x=TRUE)
    data2_wh

    data2_x*data2_sh
    rcpp_calc_multi_locus_K(data2_Pop1_present_dist_info$dist, data2_sh)

    ## Now calculate the full-loci distribution
    all_loci_SCDF_qsearch_params <- calc_qsearch_params_usingRcpp(data2_posterior_nu_list,c("Pop1","Pop2"),logten=TRUE,leave_one_out=FALSE)

    ## Calculate the quantile for individual ID1 in the full-loci distribution
    Pop1_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, data2_posterior_nu_list,
                                       "Pop1", data2_Pop1_indiv_p, logten=TRUE, leave_one_out=FALSE)

    Pop1_all_loci_nu <- lapply(data2_posterior_nu_list, function(x) x["Pop1",,drop=FALSE])
    Pop1_all_loci_info <- rcpp_calc_multi_locus_dist(Pop1_all_loci_nu,leave_one_out=FALSE)
    Pop1_all_loci_mean <- rcpp_calc_mu(Pop1_all_loci_info$dist)

    ## Reverse the quantile calculation to calculate Fhat for individual ID1 in the
    ## full-loci distribution, and check it's roughly the same as the loci-present
    ## Fhat value
    expect_equivalent(data2_Pop1_indiv_p,
                      rcpp_calc_Fhat(Pop1_q, Pop1_all_loci_info$dist,
                                     Pop1_all_loci_info$min, Pop1_all_loci_info$max,
                                     mean_dist=Pop1_all_loci_mean, logten=TRUE),
                      tolerance=1E-4)
})

test_that("Fhat and quantile calculations in saddlepoint code match up - dataset 3", {

    dat_missing3 <- data.frame(id=paste0("ID",1:12),pop=c(rep("Pop1",times=8),rep("Pop2",times=4)),
                              "Loc1.a1"=c(1,2,3,4,4,5,6,3,5,2,1,1),"Loc1.a2"=c(2,6,3,2,2,2,3,3,1,2,3,1),
                              "Loc2.a1"=c(0,0,1,1,2,1,2,2,1,1,1,1),"Loc2.a2"=c(0,0,1,1,1,1,2,2,1,1,1,1),
                              "Loc3.a1"=c(2,1,2,2,1,2,1,3,1,5,3,2),"Loc3.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc4.a1"=c(2,1,2,1,3,2,1,1,2,2,2,2),"Loc4.a2"=c(1,2,2,2,3,3,1,2,4,4,1,1),
                              "Loc5.a1"=c(2,1,1,2,3,2,1,1,2,1,1,1),"Loc5.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc6.a1"=c(1,2,1,1,1,2,1,2,2,2,2,2),"Loc6.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    data3_posterior_nu_list <- lapply(paste0("Loc",1:6),calc_posterior_locus,dat_missing3,
                                      c(dat_missing3$pop,dat_missing3$pop),
                                      prior="Rannala",refpopnames=refpopnames)
    names(data3_posterior_nu_list) <- paste0("Loc",1:6)
    posterior_nusumvec_list <- lapply(data3_posterior_nu_list,rowSums)

    ## Calculate the LGP for individual ID1 who has missing data at Loc2
    locnames_present <- c("Loc1","Loc3","Loc4","Loc5","Loc6")
    indiv_dat <- dat_missing3[1,]
    indiv_lgp <- 0
    for(loc in locnames_present){
        loc.table <- data3_posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")

        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        if(col.a1==col.a2) {
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        } else {
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote
        }

        indiv_lgp <- indiv_lgp + log10(prob.loc)
    } # end for-locus loop for loci with data intact
    data3_Pop1_lgp <- indiv_lgp['Pop1']

    ## Calculate the Fhat value for individual ID1 based on the loci they have
    data3_Pop1_present_dist_info <- rcpp_calc_multi_locus_dist(data3_posterior_nu_list[locnames_present], leave_one_out=FALSE)
    data3_Pop1_present_mean <- rcpp_calc_mu(data3_Pop1_present_dist_info$dist)

    data3_Pop1_indiv_p <- rcpp_calc_Fhat(data3_Pop1_lgp,
                                         data3_Pop1_present_dist_info$dist,
                                         data3_Pop1_present_dist_info$min,
                                         data3_Pop1_present_dist_info$max,
                                         data3_Pop1_present_mean, logten=TRUE)

    data3_x <- data3_Pop1_lgp/log10(exp(1))
    data3_x
    data3_sh <- rcpp_calc_shat(data3_x, data3_Pop1_present_dist_info$dist)
    data3_sh
    data3_wh <- rcpp_calc_what(data3_x, data3_sh, data3_Pop1_present_dist_info$dist, use_x=TRUE)
    data3_wh

    data3_x*data3_sh
    rcpp_calc_multi_locus_K(data3_Pop1_present_dist_info$dist, data3_sh)

    ## Now calculate the full-loci distribution
    all_loci_SCDF_qsearch_params <- calc_qsearch_params_usingRcpp(data3_posterior_nu_list,c("Pop1","Pop2"),logten=TRUE,leave_one_out=FALSE)

    ## Calculate the quantile for individual ID1 in the full-loci distribution
    Pop1_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, data3_posterior_nu_list,
                                       "Pop1", data3_Pop1_indiv_p, logten=TRUE, leave_one_out=FALSE)

    Pop1_all_loci_nu <- lapply(data3_posterior_nu_list, function(x) x["Pop1",,drop=FALSE])
    Pop1_all_loci_info <- rcpp_calc_multi_locus_dist(Pop1_all_loci_nu,leave_one_out=FALSE)
    Pop1_all_loci_mean <- rcpp_calc_mu(Pop1_all_loci_info$dist)

    ## Reverse the quantile calculation to calculate Fhat for individual ID1 in the
    ## full-loci distribution, and check it's roughly the same as the loci-present
    ## Fhat value
    expect_equivalent(data3_Pop1_indiv_p,
                      rcpp_calc_Fhat(Pop1_q, Pop1_all_loci_info$dist,
                                     Pop1_all_loci_info$min, Pop1_all_loci_info$max,
                                     Pop1_all_loci_mean, logten=TRUE),
                      tolerance=1E-4)
})


test_that("Fhat and quantile calculations in saddlepoint code match up - dataset 4", {

    # Note that this is only different to the 3rd dataset in that the alleles
    # of the first 2 observations are in the opposite order, but that should not
    # change the overall results
    dat_missing4 <- data.frame(id=paste0("ID",1:12),pop=c(rep("Pop1",times=8),rep("Pop2",times=4)),
                              "Loc1.a1"=c(1,2,3,4,4,5,6,3,5,2,1,1),"Loc1.a2"=c(2,6,3,2,2,2,3,3,1,2,3,1),
                              "Loc2.a1"=c(0,0,1,1,2,1,2,2,1,1,1,1),"Loc2.a2"=c(0,0,1,1,1,1,2,2,1,1,1,1),
                              "Loc3.a1"=c(1,2,2,2,1,2,1,3,1,5,3,2),"Loc3.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc4.a1"=c(1,2,2,1,3,2,1,1,2,2,2,2),"Loc4.a2"=c(1,2,2,2,3,3,1,2,4,4,1,1),
                              "Loc5.a1"=c(1,2,1,2,3,2,1,1,2,1,1,1),"Loc5.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              "Loc6.a1"=c(1,2,1,1,1,2,1,2,2,2,2,2),"Loc6.a2"=c(1,2,2,2,3,3,1,2,1,1,1,1),
                              stringsAsFactors = FALSE)
    refpopnames <- c("Pop1","Pop2")
    data4_posterior_nu_list <- lapply(paste0("Loc",1:6),calc_posterior_locus,dat_missing4,
                                      c(dat_missing4$pop,dat_missing4$pop),
                                      prior="Rannala",refpopnames=refpopnames)
    names(data4_posterior_nu_list) <- paste0("Loc",1:6)
    posterior_nusumvec_list <- lapply(data4_posterior_nu_list,rowSums)

    ## Calculate the LGP for individual ID1 who has missing data at Loc2
    locnames_present <- c("Loc1","Loc3","Loc4","Loc5","Loc6")
    indiv_dat <- dat_missing4[1,]
    indiv_lgp <- 0
    for(loc in locnames_present){
        loc.table <- data4_posterior_nu_list[[loc]]
        nusumvec <- posterior_nusumvec_list[[loc]]
        a1.name <- paste0(loc, ".a1")
        a2.name <- paste0(loc, ".a2")

        indiv_a1 <- as.character(indiv_dat[a1.name])
        indiv_a2 <- as.character(indiv_dat[a2.name])

        col.a1 <- which(colnames(loc.table)==indiv_a1)
        col.a2 <- which(colnames(loc.table)==indiv_a2)

        if(col.a1==col.a2) {
            prob.loc <- loc.table[,col.a1] * (loc.table[,col.a1] + 1) / ( nusumvec * (nusumvec+1) ) # Homozygote
        } else {
            prob.loc <- 2 * loc.table[,col.a1] * loc.table[,col.a2]  / ( nusumvec * (nusumvec+1) ) # Heterozygote
        }

        indiv_lgp <- indiv_lgp + log10(prob.loc)
    } # end for-locus loop for loci with data intact
    data4_Pop1_lgp <- indiv_lgp['Pop1']

    ## Calculate the Fhat value for individual ID1 based on the loci they have
    data4_Pop1_present_dist_info <- rcpp_calc_multi_locus_dist(data4_posterior_nu_list[locnames_present], leave_one_out=FALSE)
    data4_Pop1_present_mean <- rcpp_calc_mu(data4_Pop1_present_dist_info$dist)

    data4_Pop1_indiv_p <- rcpp_calc_Fhat(data4_Pop1_lgp,
                                         data4_Pop1_present_dist_info$dist,
                                         data4_Pop1_present_dist_info$min,
                                         data4_Pop1_present_dist_info$max,
                                         data4_Pop1_present_mean, logten=TRUE)

    data4_x <- data4_Pop1_lgp/log10(exp(1))
    data4_x
    data4_sh <- rcpp_calc_shat(data4_x, data4_Pop1_present_dist_info$dist)
    data4_sh
    data4_wh <- rcpp_calc_what(data4_x, data4_sh, data4_Pop1_present_dist_info$dist, use_x=TRUE)
    data4_wh

    data4_x*data4_sh
    rcpp_calc_multi_locus_K(data4_Pop1_present_dist_info$dist, data4_sh)

    ## Now calculate the full-loci distribution
    all_loci_SCDF_qsearch_params <- calc_qsearch_params_usingRcpp(data4_posterior_nu_list,c("Pop1","Pop2"),logten=TRUE,leave_one_out=FALSE)

    ## Calculate the quantile for individual ID1 in the full-loci distribution
    Pop1_q <- calc_quantiles_usingRcpp(all_loci_SCDF_qsearch_params, data4_posterior_nu_list,
                                     "Pop1", data4_Pop1_indiv_p, logten=TRUE, leave_one_out=FALSE)

    Pop1_all_loci_nu <- lapply(data4_posterior_nu_list, function(x) x["Pop1",,drop=FALSE])
    Pop1_all_loci_info <- rcpp_calc_multi_locus_dist(Pop1_all_loci_nu,leave_one_out=FALSE)
    Pop1_all_loci_mean <- rcpp_calc_mu(Pop1_all_loci_info$dist)

    ## Reverse the quantile calculation to calculate Fhat for individual ID1 in the
    ## full-loci distribution, and check it's roughly the same as the loci-present
    ## Fhat value
    expect_equivalent(data4_Pop1_indiv_p,
                      rcpp_calc_Fhat(Pop1_q, Pop1_all_loci_info$dist,
                                     Pop1_all_loci_info$min, Pop1_all_loci_info$max,
                                     Pop1_all_loci_mean, logten=TRUE),
                      tolerance=1E-4)
})
