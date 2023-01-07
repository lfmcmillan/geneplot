context("Overall geneplot and calc_logprob tests.")

test_that("Default inputs for geneplot and calc_logprob match and are as expected", {

    result_geneplot <- geneplot(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,plotit=FALSE)

    expect_equal(attributes(result_geneplot)$prior,"Rannala")
    expect_equal(attributes(result_geneplot)$saddlepoint,TRUE)
    expect_equal(attributes(result_geneplot)$leave_one_out,TRUE)
    expect_equal(attributes(result_geneplot)$logten,TRUE)
    expect_equal(attributes(result_geneplot)$min_loci,6)

    result_calc <- calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames)

    expect_equal(result_geneplot,result_calc)
    expect_equal(attributes(result_geneplot)$prior,attributes(result_calc)$prior)
    expect_equal(attributes(result_geneplot)$saddlepoint,attributes(result_calc)$saddlepoint)
    expect_equal(attributes(result_geneplot)$leave_one_out,attributes(result_calc)$leave_one_out)
    expect_equal(attributes(result_geneplot)$logten,attributes(result_calc)$logten)
    expect_equal(attributes(result_geneplot)$min_loci,attributes(result_calc)$min_loci)
})
