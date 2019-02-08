context("Testing validation of GenePlot inputs")

test_that("calc_logprob fails for blank data, refpopnames or locnames with all other inputs left as defaults.", {

    expect_that(calc_logprob(NULL,refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat cannot be null."))
    expect_that(calc_logprob(dat=ratData,NULL,locnames=ratLocnames), throws_error("refpopnames cannot be null."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,NULL), throws_error("locnames cannot be null."))
})

test_that("calc_logprob fails for wrong type inputs.", {

    expect_that(calc_logprob(NA,refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob("a",refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob(0,refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob(c(0,1),refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob(as.factor(c(1,2)),refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob(list(a=1,b=2),refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))
    expect_that(calc_logprob(array(1:12,dim=c(2,3,2)),refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat must be a data frame."))

    expect_that(calc_logprob(dat=ratData,NA,locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,0,locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,c(0,1),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,as.factor(c(1,2)),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,list(a=1,b=2),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,array(1:12,dim=c(2,3,2)),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,data.frame(a=c(1,2),b=c(1,2)),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))

    expect_that(calc_logprob(dat=ratData,vector(),locnames=ratLocnames), throws_error("refpopnames must be a character vector of population/sample names."))
    expect_that(calc_logprob(dat=ratData,character(),locnames=ratLocnames), throws_error("refpopnames is an empty vector."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,NA), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,0), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,c(0,1)), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,as.factor(c(1,2))), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,list(a=1,b=2)), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,array(1:12,dim=c(2,3,2))), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,data.frame(a=c(1,2),b=c(1,2))), throws_error("locnames must be a character vector of locus names."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,vector()), throws_error("locnames must be a character vector of locus names."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,character()), throws_error("locnames is an empty vector."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=NA), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=0), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=c(0,1)), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=as.factor(c(1,2))), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=list(a=1,b=2)), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=array(1:12,dim=c(2,3,2))), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=data.frame(a=c(1,2),b=c(1,2))), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=vector()), throws_error("includepopnames must be a character vector of sample names to be compared to the reference populations."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=character()), throws_error("includepopnames is an empty vector."))

})

test_that("calc_logprob fails for invalid or mismatched dat, refpopnames, locnames or includepopnames.", {

    temp1 <- ratData
    temp1$pop <- NULL
    expect_that(calc_logprob(dat=temp1,refpopnames=refpopnames,locnames=ratLocnames), throws_error("dat does not have a column named 'pop' containing the population/sample labels."))

    temp2 <- ratData
    temp2$pop <- as.factor(ratData$pop)
    expect_that(calc_logprob(dat=temp2,refpopnames=refpopnames,locnames=ratLocnames), throws_error("The 'pop' column in dat is a factor. Please convert to character and rerun."))

    temp3 <- ratData
    temp3$pop[which(temp3$pop == "Brok")] <- "Pop1"
    expect_that(calc_logprob(dat=temp3,refpopnames=refpopnames,locnames=ratLocnames), throws_error("At least one refpopname is not found in the dataset. Please check the inputs."))

    expect_that(calc_logprob(dat=ratData,refpopnames=c("Pop1","Pop2"),locnames=ratLocnames), throws_error("At least one refpopname is not found in the dataset. Please check the inputs."))

    temp4 <- ratData
    temp4$pop[which(temp2$pop =="Main")] <- "Main T"
    expect_that(calc_logprob(dat=temp4,refpopnames=c("Brok","Main T"),locnames=ratLocnames), throws_error("refpopnames cannot have spaces in them. Please edit refpopnames and the pop field in the dataset before rerunning."))

    expect_that(calc_logprob(dat=ratData,refpopnames="Brok",locnames=ratLocnames), throws_error("You must specify at least two populations in refpopnames."))

    temp5 <- ratData
    temp5[,10] <- NULL
    expect_that(calc_logprob(dat=temp5,refpopnames=refpopnames,locnames=ratLocnames), throws_error("Some of the locnames are not present in the dataset. Please check that dat has two columns for each locus, named <<locusname>>.a1 and <<locusname>>.a2."))
    temp6 <- ratData
    temp6[,9:10] <- NULL
    expect_that(calc_logprob(dat=temp6,refpopnames=refpopnames,locnames=ratLocnames), throws_error("Some of the locnames are not present in the dataset. Please check that dat has two columns for each locus, named <<locusname>>.a1 and <<locusname>>.a2."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=c(ratLocnames,"LocExtra")), throws_error("Some of the locnames are not present in the dataset. Please check that dat has two columns for each locus, named <<locusname>>.a1 and <<locusname>>.a2."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames,includepopnames = "PopExtra"), throws_error("At least one of includepopnames is not found in the dataset. Please check the inputs."))

    temp7 <- ratData
    temp7[1,] <- c("ID1","Brok",rep(0,times=20))
    expect_that(calc_logprob(dat=temp7,refpopnames=refpopnames,locnames=ratLocnames), throws_error("At least one individual is missing data at all loci. Please remove these individuals before rerunning analysis."))
    temp8 <- ratData
    temp8[24,] <- c("ID1","Kai",rep(0,times=20))
    expect_that(calc_logprob(dat=temp7,refpopnames=refpopnames,locnames=ratLocnames,includepopnames=c("Kai")), throws_error("At least one individual is missing data at all loci. Please remove these individuals before rerunning analysis."))
})


test_that("calc_logprob fails for invalid optional inputs.", {

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=NULL), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=1), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=c(1,2)), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=list(a=1,b=2)), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=matrix(c(1,2),nrow=2)), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=array(1:12,dim=c(2,3,2))), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior=data.frame(a=c(1,2),b=c(1,2))), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, prior="BadPrior"), throws_error("Invalid prior. Options are 'Rannala', 'Baudouin', 'Half', or 'Quarter'."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=NULL), throws_error("'saddlepoint' cannot be null. Please set to TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=0), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint="a"), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=c(1,2)), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=list(a=1,b=2)), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=matrix(c(1,2),nrow=2)), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=array(1:12,dim=c(2,3,2))), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=data.frame(a=c(1,2),b=c(1,2))), throws_error("'saddlepoint' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, saddlepoint=c(TRUE,FALSE)), throws_error("'saddlepoint' must be TRUE or FALSE."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=NULL), throws_error("'leave_one_out' cannot be null. Please set to TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=0), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out="a"), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=c(1,2)), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=list(a=1,b=2)), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=matrix(c(1,2),nrow=2)), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=array(1:12,dim=c(2,3,2))), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=data.frame(a=c(1,2),b=c(1,2))), throws_error("'leave_one_out' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, leave_one_out=c(TRUE,FALSE)), throws_error("'leave_one_out' must be TRUE or FALSE."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=NULL), throws_error("'logten' cannot be null. Please set to TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=0), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten="a"), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=c(1,2)), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=list(a=1,b=2)), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=matrix(c(1,2),nrow=2)), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=array(1:12,dim=c(2,3,2))), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=data.frame(a=c(1,2),b=c(1,2))), throws_error("'logten' must be TRUE or FALSE."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, logten=c(TRUE,FALSE)), throws_error("'logten' must be TRUE or FALSE."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=NULL), throws_error("'min_loci' cannot be null. Please specify a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=0), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci="a"), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=c(1,2)), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=list(a=1,b=2)), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=matrix(c(1,2),nrow=2)), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=array(1:12,dim=c(2,3,2))), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=data.frame(a=c(1,2),b=c(1,2))), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, min_loci=c(TRUE,FALSE)), throws_error("'min_loci' must be a positive integer, smaller than the length of locnames."))

    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=-1), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles="a"), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=c("a","b")), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=c(1,2)), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=list(a=1,b=2)), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=matrix(c(1,2),nrow=2)), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=array(1:12,dim=c(2,3,2))), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=data.frame(a=c(1,2),b=c(1,2))), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))
    expect_that(calc_logprob(dat=ratData,refpopnames=refpopnames,locnames=ratLocnames, quantiles=c(TRUE,FALSE)), throws_error("'quantiles' must be null, or a vector of values between 0 and 1."))

})
