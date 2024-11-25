## Function for calculating the DCM distribution at a single locus
calc.probs.func = function(nu, leave_one_out=FALSE, differenceLGPs=FALSE)
{
    calc.probs <- function(nu, leave_one_out=FALSE)
    {
        n = sum(nu)
        k = length(nu)

        if (!leave_one_out) {

            ## Non-leave-one-out version
            genotype.func <- function(nu1,nu2) nu1*nu2/(n*(n+1))

            probs.matrix <- outer(nu,nu,genotype.func)

            probs.diag <- diag(probs.matrix) + nu/(n*(n+1))
            probs.upper <- as.vector(probs.matrix[upper.tri(probs.matrix)])*2

        } else {
            ## For LOO version of SPDF saddlepoint, we want to remove copies of the
            ## alleles from the observed + prior before calculating the genotype
            ## probabilities. We remove as many copies as we can, so if there aren't
            ## any copies in the observed data then we don't remove any for calculating
            ## the LOO genotype probability.

            ## Note that nu[idx] == 1 would mean that there were no alleles of that type
            ## observed, and nu[idx] == 2 would mean fewer than 2 alleles of that type
            ## observed, because the prior will always be bigger than 0 so prior + nobserved
            ## will always be bigger than 1 if nobserved is at least 1. That's why we
            ## only remove 2 copies if nu > 2 and only remove 1 if nu > 1
            genotype.LOO.func <- function(idx1,idx2) {
                if (idx1 == idx2) {
                    if (nu[idx1] > 2) genotype.prob <- (nu[idx1]-2)*(nu[idx1]-1)/((n-2)*(n-1))
                    else if (nu[idx1] > 1) genotype.prob <- (nu[idx1]-1)*nu[idx1]/((n-1)*n)
                    else genotype.prob <- (nu[idx1])*(nu[idx1]+1)/(n*(n+1))
                } else {
                    if (nu[idx1] > 1 & nu[idx2] > 1) genotype.prob <- 2*(nu[idx1]-1)*(nu[idx2]-1)/((n-2)*(n-1))
                    else if (nu[idx1] > 1 & nu[idx2] <= 1) genotype.prob <- 2*(nu[idx1]-1)*nu[idx2]/((n-1)*n)
                    else if (nu[idx1] <= 1 & nu[idx2] > 1) genotype.prob <- 2*nu[idx1]*(nu[idx2]-1)/((n-1)*n)
                    else genotype.prob <- 2*nu[idx1]*nu[idx2]/(n*(n+1))
                }
                genotype.prob
            }

            probs.matrix <- outer(1:k,1:k,Vectorize(genotype.LOO.func))

            probs.diag <- diag(probs.matrix)
            probs.upper <- as.vector(probs.matrix[upper.tri(probs.matrix)])
        }

        probs <- c(probs.diag, probs.upper)
    }

    ## If nu has 2 rows, then treat each row separately - the first is
    ## the base population nu values - for the probabilities, the second
    ## is the comparison population nu values
    if (!is.null(dim(nu)))
    {
        if (differenceLGPs) {
            nu1 <- nu[1,] ## Base population
            nu2 <- nu[2,] ## Comparison population, we want to know how well A individuals fit into B

            ## Calculate the probabilities for populations A and B, then remove
            ## any that are zero before calculating logs. Do this BEFORE doing
            ## any sorting by size, so that the corresponding entries in the two
            ## populations are still aligned with each other

            probsA <- calc.probs(nu1)
            probsA.LOO <- calc.probs(nu1,leave_one_out=leave_one_out)
            probsB <- calc.probs(nu2)

            if (any(probsA.LOO==0|probsB == 0))
            {
                zeroIdxs <- which(probsA.LOO == 0 | probsB == 0)
                probsA <- probsA[-zeroIdxs]
                probsA.LOO <- probsA.LOO[-zeroIdxs]
                probsB <- probsB[-zeroIdxs]
            }
            values <- log(probsA.LOO) - log(probsB)
            probs <- probsA ## The probs are the NON-LOO probs for popA
        } else {
            nu1 <- nu[1,] ## Base population
            nu2 <- nu[2,] ## Comparison population, we want to know how well B fits into A

            probs <- calc.probs(nu2)

            ## Calculate the probabilities for population 1, then remove any that are
            ## zero before calculating logs. Also remove the corresponding entries
            ## from the pop2 probs -- do this BEFORE doing any sorting by size, so that
            ## the corresponding entries in the two populations are still aligned
            ## with each other
            probs1 <- calc.probs(nu1)
            if (any(probs1 == 0))
            {
                probs <- probs[-which(probs1 == 0)]
                probs1 <- probs1[-which(probs1 == 0)]
            }
            values <- log(probs1)
        }

        ## sort the probabilities in increasing order so that when the table
        ## of frequencies is created, it will be in the same order as unique(values)
        sortResult <- sort.int(values, decreasing=FALSE, index.return=TRUE)
        probs <- probs[sortResult$ix]
        values <- sortResult$x

        names(probs) <- NULL
        names(values) <- NULL
    }
    else
    {
        probs <- calc.probs(nu,leave_one_out=FALSE)

        if (leave_one_out) {
            probsForValues <- calc.probs(nu,leave_one_out=TRUE)

            ## remove any zero probabilities
            if (any(probsForValues == 0)) {
                probs <- probs[-which(probsForValues == 0)]
                probsForValues <- probsForValues[-which(probsForValues == 0)]
            }

            ## sort the probabilities in increasing order so that when the table
            ## of frequencies is created, it will be in the same order as unique(probs)
            sortResult <- sort.int(probs, decreasing=FALSE, index.return=TRUE)
            probs <- probs[sortResult$ix]
            probsForValues <- probsForValues[sortResult$ix]
            names(probs) <- NULL
            names(probsForValues) <- NULL

            values <- log(probsForValues)
        } else {
            ## remove any zero probabilities
            if (any(probs == 0)) {
                probs <- probs[-which(probs == 0)]
            }

            ## sort the probabilities in increasing order so that when the table
            ## of frequencies is created, it will be in the same order as unique(probs)
            probs <- sort(probs, decreasing=FALSE)
            names(probs) <- NULL

            values <- log(probs)
        }
    }

    rbind(values, probs)
}

## Function for calculating the minimum value of the DCM distribution at a single locus
calc.min.beta = function(nu)
{
    # If nu values are provided for multiple populations, assume that the
    # second one determines the maximum possible value i.e. reduce nu to
    # just that population
    if (!is.null(dim(nu))) nu = nu[1,]

    n <- sum(nu)
    if (length(nu) > 1)
    {
        # if the least common allele type has a frequency that is less than
        # twice the frequency of the next least common allele minus 1, then the
        # minimum log-prob will be the log-prob of the homozygote of the least
        # common allele, otherwise it will be the log-prob of the heterozygote
        # of the two least common alleles
        # (note that the latter case is only possible if neither of the least
        # common alleles is found in the sample so their posterior frequencies
        # are just the prior frequencies and thus are both less than one)
        min_nu <- min(nu)
        next_min_nu <- min(nu[-which.min(nu)])
        if (min_nu + 1 < 2*next_min_nu) return(log(min_nu*(min_nu+1)/(n*(n+1))))
        else return(log(2*min_nu*next_min_nu/(n*(n+1))))
    }
    else return(0)
}

## Function for calculating the maximum value of the DCM distribution at a single locus
calc.max.beta = function(nu)
{
    # If nu values are provided for multiple populations, assume that the
    # second one determines the maximum possible value i.e. reduce nu to
    # just that population
    if (!is.null(dim(nu))) nu = nu[1,]

    n <- sum(nu)
    if (length(nu) > 1)
    {
        # if the most common allele type has a frequency that is more than
        # twice the frequency of the next most common allele minus one, then the
        # maximum log-prob will be the log-prob of the homozygote of the most
        # common allele otherwise it will be the log-prob of the heterozygote of
        # the two most common alleles
        max_nu <- max(nu)
        next_max_nu <- max(nu[-which.max(nu)])
        if (max_nu + 1 > 2*next_max_nu) return(log(max_nu*(max_nu+1)/(n*(n+1))))
        else return(log(2*max_nu*next_max_nu/(n*(n+1))))
    }
    else return(0)
}

## Function for combining the DCM distributions from multiple loci
calc.multi.locus.probs.func = function(nutab, leave_one_out=FALSE, differenceLGPs=FALSE)
{
    multi.dist = lapply(nutab,calc.probs.func, leave_one_out=leave_one_out, differenceLGPs=differenceLGPs)

    ## when the values of the distribution are the differences between the
    ## LGPs with respect to the two populations, or leave_one_out is being used,
    ## then need to take the min and max from the calculated values instead of usual min/max calc
    if (differenceLGPs || leave_one_out) {
        min.dist = sum(sapply(multi.dist, function(dist) min(dist[1,])))
        max.dist = sum(sapply(multi.dist, function(dist) max(dist[1,])))
    } else {
        ## can simply add the locus values because they are logs of probabilities
        min.dist = sum(sapply(nutab, calc.min.beta))
        max.dist = sum(sapply(nutab, calc.max.beta))
    }

    list(dist=multi.dist, min=min.dist, max=max.dist)
}

## Function for calculating the MGF of a single-locus DCM distribution
M = function(s, probs, r=0, values, rpower.values)
{
    output = 0

    # ## r gives the level of the derivative i.e. M(0,...) is M(s), M(1,...) is Mprime(s)
    # if (is.null(rpower.values)){
    #     browser()
    #     rpower.values = values^r
    # }

    output <- sum(rpower.values*exp(s*values)*probs)

    output
}

## Functions for calculating the CGF and its derivatives for a single locus
single.K <- function(Kparams, s=0) log(M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r0power.values))

single.K1 <- function(Kparams, s=0, twoPops=FALSE)
{
    # if (twoPops) M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r1power.values)/
    #                 M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r0power.values)
    # else sum(Kparams$values*(exp((s+1)*Kparams$values)))/sum(exp((s+1)*Kparams$values))
    if (twoPops)
    {
        sum(Kparams$r1power.values*exp(s*Kparams$values)*Kparams$probs)/
            sum(Kparams$r0power.values*exp(s*Kparams$values)*Kparams$probs)
    } else sum(Kparams$values*(exp((s+1)*Kparams$values)))/sum(exp((s+1)*Kparams$values))
}

single.K2 <- function(Kparams, s=0, twoPops=FALSE)
{
    if (twoPops) M0 <- M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r0power.values)
    else M0 <- sum(exp((s+1)*Kparams$values))

    K2 = M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r2power.values)/M0 -
        (M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r1power.values)/M0)^2

    if (!is.na(K2) && K2 < 0) K2 <- 0
    K2
}

single.K3 = function(Kparams, s=0)
{
    M0 = M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r0power.values)
    M1 = M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r1power.values)
    M2 = M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r2power.values)
    M3 = M(s, Kparams$probs, values=Kparams$values, rpower.values=Kparams$r3power.values)

    M3/M0 - 3*M1*M2/(M0^2) + 2*(M1/M0)^3
}

## Functions for calculating the CGF and its derivatives for multiple loci
multi.K = function(multi.K.params, s=0) sum(sapply(multi.K.params, single.K, s=s))

multi.K1 = function(multi.K.params, s=0, twoPops=FALSE) sum(sapply(multi.K.params, single.K1, s=s, twoPops=twoPops))

multi.K2 = function(multi.K.params, s=0, twoPops=FALSE) sum(sapply(multi.K.params, single.K2, s=s, twoPops=twoPops))

multi.K3 = function(multi.K.params, s=0) sum(sapply(multi.K.params, single.K3, s=s))

## Function for obtaining interim objects for calculating the multi-locus CGF
make.K.params = function(multi.dist)
{
    multi.K.params = list()
    for (i in 1:length(multi.dist))
    {
        values = multi.dist[[i]][1,]
        probs = multi.dist[[i]][2,]

        r0power.values = rep(1,times=length(values))
        r1power.values = values
        r2power.values = values^2
        r3power.values = values^3

        multi.K.params[[i]] = list(probs=probs, values=values, r0power.values=r0power.values, r1power.values=r1power.values,
                                   r2power.values=r2power.values, r3power.values=r3power.values)
    }

    multi.K.params
}

## Function for calculating the mean of the multi-locus LGP distribution
mu <- function(multi.K.params, twoPops=FALSE) multi.K1(multi.K.params,0, twoPops = twoPops)

## Function for calculating the standard deviation of the multi-locus LGP distribution
sigma <- function(multi.K.params, twoPops=FALSE) sqrt(multi.K2(multi.K.params,0, twoPops = twoPops))

## Function for solving K'(shat) = y
shat <- function(x, multi.K.params, dist.info=NULL, mean.pop=NULL,
                 tol=.Machine$double.eps^0.25, maxiter=1000, twoPops=FALSE, verbose=FALSE)
{
    # if (twoPops) f <- function(z) multi.K1(multi.K.params, s=z, twoPops=TRUE) - x
    if (twoPops) f <- function(z) sum(sapply(multi.K.params, single.K1, s=z, twoPops=TRUE)) - x
    else
    {
        K1 <- function(Kparams,s) sum(Kparams$values*(exp((s+1)*Kparams$values)))/sum(exp((s+1)*Kparams$values))
        f <- function(z) sum(sapply(multi.K.params, K1, s=z)) - x
    }

    lower <- -50
    upper <- 150

    f.lower <- multi.K1(multi.K.params, s=lower, twoPops=twoPops) - x
    f.upper <- multi.K1(multi.K.params, s=upper, twoPops=twoPops) - x

    if (is.na(f.lower))
    {
        while(lower < upper - 10)
        {
            lower <- lower + 10
            f.lower <- multi.K1(multi.K.params, s=lower, twoPops=twoPops) - x
        }
    }

    if (is.na(f.upper))
    {
        while(upper > lower + 10)
        {
            upper <- upper - 10
            f.upper <- multi.K1(multi.K.params, s=upper, twoPops=twoPops) - x
        }
    }

    if (is.na(f.lower) || is.na(f.upper)) return(NaN)

    result <- tryCatch({

        solution <- uniroot(f, c(lower,upper), f.lower=f.lower, f.upper=f.upper,tol=tol, extendInt="yes", maxiter=maxiter)
        return(solution$root)

    }, error = function(err) {

        # error handler picks up where error was generated
        message(paste("My error:  ",err))
        if (verbose) {
            message(paste("x is ",x))
            message(paste("Dist Info max is ",dist.info$max))
            message(paste("Dist Info min is ",dist.info$min))
            message(paste("Max is ",sum(sapply(multi.K.params, m <- function(K) max(K$values)))))
            message(paste("Min is ",sum(sapply(multi.K.params, m <- function(K) min(K$values)))))
            message(paste("Initial f.upper is ",f.upper))
            message(paste("Initial f.lower is ",f.lower))
        }

        return(NaN)

    }) # END tryCatch

}

## Function for calculating the saddlepoint approximation to the CDF
Fhat <- function(x, dist.info, mean.pop, logten=FALSE, use.x.in.wh=FALSE, tol=.Machine$double.eps^0.5, twoPops=FALSE)
{
    multi.K.params <- make.K.params(dist.info$dist)
    nloci <- length(multi.K.params)

    pt.Fhat = function(pt.x)
    {
        if (pt.x >= dist.info$max)
        {
            pt.Fh <- 1
        }
        else if (pt.x <= dist.info$min)
        {
            pt.Fh <- 0
        }
        else
        {
            # Within the vicinity of mu=mean.pop, Fhat uses a different formula.
            # If this second formula is used only for pt.x == mu then Fhat
            # is smooth, but in fact due to the difficulty of accurately
            # calculating s_hat, w_hat and u_hat accurately in the vicinity
            # of mu, it is necessary to use the Fhat(mu) formula also for
            # pt.x near mu, to avoid major discontinuities in Fhat. Although
            # Fhat(mu) is not precisely accurate for x != mu, it is far
            # more accurate than the other formula for Fhat when applied to
            # pt.x close to mu. How close to mu you can get before the
            # numerical discontinuity arises depends on the number of loci
            # if (pt.x != mu)
            if (abs(pt.x - mean.pop) > nloci*1E-5)
            {
                sh = shat(pt.x,multi.K.params, mean.pop=mean.pop, tol=tol, twoPops=twoPops)

                if (is.nan(sh))
                {
                    return(NaN)
                }

                wh = what(pt.x,sh,multi.K.params, use.x=use.x.in.wh, twoPops=twoPops)
                pt.Fh <- pnorm(wh) + dnorm(wh)*(1/wh - 1/uhat(sh, multi.K.params, twoPops=twoPops))
            }
            else
            {
                pt.Fh <- 0.5 + multi.K3(multi.K.params,0)/(6*sqrt(2*pi)*(multi.K2(multi.K.params,0, twoPops=twoPops)^1.5))
            }

            ## Correct Fh if it is smaller than 0 or bigger than 1
            if (!is.na(pt.Fh) && pt.Fh < 0)
            {
                pt.Fh <- 0
            }
            else if (!is.na(pt.Fh) && pt.Fh > 1)
            {
                pt.Fh <- 1
            }
        }

        pt.Fh
    }

    if (logten) x = x/log10(exp(1))

    Fh <- sapply(x,pt.Fhat)

    Fh
}

what <- function(x,sh,multi.K.params, use.x = FALSE, twoPops=FALSE)
{
    if (use.x)  wh.inner <- sh*x - multi.K(multi.K.params, sh)
    else        wh.inner <- sh*multi.K1(multi.K.params, sh, twoPops=twoPops) - multi.K(multi.K.params, sh)

    if (wh.inner < 0)
    {
        wh.inner <- 1e-16
    }

    wh <- sign(sh)*sqrt(2*(wh.inner))

    wh
}

uhat <- function(sh,multi.K.params,twoPops=FALSE) sh*sqrt(multi.K2(multi.K.params,sh,twoPops=twoPops))

## Function for calculating the saddlepoint approximation to the PDF
fhat <- function(x, dist.info, logten=FALSE, twoPops=FALSE)
{
    multi.K.params <- make.K.params(dist.info$dist)

    pt.fhat <- function(x, maxiter=1000)
    {
        if (x >= dist.info$max || x <= dist.info$min)
        {
            pt.fh <- 0
        }
        else
        {
            mean.pop <- mu(multi.K.params, twoPops=twoPops)
            if (x == mean.pop) sh <- 0
            else {
                sh <- shat(x, multi.K.params, dist.info=dist.info, mean.pop=mean.pop,
                           maxiter=maxiter,twoPops=twoPops)
                K2 <- multi.K2(multi.K.params, sh, twoPops=twoPops)
                if (is.na(K2)) pt.fh <- NaN
                else if (K2 == 0) pt.fh <- 0
                else pt.fh <- exp(multi.K(multi.K.params,sh) - sh*x)/sqrt(2*pi*K2)
                if (is.infinite(pt.fh)) pt.fh <- exp(multi.K(multi.K.params,sh) - sh*multi.K1(multi.K.params, sh, twoPops=twoPops))/sqrt(2*pi*K2)
            }
        }
        pt.fh
    }

    if (logten) x = x/log10(exp(1))

    fh <- sapply(x, pt.fhat)

    if (any(is.infinite(fh))) stop("Argh still infinite with K1(sh) correction!")
    if (any(is.na(fh)))
    {
        if (all(is.na(fh))) stop("All fhat values are NA or NaN.")
        # if (all(is.na(fh))) browser()
        # else print("NA or NaN fhat values, trying maxiter=1E6.")
        # NA_idxs <- which(is.na(fh))
        # fh[NA_idxs] <- sapply(x[NA_idxs], pt.fhat, maxiter=1E6)
        # if (any(is.na(fh)))
        # {
        #     print("NA or NaN fhat values for maxiter=1E6.")
        # }
    }

    fh
}

Qhat.uniroot <- function(y, Fhat.qsearch.params, nutab, logten=FALSE, twoPops=FALSE)
{
    if (is.nan(y) || is.nan(Fhat.qsearch.params$Fh.max)) browser()

    dist.info <- Fhat.qsearch.params$dist.info

    if (y >= Fhat.qsearch.params$Fh.max)
    {
        Qh <- dist.info$max
        if (logten) Qh <- Qh/log(10)
    }
    else if (y <= Fhat.qsearch.params$Fh.min) ### (y <= 0) was the old version, changed on 2017-02-07
    {
        Qh <- dist.info$min
        if (logten) Qh <- Qh/log(10)
    }
    else
    {
        multi.K.params <- make.K.params(dist.info$dist)
        mean.pop <- mu(multi.K.params, twoPops=twoPops)

        f <- function(x) Fhat(x, dist.info, mean.pop, logten=FALSE, twoPops=twoPops) - y

        if (is.na(Fhat.qsearch.params$Fh.max) || is.na(Fhat.qsearch.params$Fh.min))
        {
            Qh <- NA
        }
        else if (sign(Fhat.qsearch.params$Fh.max-y) == sign(Fhat.qsearch.params$Fh.min-y))
        {
            Qh = NaN
            browser()
        }
        else
        {
            tryCatch({
                solution <- uniroot(f,
                                    lower=Fhat.qsearch.params$x.min,
                                    upper=Fhat.qsearch.params$x.max,
                                    f.upper=Fhat.qsearch.params$Fh.max-y,
                                    f.lower=Fhat.qsearch.params$Fh.min-y,
                                    extendInt="no",
                                    tol = 1e-4)
                Qh <- solution$root

                if (is.na(Qh)) browser()

                if (Qh > dist.info$max) browser() ## Need to check this BEFORE converting via logten!

                if (logten) Qh <- Qh/log(10)
            }, error_fun = function(err) {
                message(paste("y is ",y))
                message(paste("My error:  ",err))
                Qh <- NaN
                browser()
                return(NaN)

            }) # END tryCatch
        }
    }

    Qh
}

## Function to calculate the upper limit of x values to use when searching for quantiles,
## determined by where the calculation of shat starts becoming less accurate,
## and the lower limit of x values, determined when K2 becomes NA
# get.qsearch.params <- function(dist.info, multi.K.params, Fhat.delta.params, s.tol=1e-4)
get.qsearch.params <- function(dist.info, multi.K.params, s.tol=1e-4, logten=FALSE, max.quantile=0.995, twoPops=FALSE)
{
    mean.pop <- mu(multi.K.params)

    ## Find the upper limit for x by taking diff away from the distribution max
    diff.max <- 0.1
    found.max.shat.limit <- FALSE

    Fh_test <- Fhat(dist.info$max-diff.max, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)

    ## If already hitting NaNs for Fhat, try for bigger distance away from xmax.
    while(is.na(Fh_test) && diff.max < (dist.info$max - dist.info$min - 0.1))
    {
        diff.max <- diff.max + 0.1
        Fh_test <- Fhat(dist.info$max-diff.max, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)
    }

    while(!found.max.shat.limit && diff.max > 1e-12 && !is.na(Fh_test))
    {
        x_test <- multi.K1(multi.K.params, shat(dist.info$max-diff.max, multi.K.params, mean.pop=mean.pop, twoPops=twoPops), twoPops=twoPops)
        if (abs(x_test - dist.info$max + diff.max) > s.tol || Fh_test - max.quantile <= 0) found.max.shat.limit = TRUE
        else diff.max <- diff.max/10
        Fh_test <- Fhat(dist.info$max-diff.max, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)
    }
    ## go back to the previous diff value, either because for this one Fhat is NA,
    ## or because we've found the point where x and K1(s(x)) diverge, so we want to go back a step
    diff.max <- diff.max*10

    ## Find the lower limit by taking diffs away from the distribution min
    # diff.min <- 0
    # Changed on 2017-02-07 because having diff.min = 0 is pointless since
    # Fhat(dist.info$min+diff.min) = Fhat(dist.info$min) which is automatically
    # set to 0, so we never look for higher diff.min.
    diff.min <- 0.1

    Fh_test <- Fhat(dist.info$min+diff.min, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)

    ## If already hitting NaNs for Fhat, try for bigger distance away from xmax.
    while((is.na(Fh_test) | Fh_test == 1) && diff.min < (dist.info$max - diff.max - dist.info$min - 0.01))
    {
        diff.min <- diff.min + 0.01
        Fh_test <- Fhat(dist.info$min+diff.min, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)
    }

    if (diff.min > 0.1) print(paste0("diff.min is ",diff.min))

    x.max=dist.info$max-diff.max
    x.min=dist.info$min+diff.min

    # Need to calculate the upper and lower limits of Fhat for logten=FALSE,
    # because the xvalues tested have been in the space of log, not log10
    Fh.max = Fhat(x.max, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)
    Fh.min <- Fhat(x.min, dist.info, mean.pop, logten=FALSE, twoPops=twoPops)

    if (is.nan(Fh.max)) browser()

    if (x.max > dist.info$max) browser()

    list(x.max=x.max, x.min=x.min, Fh.max=Fh.max, Fh.min=Fh.min,
         dist.info=dist.info, logten=logten, max.quantile=max.quantile)
}

## Function to find all possible x values for the given LGP distribution
get.possible.x <- function(dist)
{
    x.so.far <- dist[[1]]['values',]
    for(i in 2:length(dist))
    {
        reduced <- combine.duplicates(dist[[i]])
        x.so.far <- as.vector(outer(x.so.far, reduced['values',], "+"))
    }
    x.so.far
}

## Function to combine duplicate values for the LGP distribution at a single locus
combine.duplicates <- function(loc.dist)
{
    values <- loc.dist['values',]
    probs <- loc.dist['probs',]

    if (anyDuplicated(values) > 0)
    {
        ## need to convert the values vector to character first, because that's what table does,
        ## otherwise the duplicated check uses different level of precision than table,
        ## so entries may be marked as not duplicates but be compressed into the same entry in table
        ## See nutab[[3]] for nutab = gen.many.loci(L=10,n=100,k=6) with set.seed(37)
        values <- sort(values)
        probs <- probs[which(!duplicated(as.character(values)))]*as.vector(table(values))
        values <- unique(values)

        loc.dist <- rbind(values, probs)
    }

    loc.dist
}

calc.qsearch.params <- function(posterior.nu.list, refpopnames, logten=FALSE, leave_one_out=FALSE)
{
    SCDF.qsearch.params <- list()

    for (pop in refpopnames)
    {
        ## Note: must use "lapply" specifically, in order to preserve the right format!!
        post.nu.pop <- lapply(posterior.nu.list, function(x) x[pop,])

        post.info <- calc.multi.locus.probs.func(post.nu.pop, leave_one_out=leave_one_out)
        multi.K.params <- make.K.params(post.info$dist)

        SCDF.qsearch.params[[pop]] <- get.qsearch.params(post.info, multi.K.params, logten=logten, twoPops=leave_one_out)

        if (is.na(SCDF.qsearch.params[[pop]]$Fh.max)) browser()
    }

    SCDF.qsearch.params
}

calc.quantiles.uniroot <- function(SCDF.qsearch.params, posterior.nu.list, refpopnames, quantiles.vec, logten=FALSE, leave_one_out=FALSE)
{
    quantile.mat <- matrix(nrow = length(refpopnames), ncol = length(quantiles.vec))
    rownames(quantile.mat) <- refpopnames

    for (pop in refpopnames)
    {
        SCDF.qsearch.params.pop <- SCDF.qsearch.params[[pop]]
        post.nu.pop <- sapply(posterior.nu.list, function(x) x[pop,])

        for (i in 1:length(quantiles.vec))
        {
            quantile.mat[pop,i] <- Qhat.uniroot(quantiles.vec[i], SCDF.qsearch.params.pop, post.nu.pop, logten=logten, twoPops=leave_one_out)
        }
    }

    rownames(quantile.mat) <- refpopnames

    quantile.mat
}
