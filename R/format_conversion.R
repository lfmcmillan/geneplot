trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#' Reading data from Genepop.
#'
#' Import Genepop format file and convert to format compatible with GenePlot.
#'
#' File should have locus names at the top of the file, either as single-line
#' list with commas, or as one name per line. Locus names are assumed to stop at
#' the line before the first instance of POP.
#'
#' Population names: By default, Genepop format does not include population
#' names. If the individuals in a pop within the file do not have unique IDs
#' then by default the whole population will be given the ID of the first
#' individual as the population name and then the individuals in that population
#' will be given auto-generated unique IDs. Otherwise, if the individuals in the
#' population do have unique IDs then the populations will be named Pop1, Pop2
#' etc. according to the order in which they appear in the file.
#' If there is a mixture in the file, then pops with unique ID individuals will
#' be named Popx where x is their position in the file, and pops with non-unique
#' ID individuals will be given the ID of their first individual as their popname.
#' Use \code{pop_names} to define the pop names at the point of reading in the file.
#'
#' @param file_path String definining path of the file to read, which can have any
#'     extension.
#'
#' @param header (default TRUE) Boolean, indicate whether there is an additional
#'     descriptive line at the top of the Genepop-format file (TRUE) or whether
#'     the first line is the start of the locus names (FALSE).
#'
#' @param diploid (default TRUE) Boolean, indicates whether data is diploid (TRUE)
#'     or haploid (FALSE).
#'
#' @param digits_per_allele (default 2) Indicates whether data uses 2 or 3 digits
#'     per allele.
#'
#' @param pop_names (default NULL) Character vector (optional). Define the names
#'     of the populations, in the order that they appear in the file.
#'
#' @return A list containing the following components:
#' #' \describe{
#'   \item{\code{locnames}}{Character vector of the locus names.}
#'   \item{\code{pop_data}}{The data, in a data frame, with two columns labelled as 'id' and
#'    'pop', and with two additional columns per locus, labelled in the format
#'    Loc1.a1, Loc1.a2, Loc2.a1, Loc2.a2, etc.}
#' }
#'
#' @importFrom utils flush.console read.table write.table
#'
#' @export
read_genepop_format <- function(file_path,header=TRUE,diploid=TRUE,digits_per_allele=2,
                                pop_names=NULL)
{
    if (!(digits_per_allele %in% 2:3)) stop("digits_per_allele must be 2 or 3.")

    if (header) locStartLine <- 2
    else locStartLine <- 1

    locEndLine <- NA
    con <- file(file_path)
    open(con)
    popStartLines <- vector()
    currentLine <- 1
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0)
    {
        if (trimws(tolower(line)) == 'pop') {
            if (is.na(locEndLine)) locEndLine <- currentLine-1
            popStartLines <- c(popStartLines,currentLine+1)
        }
        currentLine <- currentLine + 1
    }
    close(con)

    ### Read loci names
    if (locEndLine == locStartLine) {
        ## If loci are a single comma-separated line, read in accordingly, but
        ## can't replace failed reading with artificial locus names because don't
        ## know how many loci there are if there was an error reading them in
        rawLocnames <- read.table(file_path,header=FALSE,sep=",",skip=locStartLine-1,nrows=1,colClasses="character")
        nloci <- ncol(rawLocnames)
    } else {
        nloci <- locEndLine-locStartLine+1
        rawLocnames <- tryCatch({
            read.table(file_path,header=FALSE,sep="",skip=locStartLine-1,nrows=nloci,colClasses="character")
        }, error = function(err) { return(paste0("loc-",1:nloci)) })
    }
    locnames <- trimws(unlist(rawLocnames))
    names(locnames) <- NULL

    ### Read allelic data pop by pop
    npops <- length(popStartLines)
    popEndLines <- c(popStartLines[2:npops]-2,currentLine-1)

    rawDataList <- lapply(1:npops,function(pop){
            dat <- tryCatch({
            rawTab <- read.table(file_path, header=FALSE, sep="", colClasses="character",
                       skip = popStartLines[pop]-1, nrows = popEndLines[pop]-popStartLines[pop]+1)

            ## Remove the comma at the end of each identifier that is required in
            ## the Genepop format but not used by GenePlot -- but first check if
            ## the comma is in a separate column or within the text in column 1
            ## (Easypop generates Genepop files where there is a space between the
            ## population index and the comma in each line, and so when read in
            ## the comma appears as a separate column)
            if (rawTab[1,2] == ",") rawTab <- rawTab[,-2]
            else rawTab[,1] <- sapply(rawTab[,1],function(txt) unlist(strsplit(txt,",")))

            ## When adding ids or pop_names, note that the Genepop data is expected
            ## to turn up with pop as the FIRST column and id as the SECOND, so
            ## the next stage of this file will swap the columns around. So at
            ## this stage, add ids as column AFTER pop_names or add pop_names as
            ## the FIRST column.
            ## If IDs are non-unique,make unique ones and use the first ID as the popname
            if (length(unique(rawTab[,1])) < nrow(rawTab))
            {
                popName <- rawTab[1,1]
                if (!is.na(as.numeric(popName))) popName <- paste0("Pop",popName) ## Convert any purely numeric pop_names into "Pop#" to avoid errors in GenePlot

                rawTab <- cbind(rep(popName,times=nrow(rawTab)),paste0(popName,"-",1:nrow(rawTab)),rawTab[,2:ncol(rawTab)])
            }
            else ## If IDs are unique, make up popname
            {
                popName <- paste0("Pop",pop)
                rawTab <- cbind(rep(popName,nrow(rawTab)),rawTab)
            }

            ## Make the first two column names "pop" and "id" so that if some pops
            ## had all unique IDs and other pops didn't, the first two columns still
            ## end up being called "pop" and "id" instead of the code call used
            ## to generate them
            names(rawTab)[1] <- "pop"
            names(rawTab)[2] <- "id"

            ## Now make sure to return the rawTab object
            rawTab

        }, error = function(err) { return(NULL) })
    })

    rawData <- do.call(rbind,rawDataList)

    ### Tidy up allelic data and convert into GenePlot format
    pop_data <- data.frame(id=rawData[,2],pop=rawData[,1]) ## Note switch order of id,pop
    if (diploid){
        for (locIdx in 1:nloci)
        {
            if (locIdx+2 <= ncol(rawData))
            {
                ## Extract the data for each of the alleles, converting it into numeric format instead of factors
                rawAlleles <- rawData[,locIdx+2]

                # if (any(nchar(rawAlleles) < 2*digits_per_allele)) browser()

                allele1 <- as.numeric(sapply(rawAlleles,substr,start=1,stop=digits_per_allele))
                allele2 <- as.numeric(sapply(rawAlleles,substr,start=digits_per_allele+1,stop=2*digits_per_allele))

                pop_data <- cbind(pop_data,allele1,allele2)
            }
        }
        names(pop_data) <- c("id","pop", unlist(lapply(locnames, function(loc) c(paste0(loc,".a1"),paste0(loc,".a2")))))
    }
    else
    {
        pop_data <- cbind(pop_data,as.numeric(rawData[3:ncol(rawData)]))
        names(pop_data) <- c("id","pop", paste0(locnames,".a1"))
    }

    ## Convert the pop column into characters to avoid errors in GenePlot
    pop_data$id <- as.character(pop_data$id)
    pop_data$pop <- as.character(pop_data$pop)

    ## If the user has entered pop_names, use those
    if (!is.null(pop_names)) {
        datpop_names <- unique(pop_data$pop)
        if (length(pop_names) != length(datpop_names)) stop("The length of pop_names does not match the number of pops in the data. Please check and try again.")
        for (pop in 1:length(datpop_names)) {
            pop_data$pop[which(pop_data$pop == datpop_names[pop])] <- pop_names[pop]
        }
    }

    list(pop_data=pop_data,locnames=locnames)
}

write_genepop_format <- function(dat, outfile, locnames_all, delete_loc=NULL, min_loci_per_indiv=1){
    ## genepop.dataformat.func 9/5/12
    ## The output from this function is ready for sending to any Genepop-format
    ## application.
    ## The input data (dat) has missing code 0. It uses 2 digits per allele, and
    ## works only for diploid data.
    ##
    ## To remove specific loci from the analysis, it's easiest to keep
    ## locnames_all as the full list of loci, and specify which loci should be
    ## removed using delete_loc : e.g. delete_loc = c(1, 2) will remove the
    ## first two of the loci named in locnames_all.
    ## The primary purpose of argument locnames_all is to allow different loci
    ## sets for different species, e.g. might want to use locnames_all =
    ## locnames_rat or locnames_all = locnames_stoat.
    ##
    ## If delete_loc=NULL, no loci are removed.
    ##
    ## For Genepop format, the alleles are numbered as follows, within each locus.
    ## 00 = missing code
    ## 01, 02, 03, ...  alleles sampled at the locus.
    ## Information about the allele label is lost, e.g. an individual with
    ## genotype 96, 102 might become 0103 under this scheme, meaning that allele
    ## 96 is given label 01, and allele 102 is given label 03.
    ## The order of labels 01, 02, etc follows the order dictated by the R
    ## function table().
    ##
    ## This function assumes a maximum of 99 alleles at any single locus, which
    ## will be labelled 00 (missing code), 01, 02, ..., 99.  If there are more
    ## than 99 alleles at a locus, the code needs to be updated to 3-digit
    ## format by adding extra leading zeros: 001, 002, ..., 099.
    ##
    ##
    ## Minimum Loci:
    ##
    ## This function only includes individuals with at least min_loci_per_indiv
    ## typed loci - e.g. a reasonable minimum might be 8.  This aims to exclude
    ## altogether individuals with low-quality DNA samples.
    ## For example, if an individual has poor-quality DNA and is typed as a
    ## homozygote at one of its "successful" loci, it might actually be a
    ## heterozygote at that loci but the second allele failed to amplify.
    ##
    ## NOTE that min.loci applies to the FULL data from locnames_all: so if an
    ## individual has 8 present loci out of 10, it counts as 8 for purposes of
    ## exclusion, even if only (say) 5 of these coincide with the actual loci of
    ## interest after delete_loc have been deleted.  The rationale is that this
    ## is an individual with reasonable data quality, evidenced by 8 out of 10
    ## successful loci.
    ##
    ## min_loci_per_indiv defaults to 1, to include all individuals with any
    ## data present.
    ##
    ## EXAMPLE:
    ## The following command will prompt to accept the outfile name
    ## "genepop.gbi.dat.txt", and write results into it:
    ##          genepop.dataformat.func(gbi.dat, min_loci_per_indiv=8)


    ##--------------------------------------------------------------------------------------------
    ## Create default name for outfile if not supplied:
    if(missing(outfile)){
        outfile <- paste("genepop", substitute(dat), "txt", sep=".")
        cat("No name supplied for outfile: use", outfile, "?\n")
        cat("Enter to agree, or Ctrl-C to cancel.\n")
        readline()
    }

    ##--------------------------------------------------------------------------------------------
    ## Exclude poor-quality individuals with < min_loci_per_indiv loci from the
    ## FULL data in locnames_all:
    if(min_loci_per_indiv>1){
        ## Find number of non-missing loci (allele anything except 0) for each
        ## individual.
        ## If there is missing data at a locus, the loc.a1 entry and loc.a2
        ## entry will both be 0.
        ## Use only the loc.a1 entry.
        nloc_dat <- dat[, paste0(locnames_all, ".a1")]
        nloc <- apply(nloc_dat, 1, function(x) length(x[x!=0]))
        cat("WARNING! Excluding all individuals with <", min_loci_per_indiv,
            " loci out of ", length(locnames_all), " loci available.\n")
        ## Remove individuals from dat if they have too few loci:
        dat <- dat[nloc >= min_loci_per_indiv,]
    }

    ##--------------------------------------------------------------------------------------------
    ## Create the names of the desired loci:
    if(is.null(delete_loc)) locnames_foruse <- locnames_all
    else locnames_foruse <- locnames_all[-delete_loc]

    ##--------------------------------------------------------------------------------------------
    ## Isolate the column numbers for the required columns:
    colnames_a1 <- paste0(locnames_foruse, ".a1")
    colnames_a2 <- paste0(locnames_foruse, ".a2")
    colnames_locus <- sort(c(colnames_a1, colnames_a2))
    dat_out <- dat[,c("id", "pop", colnames_locus)]

    ##--------------------------------------------------------------------------------------------
    ## Next order the data frame by population, so that each population
    ## occupies a block of rows in the data frame:

    dat_out <- dat_out[order(dat_out$pop), ]

    ##--------------------------------------------------------------------------------------------
    ## For each locus in turn, replace the allele labels with 00 (missing), 01, 02, ...10, 11, ...
    for(loc in locnames_foruse){
        col_a1 <- paste0(loc, ".a1")
        col_a2 <- paste0(loc, ".a2")
        loc_a1 <- dat_out[, col_a1]
        loc_a2 <- dat_out[, col_a2]
        ## Identify the labels of all alleles in the data for this locus:
        locboth <- c(loc_a1, loc_a2)
        loctab <- table(locboth)
        loctab_names <- names(loctab)
        ## Number of visible alleles (non-null):
        n_visalleles <- length(loctab[loctab_names!="0"])
        new_names <- as.character(1:n_visalleles)
        ## Add leading zeros:
        new_names[1:9] <- paste("0", new_names[1:9], sep="")
        ## Add the missing code if there are any missing data at this locus.
        ## Missing (0) becomes 00.
        if(any(loctab_names=="0")) new_names <- c("00", new_names)
        ## Replace the old allele names in dat_out with the new ones:
        for(al in 1:length(loctab_names)){
            dat_out[,col_a1][as.character(loc_a1)==loctab_names[al]] <- new_names[al]
            dat_out[,col_a2][as.character(loc_a2)==loctab_names[al]] <- new_names[al]
        }
    }
    ## The results are now of the form 00, 01, 02, ... where the order of the new label
    ## is the same as the numeric order of the old label.

    ## Now add a comma after the population name for each individual, and paste
    ## together the columns for each locus.
    ## dat_new looks something like this:
    ## id            loc1     loc2     loc3      ...
    ## Awa,     1010    0611     0303     ...
    ## Awa,     0505    0711     0102     ...  etc.
    # dat_new <- data.frame(id=paste0(dat_out$pop, ","))
    dat_new <- data.frame(id=paste0(dat_out$id,","))
    for(loc in locnames_foruse){
        col_a1 <- paste(loc, ".a1", sep="")
        col_a2 <- paste(loc, ".a2", sep="")
        dat_new[,loc] <- paste0(dat_out[,col_a1], dat_out[,col_a2])
    }

    ## Split into populations:
    dat_split <- split(dat_new, dat_out$pop)


    ##--------------------------------------------------------------------------------------------
    ## Now write the results into outfile in the required Genepop format:
    cat(outfile, file=outfile, "\n", append=F)
    for(ln in paste(sort(locnames_foruse))) cat(ln, "\n", file=outfile, append=T)

    write_output <- function(pop){
        cat("POP\n", file=outfile, append=T)
        write.table(pop, outfile, row.names=F, colnames_locus=F, quote=F, sep="  ", append=T)
    }

    lapply(dat_split, write_output)
    cat("\nFile", outfile, "has been created in Genepop format.\n")
    cat("Genepop On The Web is at http://genepop.curtin.edu.au/\n")
    invisible()
}
