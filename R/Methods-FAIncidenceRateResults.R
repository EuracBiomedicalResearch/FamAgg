## This is the implementation of the FI (familial incidence ratio) and the FSIR
## (familial standardized incidence ratio) from Kerber (1995)

## calculate the time at risk
## startDate: can be either the birth date or the start date for the study.
## endDate: the end date of the study.
## incidenceDate: the date of the incidence.
## deathDate: the date of death.
## the time at risk will be the minimum (earliest) date of (endDate,
## incidenceDate, deathDate) minus the startDate
## need at least endData, in that case endDate is assumed to be "all in one"
## allowNegative: if TRUE return also negative differences. If FALSE, negative
## values are replaced by 0
## Note: we will subtract half of incidence.subtract for each affected individual.
estimateTimeAtRisk <- function(startDate = NULL, startDateFormat = "%Y-%m-%d",
                               endDate = NULL, endDateFormat = "%Y-%m-%d",
                               incidenceDate = NULL,
                               incidenceDateFormat = "%Y-%m-%d",
                               deathDate = NULL, deathDateFormat = "%Y-%m-%d",
                               allowNegative = FALSE, affected = NULL,
                               incidenceSubtract = 0.5) {
    if (is.null(startDate))
        stop("'startDate' has to be specified (either the start of the study,",
             " or the birth date).")
    if (!is(startDate, "Date")) {
        message("Transforming startDate into Date format...", appendLF=FALSE)
        startDate <- as.Date(startDate, format=startDateFormat)
        message("OK")
    }
    ## check if we've got at least one of the following:
    if (is.null(endDate)) {
        stop("'endDate' has to be specified!")
    } else {
        if(is.null(incidenceDate) & is.null(deathDate))
            message("'endDate' is assumed to specify, for each individual,",
                    " either the end year of the study, the individual's death",
                    " date or the individual's incidence date.")
        if (length(endDate) != length(startDate))
            stop("'startDate' and 'endDate' have to have the same length!")
    }
    message("Transforming endDate into Date format...", appendLF=FALSE)
    endDate <- as.Date(endDate, format=endDateFormat)
    message("OK")
    if (!is.null(affected)) {
        if (length(affected) != length(startDate)) {
            affected <- NULL
            warning("length of affected does not match length of startData; ",
                    "dropping.")
        } else {
            affected <- as.logical(affected)
        }
    }
    ## if we've got an incidenceDate, process that.
    if (!is.null(incidenceDate)) {
        message("Transforming incidenceDate into Date format...",
                appendLF=FALSE)
        incidenceDate <- as.Date(incidenceDate, format=incidenceDateFormat)
        message("OK")
        if(length(incidenceDate) != length(startDate))
            stop("'incidenceDate', 'startDate' and 'endDate' all have to have ",
                 "the same length!")
        ## now get the earlier time point: incidence or end:
        DF <- data.frame(endDate, incidenceDate)
        endDate <- as.Date(apply(DF, MARGIN=1, min, na.rm=TRUE))
        affected <- rep(FALSE, length(incidenceDate))
        affected[!is.na(incidenceDate)] <- TRUE
    }
    if (!is.null(deathDate)) {
        message("Transforming deathDate into Date format...", appendLF=FALSE)
        deathDate <- as.Date(deathDate, format=deathDateFormat)
        message("OK")
        if(length(deathDate) != length(startDate))
            stop("'deathDate', 'startDate' and 'endDate' all have to have the ",
                 "same length!")
        ## now get the earlier time point: incidence or end:
        DF <- data.frame(endDate, deathDate)
        endDate <- as.Date(apply(DF, MARGIN=1, min, na.rm=TRUE))
    }
    ## OK, now we're set.
    Diff <- endDate - startDate
    if (!is.null(affected)) {
        Diff[affected] <- Diff[affected] - incidenceSubtract
    }
    Diff <- as.numeric(Diff)
    if(!allowNegative)
        Diff[which(Diff < 0)] <- 0
    Diff
}

## takes a numeric vector of ages as input and slices it into
sliceAge <- function(x, slices=c(0, 40, Inf)){
    SliceDiff <- diff(slices)
    Res <- lapply(x, function(z){
        tmp <- z - slices
        tmp <- tmp[1:(length(tmp)-1)]
        bm <- which(tmp > SliceDiff)
        tmp[bm] <- SliceDiff[bm]
        tmp[which(tmp < 0)] <- 0
        return(tmp)
    })
    Res <- do.call(rbind, Res)
    CN <- paste0(slices[1:(length(slices)-1)], ", ", slices[2:length(slices)])
    CN <- paste0("(", CN, "]")
    colnames(Res) <- CN
    return(Res)
}

##****************************************************************************
##
##  timeAtRisk getter/setter method.
##
setMethod("timeAtRisk", "FAIncidenceRateResults",
          function(object){
              ## we're also checking that timeAtRisk, stored in its internal
              ## slot has the correct length!
              if(length(object@timeAtRisk)==0)
                  return(object@timeAtRisk)
              if(length(object@timeAtRisk) != length(object$id))
                  stop("Length of 'timeAtRisk' does not match the number of ",
                       "individuals in the pedigree!")
              return(object@timeAtRisk)
          })
setReplaceMethod("timeAtRisk", "FAIncidenceRateResults",
                 function(object, value) {
                     ## check if we've got the pedigree, otherwise does not
                     ## make sense to store the timeAtRisk.
                     if(length(object$id) == 0)
                         stop("No pedigree data in 'object'! Add the pedigree ",
                              "data before setting 'timeAtRisk'.")
                     if(length(object$id) != length(value))
                         stop("Length of 'timeAtRisk' does not match the ",
                              "number of individuals in the pedigree!")
                     object@timeAtRisk <- value
                     object
                 })


##****************************************************************************
##
##  familialIncidenceRate getter method.
##
setMethod("familialIncidenceRate", "FAIncidenceRateResults",
          function(object, ...){
              if(length(object@sim) == 0)
                  return(NULL)
              return(object@sim$fir)
          })

##****************************************************************************
##
##  $ getter method.
##
setMethod("$", "FAIncidenceRateResults", function(x, name){
    if(name == "fir" | name == "tar" | name == "timeAtRisk"
       | name == "pvalue"){
        ## extract local stuff
        if(length(x@sim) == 0){
            warning("No calculation result available!")
            return(NULL)
        }
        if(name == "fir")
            return(x@sim$fir)
        if(name == "tar")
            return(timeAtRisk(x))
        if(name == "timeAtRisk")
            return(timeAtRisk(x))
        if(name == "pvalue")
            return(x@sim$pvalue)
    }else{
        callNextMethod()
    }
})


##****************************************************************************
##
##  runSimulation
##
##  We can pass optional argument 'lowMem' to the downstream function that
##  requires less memory and runs (eventually) faster.
##  The second optional argument is 'prune' that removes unconnected
##  individuals from the pedigree prior testing.
setMethod("runSimulation", "FAIncidenceRateResults",
          function(object, nsim=50000, timeAtRisk=NULL,
                   strata=NULL, ...){
              if(length(trait(object)) == 0)
                  stop("No trait information available!")
              SimRes <- .FRSimulation(ped=pedigree(object),
                                      kin=kinship(object),
                                      trait=trait(object),
                                      timeAtRisk=timeAtRisk,
                                      nsim=nsim,
                                      strata=strata,
                                      prune=FALSE,
                                      ...)
              timeAtRisk(object) <- timeAtRisk
              object@nsim <- nsim
              object@sim <- SimRes
              return(object)
          })

##****************************************************************************
##
##  trait<-
##
##  Calling the trait replecement method from FAResult and in addition
##  reset the simulation result.
setReplaceMethod("trait", "FAIncidenceRateResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})

##****************************************************************************
##
##  familial incidence rate along with Monte Carlo simulations to
##  estimate the significance of the ratio.
##  Note: we skipped here the perFamilyTest option as that makes a) no sense,
##  and b) would make the simulations eventually problematic (sampling from
##  small sample sets).
##
##  lowMem: if TRUE a less memory demanding code is run; density information
##          is then however not available.
##****************************************************************************
.FRSimulation <- function(ped=NULL, kin=NULL, trait=NULL, timeAtRisk=NULL,
                          prune=TRUE, nsim=50000, strata=NULL,
                          lowMem=FALSE, ...){
    ## Input argument checking.
    if(is.null(ped))
        stop("No pedigree submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    if(is.null(trait))
        stop("No trait information submitted!")
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there",
             " are rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, AFF=trait)
    allIds <- as.character(ped$id)
    if(is.null(timeAtRisk)){
        stop("Argument 'timeAtRisk' has to be submitted!")
    }
    if(length(timeAtRisk) != nrow(ped))
        stop("Argument 'timeAtRisk' has to have the same length than there",
             " are rows (individuals) in the pedigree 'ped'!")
    timeAtRisk <- as.numeric(timeAtRisk)
    ped <- cbind(ped, TAR=timeAtRisk)
    if(!is.null(strata)){
        if(length(strata)!=nrow(ped))
            stop("Argument 'strata' has to have the same length than there",
                 " are individuals in the pedigree!")
        ped <- cbind(ped, STRAT=strata)
    }
    ## NOTE: We're not removing singletons here, as we are anyway removing them
    ## further down!
    ## Subset the data set to individuals with valid values for trait, timeAtRisk
    ## and strata.
    nas <- is.na(ped$AFF)
    message("Cleaning data set (got in total ", nrow(ped), " individuals):")
    message(" * not phenotyped individuals...", appendLF=FALSE)
    if(any(nas)){
        ped <- ped[!nas, , drop=FALSE]
        message(" ", sum(nas), " removed.")
    }else{
        message(" none present.")
    }
    ## time at risk
    nas <- is.na(ped$TAR)
    message(" * individuals with unknown time at risk...", appendLF=FALSE)
    if(any(nas)){
        ped <- ped[!nas, , drop=FALSE]
        message(" ", sum(nas), " removed.")
    }else{
        message(" none present.")
    }
    ## Strata
    if(!is.null(strata)){
        nas <- is.na(ped$STRAT)
        message(" * individuals without valid strata values...", appendLF=FALSE)
        if(any(nas)){
            ped <- ped[!nas, , drop=FALSE]
            message(" ", sum(nas), " removed.")
        }else{
            message(" none present.")
        }

    }
    ## Anyway removing singletons here, since they result in NA values!
    message(" * singletons (also caused by previous subsetting)...",
            appendLF=FALSE)
    origSize <- nrow(ped)
    suppressMessages(
        ped <- removeSingletons(ped)
    )
    if(nrow(ped) != origSize){
        message(" ", origSize-nrow(ped), " removed.")
    }else{
        message(" none present.")
    }
    message("Done")

    ## Calculate the kinship for the observed.
    idx <- match(as.character(ped$id), colnames(kin))
    if(any(is.na(idx)))
        stop("Some of the individuals in the pedigree are missing in the ",
             "kinship matrix!")
    kin <- as.matrix(kin[idx, idx])
    ## Remove self-self kinship and set all affected > 1 to 1.
    diag(kin) <- 0
    Denomi <- ped$TAR %*% kin
    affected <- rep(0, nrow(ped))
    affected[ped$AFF > 0] <- 1
    nAff <- sum(affected)
    nAll <- length(affected)
    sampleFrom <- 1:length(affected)
    obsFr <- as.vector((affected %*% kin) / Denomi)
    names(obsFr) <- ped$id
    ##
    ## Memory saving version
    if(lowMem){
        largerEqual <- rep(0, length(obsFr))
        if(is.null(strata)){
            for(i in 1:nsim){
                simIdx <- sample(sampleFrom, nAff)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                tmp <- as.vector((simaffs %*% kin) / Denomi)
                tmp[is.na(tmp)] <- 0
                largerEqual <- largerEqual + (tmp >= obsFr)
            }
        }else{
            affStrataCounts <- table(ped[ped$AFF > 0 , "STRAT"])
            affStrataCounts <- affStrataCounts[affStrataCounts > 0]
            for(i in 1:nsim){
                simIdx <- stratsample(ped$STRAT, counts=affStrataCounts)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                tmp <- as.vector((simaffs %*% kin) / Denomi)
                tmp[is.na(tmp)] <- 0
                largerEqual <- largerEqual + (tmp >= obsFr)
            }
        }
        PVals <- largerEqual/nsim
        denses <- NULL
    }else{
        ## full blown version
        ## Simulations: generate random sets of affected individuals, same size
        ## then the observed individuals.
        if(is.null(strata)){
            SimFr <- lapply(1:nsim, function(z){
                simIdx <- sample(sampleFrom, nAff)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                return(as.vector((simaffs %*% kin) / Denomi))
                ## Would be less memory demanding to just return whether the
                ## expected is larger than the observed
            })
        }else{
            affStrataCounts <- table(ped[ped$AFF > 0 , "STRAT"])
            affStrataCounts <- affStrataCounts[affStrataCounts > 0]
            SimFr <- lapply(1:nsim, function(z){
                simIdx <- stratsample(ped$STRAT, counts=affStrataCounts)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                return(as.vector((simaffs %*% kin) / Denomi))
            })
        }
        ## Performing the stuff on the large matrix.
        ## Want to compare the observed FR against the simulated FOR EACH
        ## INDIVIDUAL.
        ## Why: that way we keep the kinship the same and just calculate an
        ## expexted FR for random occurence of cases.
        SimFr <- do.call(cbind, SimFr)
        PVals <- rowSums(SimFr >= obsFr)/nsim
        denses <- apply(SimFr, MARGIN=1, density)
    }
    ## Result should match all ids in the (original) pedigree.
    fr <- rep(NA, length(allIds))
    names(fr) <- allIds
    idx <- match(names(obsFr), allIds)
    fr[idx] <- obsFr
    ## p-values
    pv <- rep(NA, length(allIds))
    names(pv) <- allIds
    pv[idx] <- PVals
    ## density
    if(!is.null(denses)){
        dens <- vector("list", length(allIds))
        names(dens) <- allIds
        dens[idx] <- denses
        Res <- list(fir=fr, pvalue=pv, expDensity=dens)
    }else{
        Res <- list(fir=fr, pvalue=pv)
    }
    Res
}

## prune: remove un-connected individuals from the pedigree AFTER removing
##        individuals without valid value in trait or a missing value in
##        timeAtRisk.
## UPDATE: disable perFamilyTest!!!
.FR <- function(ped=NULL, kin=NULL, trait=NULL, timeAtRisk=NULL,
                prune=TRUE, ...){
    perFamilyTest <- FALSE
    ## subset the data to all the individuals for which we do have the
    ## required data.
    ## * trait not NA
    ## * timeAtRist not NA
    ## * did not get removed after pruning for un-connected individuals.
    if(is.null(ped))
        stop("No pedigree submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    if(is.null(trait))
        stop("No trait information submitted!")
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there",
             " are rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, AFF=trait)
    if(is.null(timeAtRisk)){
        stop("Argument 'timeAtRisk' has to be submitted!")
    }
    if(length(timeAtRisk) != nrow(ped))
        stop("Argument 'timeAtRisk' has to have the same length than",
             " there are rows (individuals) in the pedigree 'ped'!")
    timeAtRisk <- as.numeric(timeAtRisk)
    ped <- cbind(ped, TAR=timeAtRisk)
    if(perFamilyTest){
        entityIs <- "family"
    }else{
        entityIs <- "pedigree"
        ped[, "family"] <- 1
    }
    ## NOTE: don't do that here, since we're going to remove them anyway
    ## further below!
    pedL <- split(ped, ped$family)
    Res <- lapply(pedL, function(z){
        allIds <- as.character(z$id)
        ## remove NAs:
        ## NA in affected/trait
        nas <- is.na(z$AFF)
        message("Cleaning data set (got in total ", nrow(z), " individuals):")
        message(" * not phenotyped individuals...", appendLF=FALSE)
        if(any(nas)){
            z <- z[!nas, , drop=FALSE]
            message(" ", sum(nas), " removed.")
        }else{
            message(" none present.")
        }
        ## NA in time at risk:
        nas <- is.na(z$TAR)
        message(" * individuals with unknown time at risk...", appendLF=FALSE)
        if(any(nas)){
            z <- z[!nas, , drop=FALSE]
            message(" ", sum(nas), " removed.")
        }else{
            message(" none present.")
        }
        message(" * singletons (also caused by previous subsetting)...", appendLF=FALSE)
        origSize <- nrow(z)
        suppressMessages(
            z <- removeSingletons(z)
        )
        if(nrow(z) != origSize){
            message(" ", origSize-nrow(z), " removed.")
        }else{
            message(" none present.")
        }
        message("Done")

        ## Subset the kinship matrix and ensure correct ordering.
        idx <- match(as.character(z$id), colnames(kin))
        if(any(is.na(idx)))
            stop("Some of the individuals in the pedigree are missing in the ",
                 "kinship matrix!")
        kinSub <- as.matrix(kin[idx, idx])
        ## remove self-self kinship.
        diag(kinSub) <- 0
        affected <- rep(0, nrow(z))
        affected[z$AFF > 0] <- 1
        ## The magic formula (3) in Kerber 1995 as an R matrix multiplication.
        FR <- as.vector((affected %*% kinSub)/(z$TAR %*% kinSub))
        names(FR) <- z$id
        retFR <- rep(NA, length(allIds))
        names(retFR) <- allIds
        retFR[names(FR)] <- FR
        return(retFR)
    })
    Res
}


setMethod("[", "FAIncidenceRateResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAIncidenceRateResults object is not supported!")
})

##****************************************************************************
##
##  plotPed
##
##  the plot ped can be done for individuals and for families. In essence we
##  are just calling the plotPed from the FAData object, eventually highlighting
##  the individual and providing the FR as numeric values below the id.
setMethod("plotPed", "FAIncidenceRateResults",
          function(object, id=NULL, family=NULL, filename=NULL,
                   device="plot", only.phenotyped=FALSE, ...){
              if(is.null(id) & is.null(family))
                  stop("Either id or family has to be provided!")
              if(length(object@sim) == 0){
                  ## Means we don't have any results...
                  FR <- rep(NA, length(object$id))
              }else{
                  FR <- object@sim$fir
              }
              callNextMethod(object=object, id=id, family=family,
                             filename=filename,
                             device=device, label1=FR, proband.id=id,
                             only.phenotyped=only.phenotyped, ...)
          })


##****************************************************************************
##
##  show
##
setMethod("show", "FAIncidenceRateResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    if(length(object@sim) > 0){
        cat(paste0(" * Mean familial incidence ratio: ",
                   format(mean(object@sim$fir, na.rm=TRUE), digits=2), ".\n"))
        cat(paste0(" * Mean time at risk: ",
                   format(mean(timeAtRisk(object), na.rm=TRUE), digits=2), ".\n"))
    }else{
        cat(" * No simulation results available yet.\n")
    }
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})

##****************************************************************************
##
##  plotRes
##
setMethod("plotRes", "FAIncidenceRateResults",
          function(object, id=NULL, family=NULL, addLegend=TRUE,
                   type="density", ...){
              type <- match.arg(type, c("density", "hist"))
              if(type == "hist")
                  stop("Type 'hist' is not supported (yet).")
              if(length(object@sim) == 0)
                  stop("No analysis performed yet!")
              if(is.null(id))
                  stop("Argument 'id' is required.")
              id <- as.character(id)
              if(!is.null(family))
                  stop("Argument 'family' is not supported.")
              ## check if we've got id at all.
              fr <- object@sim$fir
              if(!any(names(fr) == id))
                  stop("Individual with id ", id, " not found ",
                       "in the pedigree.")
              obsFr <- fr[id]
              if(is.na(obsFr))
                  stop(paste0("No Familial Incidence Rate value for individual ",
                              id, " calculated."))
              fam <- family(object, id=id)[1, "family"]
              if(type == "density"){
                  ## Let's see whether we have the required information available.
                  if(!any(names(object@sim) == "expDensity"))
                      stop("Distribution of familial incidence rates from ",
                           "the simulation runs not available. You need to",
                           " run 'runSimulation' without optional argument",
                           " 'lowMem=TRUE'.")
                  dens <- object@sim$expDensity[[id]]
                  XL <- range(c(range(dens$x, na.rm=TRUE), obsFr), na.rm=TRUE)
                  plot(dens, main=paste0("Individual: ", id, ", family: ",
                                         fam), xlab="Familial incidence rate",
                       type="h", lwd=3, col="lightgrey", xlim=XL)
                  points(dens, col="grey", type="l", lwd=2)
              }
              Blue <- "#377EB8"
              abline(v=obsFr, col=Blue)
              if(addLegend){
                  legend("topright",
                         legend=c(paste0("fir     : ", format(obsFr, digits=2)),
                                  paste0("p-value : ",
                                         format(object@sim$pvalue[id],
                                                digits=2))
                                  ))
              }
          })

##****************************************************************************
##
##  result
##
setMethod("result", "FAIncidenceRateResults", function(object, method="BH"){
    method <- match.arg(method, p.adjust.methods)
    Trait <- trait(object, na.rm=TRUE)
    TraitName <- object@traitname
    if(length(TraitName) == 0)
        TraitName <- NA
    if(length(object@sim) == 0){
        MyRes <- data.frame(trait_name=TraitName,
                            total_phenotyped=length(Trait),
                            total_affected=sum(Trait!=0),
                            total_tested=0,
                            id=NA,
                            family=NA,
                            fir=NA,
                            pvalue=NA,
                            padj=NA,
                            check.names=FALSE, stringsAsFactors=FALSE)
        warning(paste0("No simulation data available! Please run a simulation ",
                       "using the 'runSimulation' method on this, or the ",
                       "'familialIncidenceRateTest' on a 'FAData' object."))
        return(MyRes)
    }
    ## @sim$fir has to be in the same order than the pedigree.
    nind <- length(object$id)
    MyRes <- data.frame(trait_name=rep(TraitName, nind),
                        total_phenotyped=rep(length(Trait), nind),
                        total_affected=rep(sum(Trait!=0), nind),
                        total_tested=rep(sum(!is.na(object@sim$fir)), nind),
                        id=as.character(object$id),
                        family=as.character(object$family),
                        fir=object@sim$fir,
                        pvalue=object@sim$pvalue,
                        padj=p.adjust(object@sim$pvalue, method=method),
                        check.names=FALSE, stringsAsFactors=FALSE)
    MyRes <- MyRes[order(MyRes$pvalue, -MyRes$fir), ]
    rownames(MyRes) <- as.character(MyRes$id)
    MyRes
})



