## Methods and functions related to the FSIR from Kerber.





##************************************************
##
## FSIR Kerber 1995
## the familial standardized incidence ratio.
## What do we need to calculate it:
## * affected vector: can get the cases from that.
## * the kinship matrix: setting diag(kin) <- 0, so we exclude self-self kinship.
## * strata: some sort of stratification of cases.
## * time in strata: related to the strata: the time each individual spend in strata.
## * the population incidence rate per stratum: that's really tricky! how do we get the
##   population incidence
##   if we have only the cohort available... sampling? that will not work, if, only with
##   replacement... but still
##
## Side note: what if we use the incidence ratio from the whole pedigree? the FSIR is anyway
##            just calculated based
##            on individuals from the same family (got a kinship value of 0 for other
##            individuals). Such lambda would
##            represent the expected number of affected individuals in the family given the
##            data from the whole pedigree.
##************************************************
## Note: this function assumes that all elements are in the correct order!
## Arguments:
## * affected: 0, 1 coding, length corresponds to the phenotyped; length n.
## * kin: nxn kinship matrix. ordering has to match individuals in affected.
## * lambda: population incidence rate. numeric vector, length has to match the number of strata.
## * timeInStrata: ncol=length(lambda), nrow=n. has to have the same unit than lambda.
doFsir <- function(affected, kin, lambda, timeInStrata){
    if(missing(affected))
        stop("affected missing!")
    if(missing(kin))
        stop("kinship matrix kin missing!")
    if(missing(lambda))
        stop("lambda missing!")
    if(missing(timeInStrata))
        stop("timeInStrata missing!")
    if(length(unique(c(length(affected), nrow(kin), ncol(kin), nrow(timeInStrata)))) != 1)
        stop("Length of affected has to match the number of rows and cols of",
             " kin and the number of rows of timeInStrata!")

    ## match names of lambda with colnames of timeInStrata:
    if(is.null(names(lambda)))
        stop("lambda has to be a named numeric vector with the names corresponding",
             " to the colnames of timeInStrata!")
    if(is.null(colnames(timeInStrata)))
        stop("timeInStrata has to be a matrix with column names corresponding to",
             " the names of argument lambda!")
    if(length(lambda) != ncol(timeInStrata))
        stop("length of lambda does not match number of columns of timeInStrata!")
    if(!all(names(lambda) %in% colnames(timeInStrata)))
        stop("names of lambda do not match with colnames of timeInStrata!")
    lambda <- lambda[colnames(timeInStrata)]

    ## setting the diagonal of kin to 0 to avoid self-self kinship
    diag(kin) <- 0
    ## calculate nominator of Kerber 1995 formula (4): the observed cases
    ## this results in a vector length affected, with 0 if the individual is not related to any
    ## affected.
    nomi <- as.numeric(affected %*% kin)
    ## funny thing is the denominator that lists the expected cases.
    ## lambda should be a vector of length k (k=number of strata), timeInStrat a matrix
    ## with k columns, n rows with n being the number of individuals and kin a nxn matrix.
    ## lambda and timeInStrat should be in the same unit,
    ## i.e. either proportion, or proportion per person-year.
    denomi <- lambda %*% t(timeInStrata) %*% kin
    if(nrow(denomi) > 1)
        stop("Did not get the expected single-row matrix!")
    denomi <- as.numeric(denomi)
    res <- nomi/denomi
    names(res) <- names(affected)
    ##return(list(nomi=nomi, denomi=denomi))
    return(res)
}


##****************************************************************************
##
##  runSimulation
##
setMethod("runSimulation", "FAStdIncidenceRateResults",
          function(object, nsim=50000, lambda=NULL, timeInStrata=NULL,
                   strata=NULL, ...){
              if(length(trait(object)) == 0)
                  stop("No trait information available!")
              ## Check input parameter:
              if(is.null(timeInStrata))
                  stop("Required argument 'timeInStrata' missing!")
              if(is.null(lambda))
                  stop("Required argument 'lambda' missing!")
              if(is.null(ncol(timeInStrata))){
                  ## Special case: numeric vector, no strata.
                  timeInStrata <- matrix(timeInStrata, ncol=1)
                  colnames(timeInStrata) <- "A"
                  names(lambda)[1] <- "A"
              }
              ## Start defining input arguments for the downstream function
              affected <- trait(object)
              kin <- kinship(object)
              ## Check correct dimensions etc.
              if(nrow(timeInStrata) != length(affected))
                  stop("Number of rows of 'timeInStrata' (",
                       nrow(timeInStrata), ") has to match the number ",
                       "of individuals in the pedigree (",
                       length(affected), ")!")
              if(ncol(timeInStrata) != length(lambda))
                  stop("Number of columns of 'timeInStrata' (",
                       ncol(timeInStrata), ") has to match the length ",
                       "of 'lambda' (", length(lambda), ")!")
              if(any(colnames(timeInStrata) != names(lambda)))
                  stop("Column names of 'timeInStrata' have to match ",
                       "names of lambda!")
              if(!is.null(strata)){
                  if(length(strata) != length(affected))
                      stop("Length of 'strata' (", length(strata), ") ",
                           "has to match the number of individuals in ",
                           "the pedigree (", length(affected), ")!")
              }
              object@timeInStrata <- timeInStrata
              object@lambda <- lambda
              ## Remove entries with NAs
              ## trait
              trait <- trait(object)   # that way we ensure that we have the same ordering.
              kin <- kinship(object)
              ## Just to be on the save side... ensure that the ordering of id/Trait
              ## matches the kin
              kin <- kin[names(trait), names(trait)]
              diag(kin) <- 0
              ## Now start subsetting the data:
              message("Cleaning data set (got in total ", nrow(kin), " individuals):")
              message(" * not phenotyped individuals...", appendLF=FALSE)
              ## * NA in trait
              nas <- is.na(trait)
              if(any(nas)){
                  ## warning(paste0("Excluding ", sum(nas),
                  ##                " individuals because of a missing value in the trait."))
                  trait <- trait[!nas]
                  kin <- kin[!nas, !nas]
                  timeInStrata <- timeInStrata[!nas, , drop=FALSE]
                  if(!is.null(strata))
                      strata <- strata[!nas]
                  message(" ", sum(nas), " removed.")
              }else{
                  message(" none present.")
              }
              ## * NA in timeInStrata.
              nas <- apply(timeInStrata, MARGIN=1, function(z){
                  any(is.na(z))
              })
              message(" * individuals with missing time in strata...", appendLF=FALSE)
              if(any(nas)){
                  ## warning(paste0("Excluding ", sum(nas),
                  ##                " individuals because of a missing value in timeInStrata."))
                  trait <- trait[!nas]
                  kin <- kin[!nas, !nas]
                  timeInStrata <- timeInStrata[!nas, , drop=FALSE]
                  if(!is.null(strata))
                      strata <- strata[!nas]
                  message(" ", sum(nas), " removed.")
              }else{
                  message(" none present.")
              }
              ## * NA in strata
              if(!is.null(strata)){
                  nas <- is.na(strata)
                  message(" * individuals without valid strata values...", appendLF=FALSE)
                  if(any(nas)){
                      ## warning(paste0("Excluding ", sum(nas),
                      ##                " individuals because of a missing value in strata."))
                      trait <- trait[!nas]
                      kin <- kin[!nas, !nas]
                      timeInStrata <- timeInStrata[!nas, , drop=FALSE]
                      strata <- strata[!nas]
                      message(" ", sum(nas), " removed.")
                  }else{
                      message(" none present.")
                  }
              }
              ## Anyway removing singletons here, since they result in NA values!
              message(" * singletons (also caused by previous subsetting)...", appendLF=FALSE)
              ## * Not related, i.e. individuals with a kinship sum of 0
              nas <- colSums(kin) == 0
              if(any(nas)){
                  ## warning(paste0("Excluding ", sum(nas),
                  ##                " individuals because they do not share kinship with",
                  ##                " any individual in the pedigree."))
                  trait <- trait[!nas]
                  kin <- kin[!nas, !nas]
                  timeInStrata <- timeInStrata[!nas, , drop=FALSE]
                  if(!is.null(strata))
                      strata <- strata[!nas]
                  message(" ", sum(nas), " removed.")
              }else{
                  message(" none present.")
              }
              message("Done")

              ## OK, now run the test...
              Sim <- doFsirSimulation(affected=trait, kin=kin, lambda=lambda,
                                      timeInStrata=timeInStrata,
                                      nsim=nsim, strata=strata, ...)
              fsirs <- Sim$fsir
              allFsirs <- rep(NA, length(object$id))
              names(allFsirs) <- object$id
              idx <- match(names(fsirs), object$id)
              allFsirs[idx] <- fsirs
              allPvals <- rep(NA, length(object$id))
              names(allPvals) <- object$id
              allPvals[idx] <- Sim$pvalue
              if(length(Sim$expDensity) > 0){
                  allExps <- vector("list", length(object$id))
                  names(allExps) <- object$id
                  allExps[idx] <- Sim$expDensity
                  Res <- list(fsir=allFsirs, pvalue=allPvals, expDensity=allExps)
              }else{
                  Res <- list(fsir=allFsirs, pvalue=allPvals)
              }
              object@sim <- Res
              object@nsim <- nsim
              return(object)
          })



##****************************************************************************
##
##  FSIR along with Monte Carlo simulations to
##  estimate significances.
##  We are NOT checking for NAs here, so that has to be done upstream!
##  Also, we assume that everything is in the right order.
##
##  lowMem: if TRUE a less memory demanding code is run; density information
##          is then however not available.
##****************************************************************************
doFsirSimulation <- function(affected, kin, lambda, timeInStrata,
                             nsim=50000, strata=NULL, lowMem=FALSE){
    if(missing(affected))
        stop("affected missing!")
    if(missing(kin))
        stop("kinship matrix kin missing!")
    if(missing(lambda))
        stop("lambda missing!")
    if(missing(timeInStrata))
        stop("timeInStrata missing!")
    if(length(unique(c(length(affected), nrow(kin), ncol(kin), nrow(timeInStrata)))) != 1)
        stop("Length of affected has to match the number of rows ",
             "and cols of kin and the number of rows of timeInStrata!")
    if(!is.null(strata)){
        if(length(strata) != length(affected))
            stop("Length of arguments 'strata' and 'affected' have to match!")
    }
    ## match names of lambda with colnames of timeInStrata:
    if(is.null(names(lambda)))
        stop("lambda has to be a named numeric vector with the names",
             " corresponding to the colnames of timeInStrata!")
    if(is.null(colnames(timeInStrata)))
        stop("timeInStrata has to be a matrix with column names",
             " corresponding to the names of argument lambda!")
    if(length(lambda) != ncol(timeInStrata))
        stop("length of lambda does not match number of columns of timeInStrata!")
    if(!all(names(lambda) %in% colnames(timeInStrata)))
        stop("names of lambda do not match with colnames of timeInStrata!")
    lambda <- lambda[colnames(timeInStrata)]

    ## setting the diagonal of kin to 0 to avoid self-self kinship
    diag(kin) <- 0
    ## calculate nominator of Kerber 1995 formula (4): the observed cases
    ## this results in a vector length affected, with 0 if the individual is not
    ## related to any affected.
    nomi <- as.numeric(affected %*% kin)
    ## funny thing is the denominator that lists the expected cases.
    ## lambda should be a vector of length k (k=number of strata), timeInStrat a
    ## matrix with k columns, n rows with n being the number of individuals and
    ## kin a nxn matrix. lambda and timeInStrat should be in the same unit,
    ## i.e. either proportion, or proportion per person-year.
    denomi <- lambda %*% t(timeInStrata) %*% kin
    if(nrow(denomi) > 1)
        stop("Did not get the expected single-row matrix!")
    denomi <- as.numeric(denomi)
    obsFsir <- nomi/denomi
    names(obsFsir) <- names(affected)
    ##
    ## Performing the simulation.
    sampleFrom <- 1:length(affected)
    nAff <- sum(affected)   ## affected has to be encoded 1,0 or TRUE,FALSE
    nAll <- length(affected)
    if(lowMem){
        largerEqual <- rep(0, nAll)
        if(is.null(strata)){
            for(i in 1:nsim){
                simIdx <- sample(sampleFrom, nAff)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                tmp <- as.vector((simaffs %*% kin) / denomi)
                tmp[is.na(tmp)] <- 0
                largerEqual <- largerEqual + (tmp >= obsFsir)
            }
        }else{
            affStrataCounts <- table(strata[affected > 0])
            affStrataCounts <- affStrataCounts[affStrataCounts > 0]
            for(i in 1:nsim){
                simIdx <- stratsample(strata, counts=affStrataCounts)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                tmp <- as.vector((simaffs %*% kin) / denomi)
                tmp[is.na(tmp)] <- 0
                largerEqual <- largerEqual + (tmp >= obsFsir)
            }
        }
        PVals <- largerEqual/nsim
        denses <- NULL
    }else{
        ## memory demanding and "slower" version
        if(is.null(strata)){
            SimFsir <- lapply(1:nsim, function(z){
                simIdx <- sample(sampleFrom, nAff)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                return(as.vector((simaffs %*% kin) / denomi))
            })
        }else{
            affStrataCounts <- table(strata[affected > 0])
            affStrataCounts <- affStrataCounts[affStrataCounts > 0]
            SimFsir <- lapply(1:nsim, function(z){
                simIdx <- stratsample(strata, counts=affStrataCounts)
                simaffs <- rep(0, nAll)
                simaffs[simIdx] <- 1
                return(as.vector((simaffs %*% kin) / denomi))
            })
        }
        SimFsir <- do.call(cbind, SimFsir)
        PVals <- rowSums(SimFsir >= obsFsir)/nsim
        denses <- apply(SimFsir, MARGIN=1, density)
    }
    return(list(fsir=obsFsir, pvalue=PVals, expDensity=denses))
}


## convert a factor to a matrix, levels
factor2matrix <- function(x){
    if(class(x)!="factor")
        x <- factor(x)
    Ls <- lapply(levels(x), function(z){
        return(as.numeric(x == z))
    })
    Mat <- do.call(cbind, Ls)
    colnames(Mat) <- levels(x)
    if(!is.null(names(x)))
        rownames(Mat) <- names(x)
    return(Mat)
}



##****************************************************************************
##
##  show
##
setMethod("show", "FAStdIncidenceRateResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    if(length(object@sim) > 0){
        cat(paste0(" * Mean familial standardized incidence ratio: ",
                   format(mean(object@sim$fsir, na.rm=TRUE), digits=2), ".\n"))
        cat(paste0(" * lambda: \n"))
        Lam <- object@lambda
        for(i in 1:length(lambda)){
            cat(paste0("   - ", names(lambda)[i], ": ", lambda[i], "\n"))
        }
        cat(paste0())
    }else{
        cat(" * No simulation results available yet.\n")
    }
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})

##****************************************************************************
##
##  plotRes
##
setMethod("plotRes", "FAStdIncidenceRateResults",
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
              ## check if id is valid
              fsirs <- object$fsir
              if(!any(names(fsirs) == id))
                  stop("Individual with id ", id, " not found ",
                       "in the pedigree.")
              obsFsir <- fsirs[id]
              if(is.na(obsFsir))
                  stop("No Familial Standardized Incidence Rate (FSIR) ",
                       "calculated for individual ", id, ".")
              fam <- family(object, id=id)[1, "family"]
              if(type == "density"){
                  ## Let's see whether we have the required information available.
                  if(!any(names(object@sim) == "expDensity"))
                      stop("Distribution of familial standardiced incidence rates ",
                           "from the simulation runs not available. You need to",
                           " run 'runSimulation' without optional argument",
                           " 'lowMem=TRUE'.")
                  dens <- object@sim$expDensity[[id]]
                  XL <- range(c(range(dens$x, na.rm=TRUE), obsFsir), na.rm=TRUE)
                  plot(dens, main=paste0("Individual: ", id, ", family: ",
                                         fam), xlab="Familial standardized incidence rate",
                       type="h", lwd=3, col="lightgrey", xlim=XL)
                  points(dens, col="grey", type="l", lwd=2)
              }
              Blue <- "#377EB8"
              abline(v=obsFsir, col=Blue)
              if(addLegend){
                  legend("topright",
                         legend=c(paste0("fsir     : ", format(obsFsir, digits=2)),
                                  paste0("p-value : ", format(object@sim$pvalue[id],
                                                              digits=2))
                                  ))
              }
          })

##****************************************************************************
##
##  plotPed
##
##  the plot ped can be done for individuals and for families. In essence we
##  are just calling the plotPed from the FAData object, eventually highlighting
##  the individual and providing the FSIR as numeric values below the id.
setMethod("plotPed", "FAStdIncidenceRateResults",
          function(object, id=NULL, family=NULL, filename=NULL,
                   device="plot", only.phenotyped=FALSE, ...){
              if(is.null(id) & is.null(family))
                  stop("Either id or family has to be provided!")
              if(length(object@sim) == 0){
                  ## Means we don't have any results...
                  FSIR <- rep(NA, length(object$id))
              }else{
                  FSIR <- object@sim$fsir
              }
              callNextMethod(object=object, id=id, family=family, filename=filename,
                             device=device, label1=FSIR, proband.id=id,
                             only.phenotyped=only.phenotyped, ...)
          })

##****************************************************************************
##
##  result
##
setMethod("result", "FAStdIncidenceRateResults", function(object, method="BH"){
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
                            fsir=NA,
                            pvalue=NA,
                            padj=NA,
                            check.names=FALSE, stringsAsFactors=FALSE)
        warning(paste0("No simulation data available! Please run a simulation ",
                       "using the 'runSimulation' method on this, or the ",
                       "'fsirTest' on a 'FAData' object."))
        return(MyRes)
    }
    nind <- length(object$id)
    MyRes <- data.frame(trait_name=rep(TraitName, nind),
                        total_phenotyped=rep(length(Trait), nind),
                        total_affected=rep(sum(Trait!=0), nind),
                        total_tested=rep(sum(!is.na(object@sim$fsir)), nind),
                        id=as.character(object$id),
                        family=as.character(object$family),
                        fsir=object@sim$fsir,
                        pvalue=object@sim$pvalue,
                        padj=p.adjust(object@sim$pvalue, method=method),
                        check.names=FALSE, stringsAsFactors=FALSE)
    MyRes <- MyRes[order(MyRes$pvalue, -MyRes$fsir), ]
    rownames(MyRes) <- as.character(MyRes$id)
    return(MyRes)
})

##****************************************************************************
##
##  fsir getter method.
##
setMethod("fsir", "FAStdIncidenceRateResults",
          function(object, ...){
              if(length(object@sim) == 0)
                  return(NULL)
              return(object@sim$fsir)
          })

##****************************************************************************
##
##  lambda getter method.
##
setMethod("lambda", "FAStdIncidenceRateResults",
          function(object, ...){
              return(object@lambda)
          })

##****************************************************************************
##
##  $ getter method.
##
setMethod("$", "FAStdIncidenceRateResults", function(x, name){
    if(name == "fsir" | name == "lambda" | name == "pvalue"
       | name == "timeInStrata"){
        ## extract local stuff
        if(length(x@sim) == 0){
            warning("No calculation result available!")
            return(NULL)
        }
        if(name == "fsir")
            return(fsir(x))
        if(name == "lambda")
            return(lambda(x))
        if(name == "pvalue")
            return(x@sim$pvalue)
        if(name == "timeInStrata")
            return(timeInStrata(x))
    }else{
        callNextMethod()
    }
})


##****************************************************************************
##
##  [ method.
##
setMethod("[", "FAStdIncidenceRateResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAStdIncidenceRateResults object is not supported!")
})

##****************************************************************************
##
##  timeInStrata getter/setter method.
##
setMethod("timeInStrata", "FAStdIncidenceRateResults",
          function(object){
              ## we're also checking that timeAtRisk, stored in its internal
              ## slot has the correct length!
              if(nrow(object@timeInStrata)==0)
                  return(object@timeInStrata)
              if(nrow(object@timeInStrata) != length(object$id))
                  stop("Length of 'timeInStrata' does not match the ",
                       "number of individuals in ",
                       "the pedigree!")
              return(object@timeInStrata)
          })
## setReplaceMethod("timeAtRisk", "FAIncidenceRateResults", function(object, value){
##     ## check if we've got the pedigree, otherwise does not make sense to store
##     ## the timeAtRisk.
##     if(length(object$id) == 0)
##         stop(paste0("No pedigree data in 'object'! Add the pedigree data ",
##                     "before setting 'timeAtRisk'."))
##     if(length(object$id) != length(value))
##         stop(paste0("Length of 'timeAtRisk' does not match the number of ",
##                     "individuals in the pedigree!"))
##     object@timeAtRisk <- value
##     return(object)
## })


##****************************************************************************
##
##  resultForId.
##
##  extract data for an individual.
setMethod("resultForId", "FAStdIncidenceRateResults",
          function(object, id=NULL){
              if(is.null(id))
                  stop("'id' has to be provided!")
              ## check if id is there.
              if(!any(object$id == id))
                  stop("Individual with id ", id, " not in the pedigree.")
              ## check if we have any results:
              if(length(object@sim) == 0)
                  stop("No simulation results available. Please run 'runSimulation' first.")
              idx <- which(object$id == id)
              return(list(id=id, fsir=object$fsir[id], pvalue=object@sim$pvalue[idx],
                          timeInStrata=object@timeInStrata[idx, ],
                          lambda=object@lambda))
          })

##****************************************************************************
##
##  trait<-
##
##  Calling the trait replecement method from FAResult and in addition
##  reset the simulation result.
setReplaceMethod("trait", "FAStdIncidenceRateResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})



