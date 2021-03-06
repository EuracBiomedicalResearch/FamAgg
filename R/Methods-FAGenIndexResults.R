## These are methods for the "genealogical index" for familial clustering [Hill, 1980]
##************************************************
##
##       FAGenIndexResult
##       Methods for result objects from the genealogical index test.
##
##************************************************
setMethod("show", "FAGenIndexResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    cat(paste0(" * Dimension of result data.frame: ",
               length(object@sim),".\n"))
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
    cat(paste0(" * Method for control set definition: ",
               object@controlSetMethod, ".\n"))
    if(object@perFamilyTest){
        cat(" * Analysis performed for each family.\n")
    }else{
        cat(" * Analysis performed for the whole pedigree.\n")
    }
})
## calling the trait replecement method from FAResult and in addition
## reset the simulation result.
setReplaceMethod("trait", "FAGenIndexResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    object@controlSetMethod <- "getAll"
    object@perFamilyTest <- FALSE
    return(object)
})

####
## the analysis method.
setMethod("runSimulation", "FAGenIndexResults",
          function(object, nsim = 50000, perFamilyTest = FALSE,
                   controlSetMethod = "getAll", rm.singletons = TRUE,
                   strata = NULL, ...) {
              if(length(trait(object)) == 0)
                  stop("No trait information available!")
              kin <- kinship(object)
              trait <- trait(object)
              ResList <- .genIndex(ped = pedigree(object),
                                   kin = kinship(object),
                                   trait = trait(object),
                                   perFamilyTest = perFamilyTest,
                                   controlSetMethod = controlSetMethod,
                                   nsim = nsim,
                                   prune = rm.singletons, strata = strata, ...)
              object@nsim <- nsim
              object@controlSetMethod <- controlSetMethod
              object@perFamilyTest <- perFamilyTest
              object@sim <- ResList
              object
          })

## results; get the results table.
setMethod("result", "FAGenIndexResults", function(object, method="BH"){
    method <- match.arg(method, p.adjust.methods)
    TraitName <- object@traitname
    if(length(TraitName)==0)
        TraitName <- NA
    if(length(object@sim) == 0){
        ## generate a dummy matrix...
        MyRes <- data.frame(trait_name=TraitName,
                            total_phenotyped=length(phenotypedIndividuals(object)),
                            total_affected=length(affectedIndividuals(object)),
                            entity_id=NA,
                            entity_ctrls=NA,
                            entity_affected=NA,
                            genealogical_index=NA,
                            pvalue=NA,
                            padj=NA,
                            check.names=FALSE, stringsAsFactors=FALSE
                            )
        warning("No simulation data available! This means that either no ",
                "simulation was run yet (using the genealogicalIndex function ",
                "or runSimulation method) or that the simulation returned no ",
                "result (i.e. too few affected individuals in the trait).")
        return(MyRes)
    }
    ## otherwise compile the result...
    MyRes <- do.call(rbind, lapply(object@sim, function(z){
        return(c(entity_ctrls=length(z$ctrls),
                 entity_affected=length(z$affected),
                 genealogical_index=z$meanKinship * 100000,
                 pvalue=z$pvalueKinship))
    }))
    MyRes <- data.frame(trait_name=rep(TraitName, nrow(MyRes)),
                        total_phenotyped=rep(length(phenotypedIndividuals(object)),
                                             nrow(MyRes)),
                        total_affected=rep(length(affectedIndividuals(object)),
                                           nrow(MyRes)),
                        entity_id=as.character(names(object@sim)),
                        MyRes,
                        padj=p.adjust(MyRes[, "pvalue"], method=method),
                        stringsAsFactors=FALSE, check.names=FALSE)
    MyRes <- MyRes[order(MyRes$pvalue), ]
    rownames(MyRes) <- as.character(MyRes$entity_id)
    MyRes
})


## approach:
## * per family: do the analysis family-wise: find families with.
## * "normal" analysis: whole pedigree: disable getGenerationMatched!
## prepare: prune the pedigree/family, i.e. remove unconnected individuals.
## 1) select affected, calculate mean kinship.
## 2) select matched controls: using a control selection method.
##    sample from them, using eventually a strata sampler.
## ADDON: use a stratafun? something like strata by generation, strata by generation
## and sex or alike.
##
## ped: data.frame with pedigree information (pedigree(fad), fad being a FAData)
## kin: kinship matrix (kinship(fad))
## trait: numeric or logical vector, length=nrow(ped)
## perFamilyTest: perform the test for each family separately. That way we determine
##                whether the kinship between affected in that family is higher than
##                we might expect by chance, given the kinship between individuals in
##                the same family. Thus, we test for the presence of a cluster of close
##                relatives in the family affected by the trait.
## controlSetMethod: function to select matched controls.
## nsim: number of simulations.
## prune: remove un-connected individuals from the pedigree; should always be TRUE
## strata: to do strata sampling.
## aggFun: the function to be used to aggregate the kinship of affected. Defaults to
##         mean, as in the original definition by Hill, but a kinship sum might be
##         more informative.
##
## return:
## a list with:
## + meanKinship: the mean kinship between affected.
## + pvalueKinship: the p-value.
## + expDensity: the density distribution of expected mean kinships of random samples.
## + affected: the ids of the affected.
## + ctrls: the ids of the (matched) controls.
.genIndex <- function(ped=NULL, kin=NULL, trait=NULL, perFamilyTest=FALSE,
                      controlSetMethod="getAll", nsim=50000, prune=TRUE,
                      strata=NULL, aggFun=mean, ...){
    controlSetMethod <- match.arg(controlSetMethod, c("getAll",
                                                      "getSexMatched",
                                                      "getGenerationMatched",
                                                      "getGenerationSexMatched",
                                                      "getExternalMatched"))
    if(!perFamilyTest){
        ## don't allow some matched control methods.
        if((controlSetMethod == "getGenerationMatched") |
           (controlSetMethod == "getGenerationSexMatched")){
            stop("Can not use 'getGenerationMatched' or ",
                 "'getGenerationSexMatched' functions for ",
                 "'perFamilyTest=FALSE' since the generation estimation ",
                 "procedure only allows to define within-family generations, ",
                 "not generations across families in a full pedigree.")
        }
        entityIs <- "pedigree"
        ## replace the family id with 1 for all... thus performing the test for
        ## the whole pedigree in one go.
        ped[, "family"] <- 1
    }else{
        entityIs <- "family"
    }
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(trait))
        stop("No trait information submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    ## add the trait information to the pedigree data.frame
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there are",
             " rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, AFF=trait)
    ## after all we shouldn't need this; just in case. it's for prune=TRUE to
    ## add eventually removed mates of a childless couple; but in the end we
    ## don't want them anyway, and they are NOT defined in a pedigree.
    dotList <- list(...)
    if(any(names(dotList) == "addMissingMates")){
        addMisMate <- dotList$addMissingMates
    }else{
        addMisMate <- FALSE
    }
    ## check strata
    haveStrata <- FALSE
    if(!is.null(strata)){
        if(length(strata)!=nrow(ped))
            stop("Argument 'strata' has to have the same length than there",
                 " are rows (individuals) in the pedigree 'ped'!")
        ## add strata.
        ped <- cbind(ped, STRATA=strata)
        haveStrata <- TRUE
    }
    ## First remove singletons
    if(prune){
        ped <- removeSingletons(ped)
    }
    ## OK, now we can go on: lapply on the splitted ped by family
    pedL <- split(ped, ped$family)
    Res <- lapply(pedL, function(z){
        ## Now get the controls and the affected ids. We have to eventually
        ## further subset these if we have non-phenotyped controls etc.
        affIds <- as.character(z[which(z$AFF > 0), "id"])
        ## Eventually remove affected without a value in strata:
        message("Cleaning data set (got in total ", nrow(z), " individuals):")
        if(haveStrata){
            message(" * affected individuals without valid strata values...",
                    appendLF=FALSE)
            nas <- is.na(z[affIds, "STRATA"])
            if(any(nas)){
                affIds <- affIds[!nas]
                message(" ", sum(nas), " removed.")
            }else{
                message(" none present.")
            }
        }
        ## Return if there are not enough affected.
        if(length(affIds) < 2){
            warning(paste0("Can not perform test for ",
                         entityIs, " ", z[1, "family"],
                         ": less than 2 affected individuals in the trait."))
            return(list(meanKinship=NA,
                        pvalueKinship=NA,
                        expDensity=NULL,
                        affected=affIds,
                        ctrls=NA))
        }
        ## Get the matched controls for these guys...
        ctrls <- as.character(
            unlist(do.call(what=controlSetMethod, args=list(object=z,
                                                            id=affIds, ...)),
                   use.names=FALSE)
            )
        ## Subset to phenotyped individuals, i.e. those with a valid value in
        ## "affected"
        nas <- is.na(z[ctrls, "affected"])
        message(" * not phenotyped individuals among selected controls...",
                appendLF=FALSE)
        if(any(nas)){
            ctrls <- ctrls[!nas]
            message(" ", sum(nas), " removed.")
        }else{
            message(" none present.")
        }
        if(haveStrata){
            ## Subsetting the data.frame to elements that have a non-NA value
            ## in strata.
            nas <- is.na(z[ctrls, "STRATA"])
            message(" * control individuals without valid strata values...",
                    appendLF=FALSE)
            if(any(nas)){
                ctrls <- ctrls[!nas]
                message(" ", sum(nas), " removed.")
            }else{
                message(" none present.")
            }
        }
        message("Done")
        if(length(ctrls) < 4*length(affIds))
            warning("The proportion of affected individuals among the ",
                    "matched controls is > 25% in ", entityIs, " ",
                    z[1, "family"], ".")
        ## Subset the kinship matrix and transform into a "normal" matrix,
        ## assuming that subsetting and other functions run faster on that.
        bm <- colnames(kin) %in% unique(c(affIds, ctrls))
        localKin <- as.matrix(kin[bm, bm])
        ## Remove self-self kinships; we don't want to count them.
        diag(localKin) <- NA
        ## Calculate the (observed) mean kinship between the real affected.
        localKinColNames <- colnames(localKin)
        bm <- localKinColNames %in% affIds
        obsMeanKin <- aggFun(as.vector(localKin[bm, bm]), na.rm=TRUE)
        ## OK, now that we've got the controls, do the test!
        if(!haveStrata){
            affSize <- length(affIds)
            expMeanKins <- sapply(1:nsim, function(y){
                ## pick random samples, same size than affIds
                expIds <- sample(ctrls, size=affSize)
                bm <- localKinColNames %in% expIds
                return(aggFun(as.vector(localKin[bm, bm]), na.rm=TRUE))
            })
        }else{
            affStrataCounts <- table(z[affIds, "STRATA"])
            affStrataCounts <- affStrataCounts[affStrataCounts > 0]
            phenoStrata <- z[ctrls, "STRATA"]
            expMeanKins <- sapply(1:nsim, function(y){
                ## pick random strata sample using the stratsample from the
                ## survey package
                expIdx <- stratsample(phenoStrata, affStrataCounts)
                bm <- localKinColNames %in% ctrls[expIdx]
                return(aggFun(as.vector(localKin[bm, bm]), na.rm=TRUE))
            })
        }
        return(list(meanKinship=obsMeanKin,
                    pvalueKinship=sum(expMeanKins >= obsMeanKin)/nsim,
                    expDensity=density(expMeanKins),
                    expHist=hist(expMeanKins, plot=FALSE),
                    affected=affIds, ctrls=ctrls))
    })
    Res
}

#############
## plotting method...
## plotPed representing the results from the kinship clustering test.
## plotPed for FAKinClustResult does only support plotting for id.
setMethod("plotPed", "FAGenIndexResults",
          function(object, id=NULL,
                   family=NULL, filename=NULL,
                   device="plot", ...){
              if(!is.null(id)){
                  ## get the family id for this id...
                  if(!any(as.character(object$id) == id))
                      stop("Can not find any individual with id ", id,
                           " in the pedigree!")
                  family <- as.character(object$family)[as.character(object$id) == id]
                  family <- family[1]
              }
              if(length(object@sim) > 0){
                  if(!object@perFamilyTest){
                      warning("The genealogical index test has been performed ",
                              "on the full pedigree, but the plot was only ",
                              "generated for family ", family, "!")
                      ctrls <- object@sim[[1]]$ctrls
                  }else{
                      ctrls <- object@sim[[family]]$ctrls
                  }
              }else{
                  ctrls <- NA
              }

              Blue <- "#377EB8"
              callNextMethod(object = object, id = id, family = family,
                             filename = filename, device = device,
                             proband.id = ctrls, ...)
              ## alternatively use highlight.ids
          })


setMethod("[", "FAGenIndexResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAGenIndexResults object is not supported!")
})

## This will be a crazy funky method to plot the simulation results.
setMethod("plotRes", "FAGenIndexResults",
          function(object, id=NULL,
                   family=NULL, addLegend=TRUE, type="density", ...){
              type <- match.arg(type, c("density", "hist"))
              if(length(object@sim) == 0)
                  stop("No analysis performed yet!")
              if(!object@perFamilyTest){
                  ## Just ignore id and family...
                  meanKin <- object@sim[[1]]$meanKinship
                  pval <- object@sim[[1]]$pvalueKinship
                  dens <- object@sim[[1]]$expDensity
                  hist <- object@sim[[1]]$expHist
                  affCount <- length(object@sim[[1]]$affected)
                  ctrlsCount <- length(object@sim[[1]]$ctrls)
                  family <- 1
              }else{
                  if(!is.null(id)){
                      bm <- as.character(object$id) == id
                      if(!any(bm))
                          stop("No individual with id ", id,
                               " present in the pedigree!")
                      family <- as.character(object$family)[bm][1]
                  }
                  if(is.null(family))
                      stop("With per family genealogical index ",
                           "(perFamilyTest=TRUE) either the id of an ",
                           "individual or the id of a family has to be ",
                           "specified!")
                  if(!any(names(object@sim) == family))
                      stop("No family with id ", family,
                           " present in the pedigree!")
                  meanKin <- object@sim[[as.character(family)]]$meanKinship
                  pval <- object@sim[[as.character(family)]]$pvalueKinship
                  dens <- object@sim[[as.character(family)]]$expDensity
                  hist <- object@sim[[as.character(family)]]$expHist
                  affCount <- length(object@sim[[as.character(family)]]$affected)
                  ctrlsCount <- length(object@sim[[as.character(family)]]$ctrls)
                  if(is.na(meanKin))
                      stop("No results for family ", family, " available!")
              }

              ## plot it...
              par(xpd=FALSE)
              entity <- ifelse(object@perFamilyTest, yes="Family", no="Pedigree")
              if(type == "density"){
                  plot(dens, main = paste0(entity, " ", family, "; ",
                                           object@traitname),
                       xlab="Mean kinship", type="h", lwd=3, col="lightgrey",
                       xlim=range(c(range(dens$x), meanKin)))
                  points(dens, col="grey", type="l", lwd=2)
              }
              if(type == "hist"){
                  plot(hist, main=paste0(entity, " ", family, "; ",
                                         object@traitname), xlab="Mean kinship",
                       col="lightgrey", border="grey",
                       xlim=range(c(range(hist$x), meanKin)))
              }
              Blue <- "#377EB8"
              abline(v=meanKin, col=Blue)
              if(addLegend){
                  legend(
                      "topright",
                      legend=c(paste0("GIF: ", format(100000*meanKin, digits=4)),
                               paste0("mean kinship: ", format(meanKin, digits=2)),
                               paste0("p-value     : ", format(pval, digits=3)),
                               paste0("affect count: ", affCount),
                               paste0("ctrls count : ", ctrlsCount)))
              }
          })

