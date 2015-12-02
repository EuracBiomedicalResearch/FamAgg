## Methods for the classes
setMethod("cliques", "missing", function(...){
    return(igraph::cliques(...))
})

##************************************************
##
##       FAKinGroupResults
##
##
##************************************************
setMethod("show", "FAKinGroupResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    cat(paste0(" * Dimension of result data.frame: ",
               dim(result(object)),".\n"))
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})
setMethod("affectedKinshipGroups", "FAKinGroupResults", function(object){
    return(object@affectedKinshipGroups)
})
setMethod("independentGroupCount", "FAKinGroupResults", function(object, method="jo"){
    method <- match.arg(method, c("jo", "daniel"))
    if(length(object@affectedKinshipGroups) == 0)
        stop("No affected kinship groups present!")
    affectedKins <- object@affectedKinshipGroups
    indGroupCount <- 0
    ## on the long run we want to get rid of this...
    if(method=="jo"){
        AffIds <- lapply(affectedKins, function(z){return(z$aff)})
        for(i in 1:length(AffIds)){
            ## check: any of the ids of the group also in any of the other groups?
            if(length(intersect(unlist(AffIds[i]), unlist(AffIds[-i]))) == 0)
                indGroupCount <- indGroupCount + 1
        }
    }
    if(method=="daniel"){
        temp <- c()
        for(i in 1:length(affectedKins)) {
            if(length(intersect(temp, affectedKins[[i]]$aff)) == 0) {
                indGroupCount <- indGroupCount + 1
            }
            temp <- c(temp, affectedKins[[i]]$aff)
        }
    }
    return(indGroupCount)
})
## calling the trait replecement method from FAResult and in addition reset the simulation result.
setReplaceMethod("trait", "FAKinGroupResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@affectedKinshipGroups <- list()
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})


####
## The analysis method.
setMethod("runSimulation", "FAKinGroupResults", function(object, nsim=50000, strata=NULL){
    if(length(trait(object)) == 0)
        stop("No trait information available!")
    ## Eventually use the BiocParallel capabilities...
    ## ## check the BPPARAM:
    ## if(!is.null(BPPARAM) & !is(BPPARAM, "BiocParallelParam")){
    ##     warning("Submitted BPPARAM is not a BiocParallelParam object! Not performing parallel processing")
    ##     BPPARAM <- NULL
    ## }

    ## Get the kinship matrix: uses the pre-calculated if present or calculates one if none present.
    kin <- kinship(object)

    ## Get the ids of the affected and of the phenotyped.
    affectedIds <- affectedIndividuals(object)
    phenotypedIds <- phenotypedIndividuals(object)
    ## Ensure again that we do have pedigree information for these.
    affectedIds <- intersect(colnames(kin), affectedIds)
    phenotypedIds <- intersect(colnames(kin), phenotypedIds)

    ## break...
    if(length(affectedIds) <= 1){
        warning(paste0("Trait ",
                       object@traitname,
                       ": only one or none affected individual!"))
        return(object)
    }

    ## Strata. First check that strata has the correct length, then ensure that all of
    ## the affected and phenotyped ids have a non-NA value in strata.
    if(!is.null(strata)){
        if(length(strata) != length(object$id))
            stop("Argument 'strata' has to have the same length than there are individuals in the pedigree!")
        nas <- is.na(strata)
        names(strata) <- object$id
        newAffIds <- affectedIds[affectedIds %in% names(strata)[!nas]]
        if(length(newAffIds) != length(affectedIds)){
            warning(paste0("Removed ", length(affectedIds) - length(newAffIds),
                           " affected individuls as they have a missing value in 'strata'."))
            affectedIds <- newAffIds
        }
        newPheIds <- phenotypedIds[phenotypedIds %in% names(strata)[!nas]]
        if(length(newPheIds) != length(phenotypedIds)){
            warning(paste0("Removed ", length(phenotypedIds) - length(newPheIds),
                           " phenotyped individuls as they have a missing value in 'strata'."))
            phenotypedIds <- newPheIds
        }
        ## Subset strata and ensure same ordering than phenotypedIds
        strata <- strata[phenotypedIds]
    }

    ## Subset the kinship matrix
    ## Transform the Matrix into a matrix... needs more memory, but runs faster.
    affectedIds <- as.character(affectedIds)
    phenotypedIds <- as.character(phenotypedIds)
    kinAffected <- as.matrix(kin[affectedIds, affectedIds])
    kinPhenotyped <- as.matrix(kin[phenotypedIds, phenotypedIds])
    ## Remove self-self kinship
    diag(kinAffected) <- NA
    diag(kinPhenotyped) <- NA

    ## Get groups of affected related individuals around each affected individual
    ## affectedKins is then a list one element for each affected individual with
    ## $aff being all affected individuals having kinship with the individual and
    ## $phe all phenotyped individuals up to a kinship being smaller or equal to the
    ## smallest kinship of the individual with any other affected.
    ## $kinfreq: a table with the frequency (counts) of kinship values (smaller 0.5).
    ##           the names of the table correspond to the kinship values increasingly
    ##           ordered.
    affectedKins <- lapply(as.list(affectedIds), function(z){
        ## process all affected around the current affected
        currentKin <- kinAffected[, z]   ## that's a numeric now...
        ## get all individuals that share kinship...
        currentKin <- currentKin[which(currentKin > 0)]
        affIds <- c(z, names(currentKin))
        ## remove the self-self kinship
        ## currentKin <- currentKin[-(match(z, affIds))]
        minKin <- ifelse(length(currentKin) > 0, yes=min(currentKin, na.rm=TRUE), no=NA)
        meanKin <- ifelse(length(currentKin) > 0, yes=mean(currentKin, na.rm=TRUE), no=NA)
        ## process all phenotyped around the current affected
        currentKin <- kinPhenotyped[, z]
        ## Note: which works nicely here even if minMin is NA
        phenIds <- names(currentKin[which(currentKin >= minKin)])
        ## the kinship frequencies:
        ## returning a table with the count of all kinship values < 0.5; names of the table
        ## are the kinship values, increasingly ordered.
        if(length(affIds) >= 2){
            idx <- match(affIds, colnames(kinPhenotyped))
            temp <- kinPhenotyped[idx, idx]
            kinFreq <- table(as.vector(temp))
        }else{
            kinFreq <- NULL
        }
        return(list(aff=sort(affIds), pheno=sort(phenIds), kinfreq=kinFreq, meankin=meanKin,
                    minkin=minKin))
    })
    names(affectedKins) <- affectedIds

    ## remove duplicated entries...
    ## re-order affectedKins by phenotyped size... that way we select below the smaller kinship
    ## group if two groups have the same affected individuals. We select the smaller group to
    ## get more conservative p-values in the permutation below <- CHRIS
    affectedKins <- affectedKins[order(unlist(lapply(affectedKins, function(z){
        length(z$phe)
    })), decreasing=FALSE)]
    tmp <- lapply(affectedKins, function(z){z$aff})
    affectedKins <- affectedKins[match(unique(tmp), tmp)]  ## match returns only the first hit
    ## remove groups with a single affected
    affectedKins <- affectedKins[unlist(lapply(affectedKins, function(z){
        return(length(z$aff))
    })) > 1]
    ## re-order the affectedKins
    commonIds <- affectedIds[affectedIds %in% names(affectedKins)]
    affectedKins <- affectedKins[commonIds]

    ## break...
    if(length(affectedKins) == 0){
        warning(paste0("Trait ",
                       object@traitname,
                       ": No family with more than one affected member found!"))
        return(object)
    }

    ## Finally, do the simulations.
    ## 1) Sample number of affected individuals (length(affectedIds)) from all
    ##    phenotyped (phenotypedIds).
    ## 2) Loop over the affected kin groups (affectedKins)
    ##    2.1) get those phenotyped in the affected kinship group that match the sampled ids.
    ##    2.2) if the expected overlap is >= than the observed -> +1 to frequency ratio.
    ##    2.3) get the kinship matrix for the overlap (if overlap >=2)
    ##    2.4) make a frequency table of the expected kinship values and check:
    ##         a) if the largest expected kinship value is > any of the observed kinship values: -> +1
    ##            to the kinship frequency.
    ##         b) if the largest expected kinship value == the largest observed kinship value: compare
    ##            their counts and if the expected count >= the observed count: -> +1 to the kinship frequency.
    ## In other words:
    ## Define random sets of affected individuals by randomly selecting X individuals among all phenotyped
    ## individuals, with X being the total number of affected cases. Then, for each group, iteratively determine
    ## the number of individuals from the group which have been labeled as affected in the simulation and their
    ## kinship value.
    ExpRatio <- rep(0, length(affectedKins))
    ExpKin <- rep(0, length(affectedKins))
    names(ExpRatio) <- names(affectedKins)
    names(ExpKin) <- names(affectedKins)
    ## parallelize the simulations?
    ## in case: generate a list with the sampled affected ids and run a mclapply over than.
    ## the speed gain might however not be that large... so we don't do it yet.
    ## for parallelization its better to perform the sampling once!
    if(is.null(strata)){
        Sims <- lapply(1:nsim, function(z){
            return(sample(phenotypedIds, length(affectedIds)))
        })
    }else{
        affStrataCounts <- table(strata[affectedIds])
        affStrataCounts <- affStrataCounts[affStrataCounts > 0]
        Sims <- lapply(1:nsim, function(z){
            idx <- stratsample(strata, counts=affStrataCounts)
            return(phenotypedIds[idx])
        })
    }
    ##        BigRes <- bplapply(Sims, BPPARAM=BPPARAM, FUN=function(x){
    BigRes <- lapply(Sims, FUN=function(x){
        ExpRes <- do.call(rbind, lapply(affectedKins, function(z){
            Res <- c(0, 0, 0)
            expectedAffInKin <- intersect(x, z$phe)
            Res[3] <- length(expectedAffInKin)
            if(length(expectedAffInKin) >= length(z$aff)){
                Res[1] <- 1
            }
            ## check frequency...
            if(length(expectedAffInKin) >= 2){
                ## that's the slowest part of the function... subsetting of the kinship matrix...
                ## we get some speedup if we subset a matrix, not a Matrix!
                expectedAffInKin <- match(expectedAffInKin, colnames(kinPhenotyped))
                expectedKin <- kinPhenotyped[expectedAffInKin, expectedAffInKin]
                diag(expectedKin) <- NA
                ## get the table with the frequencies of kinship values...
                expectedKinFreqs <- table(as.vector(expectedKin))
                observedKinFreqs <- z$kinfreq
                if(length(expectedKinFreqs) > 0 & length(observedKinFreqs) > 0){
                    ## names of the table are the kinship values, ordered increasingly by size.
                    expKinShipValue <- names(expectedKinFreqs[length(expectedKinFreqs)])
                    obsKinShipValue <- names(observedKinFreqs[length(observedKinFreqs)])
                    if(as.numeric(expKinShipValue) >= as.numeric(obsKinShipValue)){
                        ## if the kinship value is the same, compare the counts.
                        if(expKinShipValue==obsKinShipValue){
                            if(expectedKinFreqs[expKinShipValue] >= observedKinFreqs[obsKinShipValue])
                                Res[2] <- 1
                        }else{
                            ## means its larger, so add +1
                            Res[2] <- 1
                        }
                    }
                }
            }
            return(Res)
        }))
        return(ExpRes)
    })
    ## BigRes is now a list, length nsim, each element being a matrix, ncol=2, nrow=length affectedKins.
    ExpRatio <- colSums(do.call(rbind, lapply(BigRes, function(simRes){simRes[, 1]})))
    ExpKin <- colSums(do.call(rbind, lapply(BigRes, function(simRes){simRes[, 2]})))
    ExpVals <- do.call(rbind, lapply(BigRes, function(simRes){simRes[, 3]}))
    ExpVals <- split(t(ExpVals), f=names(affectedKins))
    expDensity <- lapply(ExpVals, density)
    expHist <- lapply(ExpVals, hist, plot=FALSE)
    expDensity <- expDensity[names(affectedKins)]
    expHist <- expHist[names(affectedKins)]
    names(ExpRatio) <- names(affectedKins)
    names(ExpKin) <- names(affectedKins)
    ExpRatio <- ExpRatio / nsim
    ExpKin <- ExpKin / nsim

    object@affectedKinshipGroups <- affectedKins
    object@sim <- list(pvalueKinship=ExpKin, pvalueRatio=ExpRatio,
                       expDensity=expDensity, expHist=expHist)
    object@nsim <- nsim

    return(object)
})



setMethod("result", "FAKinGroupResults", function(object, method="BH"){
    ## generate the result table...
    method <- match.arg(method, p.adjust.methods)
    Trait <- trait(object, na.rm=TRUE)
    TraitName <- object@traitname
    if(length(TraitName)==0)
        TraitName <- NA
    if(length(object@sim)==0){
        ## generate a dummy matrix...
        MyRes <- data.frame(trait_name=TraitName,
                            total_phenotyped=length(Trait),
                            total_affected=sum(Trait!=0),
                            phenotyped=0,
                            affected=0,
                            group_id=NA,
                            family=NA,
                            group_phenotyped=NA,
                            group_affected=NA,
                            ratio_pvalue=NA,
                            ratio_padj=NA,
                            mean_kinship=NA,
                            kinship_pvalue=NA,
                            kinship_padj=NA,
                            check.names=FALSE, stringsAsFactors=FALSE)
        warning(paste0("No simulation data available! This means that",
                       " either no simulation was run yet (using the ",
                       "kinshipTest function or runSimulation method) ",
                       "or that the simulation returned no result."))
        return(MyRes)
    }
    affKin <- affectedKinshipGroups(object)
    ped <- pedigree(object)
    ped <- ped[ped$id %in% names(affKin), c("id", "family")]
    pedL <- split(ped, f=ped$id)
    fams <- unlist(lapply(pedL, function(z){
        return(paste(unique(z[, "family"]), collapse=", "))
    }))
    fams <- fams[names(affKin)]
    MyRes <- data.frame(trait_name=rep(TraitName, length(affKin)),
                        total_phenotyped=rep(length(Trait), length(affKin)),
                        total_affected=rep(sum(Trait!=0), length(affKin)),
                        phenotyped=rep(
                            length(unique(unlist(lapply(affKin, function(z)z$phe)))),
                            length(affKin)),
                        affected=rep(
                            length(unique(unlist(lapply(affKin, function(z)z$aff)))),
                            length(affKin)),
                        group_id=names(affKin),
                        family=fams,
                        group_phenotyped=unlist(lapply(affKin, function(z){length(z$phe)})),
                        group_affected=unlist(lapply(affKin, function(z){length(z$aff)})),
                        ratio_pvalue=object@sim$pvalueRatio,
                        ratio_padj=p.adjust(object@sim$pvalueRatio, method=method),
                        mean_kinship=unlist(lapply(affKin, function(z)z$meankin)),
                        kinship_pvalue=object@sim$pvalueKinship,
                        kinship_padj=p.adjust(object@sim$pvalueKinship, method=method),
                        check.names=FALSE, stringsAsFactors=FALSE
                        )
    MyRes <- MyRes[order(MyRes$ratio_pvalue, MyRes$kinship_pvalue, -MyRes$mean_kinship), ]
    rownames(MyRes) <- as.character(MyRes$group_id)
    return(MyRes)
})


#############
## plotting method...
## plotPed representing the results from the kinship test.
## plotPed for FAKinGroupResults does only support plotting for id.
setMethod("plotPed", "FAKinGroupResults",
          function(object, id=NULL,
                   family=NULL, filename=NULL,
                   device="plot", ...){
                       if(!is.null(family))
                           stop("Generating a pedigree for a family is not supported for FAKinGroupResults. See help for more information.")
                       affGroup <- affectedKinshipGroups(object)
                       if(!any(names(affGroup) == id))
                           stop("The id should be one of the names of the affected kinship groups (i.e. names(affectedKinshipGroups(object)))!")
                       affGroup <- affGroup[[id]]
                       callNextMethod(object=object, id=id, family=family, filename=filename,
                                      device=device, proband.id=unique(c(affGroup$aff, affGroup$pheno)), ...)
                   })


## this buildPed simply ensures that all phenotyped in the group will be included in the pedigree.
setMethod("buildPed", "FAKinGroupResults",
          function(object, id=NULL, max.generations.up=3, max.generations.down=16,
                   prune=FALSE){
              if(is.null(id))
                  stop("The id of the group has to be speficied!")
              affGroup <- affectedKinshipGroups(object)
              if(!any(names(affGroup) == id))
                  stop("The id should be one of the names of the affected kinship groups (i.e. names(affectedKinshipGroups(object)))!")
              affGroup <- affGroup[[id]]
              phenoAndAff <- as.character(unique(c(affGroup$aff, affGroup$pheno)))
              kin <- kinship(object)
              ped <- family(object, id=id)
              whichs <- which(colnames(kin) %in% ped$id)
              ## get the kinship of the kin with all other in the pedigree
              kin <- kin[id, whichs]
              keep <- names(kin)[kin >= affGroup$minkin]
              keep <- unique(c(keep, affGroup$pheno, affGroup$aff))
              ped <- callNextMethod(object=object, id=keep, prune=prune,
                                    max.generations.up=max.generations.up,
                                    max.generations.down=max.generations.down)
              return(ped)
          })

## this is to get all those that are related with any individuals of the group!!!
setMethod("shareKinship", "FAKinGroupResults",
          function(object, id=NULL){
              if(is.null(id))
                  stop("The id of the group has to be speficied!")
              affGroup <- affectedKinshipGroups(object)
              if(!any(names(affGroup) == id))
                  stop("The id should be one of the names of the affected kinship groups (i.e. names(affectedKinshipGroups(object)))!")
              affGroup <- affGroup[[id]]
              allids <- unique(c(affGroup$aff, affGroup$pheno))
              ## get all individuals of the group...
              return(doShareKinship(kin=kinship(object), id=allids))
          })



## plot the pedigree for the kinship group...
## the function does add eventually missing parents.
plotPedForKinshipResult <- function(object, groupid=NULL, filename=NULL,
                                    device="pdf", ...){
    if(is.null(groupid))
        stop("No id specified!")
    if(length(groupid) > 1){
        groupid <- groupid[1]
        warning("length groupid > 1; using only the first element!")
    }
    groupid <- as.character(groupid)
    if(!any(names(affectedKinshipGroups(object)) == groupid))
        stop("groupid should be one of the names of the affected kinship groups (i.e. names(affectedKinshipGroups(object)))!")
    affs <- affectedKinshipGroups(object)[[groupid]]$aff
    phens <- affectedKinshipGroups(object)[[groupid]]$pheno
    fam <- buildPedigree(family(object, id=groupid), ids=unique(c(affs, phens)))
    cat("nrow fam ", nrow(fam), " length affs phens: ", length(unique(c(affs, phens))), "\n")
    ## check if we've got affected
    if(any(colnames(fam)=="affected")){
        affected <- fam[, "affected"]
    }else{
        affected <- rep(NA, nrow(fam))
        names(affected) <- fam$id
    }
    ## check ages...
    ages <- age(object)[fam$id]
    ages[!is.na(ages)] <- paste(format(ages[!is.na(ages)], digits=3), "years")
    ## ages[!is.na(ages)] <- sprintf("%d years", ages[!is.na(ages)])
    ages[is.na(ages)] <- ""
    ## do not have affected... obviously...
    is.proband <- rep(FALSE, nrow(fam))
    names(is.proband) <- fam$id
    if(!is.null(groupid))
        is.proband[as.character(groupid)] <- TRUE
    ## how is the code from Daniel labeling the affected:
    ## - any phenotyped: 1
    ## - any affected: 4
    ## - any affected in the affected group: 2.
    ## OK, now plot!!!
    res <- doPlotPed(family=fam$family, individual=fam$id, father=fam$father,
                     mother=fam$mother, gender=fam$sex, is.proband=is.proband,
                     text1.below.symbol=ages, affected=affected,
                     filename=filename, device=device, ...)
    res
}

setMethod("[", "FAKinGroupResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAKinGroupResults object is not supported!")
})


## This will be a crazy funky method to plot the simulation results.
setMethod("plotRes", "FAKinGroupResults",
          function(object, id=NULL,
                   family=NULL, addLegend=TRUE, type="density", ...){
              type <- match.arg(type, c("density", "hist"))
              if(length(object@sim) == 0)
                  stop("No analysis performed yet!")
              if(is.null(id))
                  stop("The id of the affected individual for whom the simulation results should be displayed has to be specified!")
              id <- as.character(id)
              ## result on family is not allowed.
              if(!is.null(family) | length(id) > 1)
                  stop("plotRes for FAKinGroupResults does only support specifying a single id!")
              ## check if the id is an affected for which we tested...
              if(!any(names(object@sim$pvalueKinship) == id))
                  stop(paste0("No simulation result is available for individual ", id,
                              ". id should be one of result(object)$group_id."))
              ## OK, now we can really start doing things...
              affGroup <- affectedKinshipGroups(object)[[id]]
              meanKin <- affGroup$meankin
              affCount <- length(affGroup$aff)
              pheCount <- length(affGroup$phe)
              kinPval <- object@sim$pvalueKinship[id]
              fam <- family(object, id=id)[1, "family"]
              ratioP <- object@sim$pvalueRatio[id]
              kinP <- object@sim$pvalueKinship[id]
              par(xpd=FALSE)
              if(type == "density"){
                  toplot <- object@sim$expDensity[[id]]
                  XL <- range(c(range(toplot$x), affCount))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Affected count", type="h",
                       lwd=3, col="lightgrey", xlim=XL)
                  points(toplot, col="grey", type="l", lwd=2)
              }
              if(type == "hist"){
                  toplot <- object@sim$expHist[[id]]
                  XL <- range(c(range(toplot$mids), affCount))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Affected count",
                       col="lightgrey", border="grey", xlim=XL)
              }
              Blue <- "#377EB8"
              abline(v=affCount, col=Blue)
              if(addLegend){
                  legend("topright", legend=c(paste0("affect. count: ", affCount),
                                              paste0("pheno. count:  ", pheCount),
                                              paste0("ratio p-value: ", format(ratioP, digits=3)),
                                              paste0("mean kinship:  ", format(meanKin, digits=3)),
                                              paste0("kinhip p-value: ", format(kinPval, digits=3))
                                             ))
              }
})



