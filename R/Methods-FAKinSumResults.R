##************************************************
##
##       FAKinSumResults
##       Most of the code here is from Chris X. Weichenberger
##
##************************************************
setMethod("show", "FAKinSumResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    cat(paste0(" * Dimension of result data.frame: ",
               dim(result(object)),".\n"))
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})
## calling the trait replecement method from FAResult and in addition
## reset the simulation result.
setReplaceMethod("trait", "FAKinSumResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})

####
## the analysis method.
setMethod("runSimulation", "FAKinSumResults",
          function(object, nsim=50000, strata=NULL, ...){
              if(length(trait(object)) == 0)
                  stop("No trait information available!")
              ind <- TRUE
              kin <- kinship(object)
              ## check strata
              ped <- pedigree(object)
              affectedIds <- as.character(ped[which(ped$affected > 0), "id"])
              phenotypedIds <-as.character(ped[!is.na(ped$affected), "id"])
              message("Cleaning data set (got in total ", nrow(ped), " individuals):")
              ## Dummy for removing unaffected:
              message(" * not phenotyped individuals...")
              if(length(phenotypedIds) != nrow(ped)){
                  message(" ", nrow(ped)-length(phenotypedIds), " removed")
              }else{
                  message(" none present.")
              }
              if(!is.null(strata)){
                  if(length(strata)!=nrow(ped))
                      stop("Argument 'strata' has to have the same length than there are individuals in the pedigree!")
                  ## add strata.
                  ped <- cbind(ped, STRATA=strata)
                  ## Subset the affectedIds and phenotypedIds to those for which strata
                  ## is not NA.
                  ped <- ped[!is.na(ped$STRATA) & !is.na(ped$affected), , drop=FALSE]
                  strata <- as.character(ped$STRATA)
                  names(strata) <- as.character(ped$id)
                  newAffIds <- as.character(ped[which(ped$affected > 0), "id"])
                  newPhenoIds <- as.character(ped[, "id"])
                  message(" * unaffected individuals without valid strata values...", appendLF=FALSE)
                  if(length(newPhenoIds)!=length(phenotypedIds)){
                      ## warning(paste0("Removed ", length(phenotypedIds) - length(newPhenoIds),
                      ##                " phenotyped individuals, as they have",
                      ##                " a missing value in strata."))
                      message(" ", length(phenotypedIds)-length(newPhenoIds), " removed.")
                      phenotypedIds <- newPhenoIds
                  }else{
                      message(" none present.")
                  }
                  message(" * affected individuals without valid strata values...", appendLF=FALSE)
                  if(length(newAffIds)!=length(affectedIds)){
                      ## warning(paste0("Removed ", length(affectedIds) - length(newAffIds),
                      ##                " affected individuals, as they have",
                      ##                " a missing value in strata."))
                      message(" ", length(affectedIds)-length(newAffIds), " removed.")
                      affectedIds <- newAffIds
                  }else{
                      message(" none present.")
                  }
              }
              message("Done")
              ## now forcing the kinship matrix to match ordering of phenotypedIds and
              ## eventually strata.
              kin <- kin[phenotypedIds, phenotypedIds]
              ## removing self-self kinship:
              diag(kin) <- NA
              object@nsim <- nsim
              ## break...
              if(length(affectedIds) <= 1){
                  warning(paste0("Trait ",
                                 object@traitname,
                                 ": only one or none affected individual!"))
                  return(object)
              }
              if(ind){
                  res <- significant.individuals.althyp(kin, affectedIds,
                                                        phenotypedIds,
                                                        nr.sim=nsim,
                                                        strata=strata,
                                                        ...)
              }else{
                  ## res <- significan.population(kin, affectedIds, phenotypedIds,
                  ##                              nr.sim=nsim)
              }
              object@sim <- res
              return(object)
          })

## results; get the results table.
setMethod("result", "FAKinSumResults", function(object, method="BH"){
    method <- match.arg(method, p.adjust.methods)
    TraitName <- object@traitname
    if(length(TraitName)==0)
        TraitName <- NA
    if(length(object@sim) == 0){
        ## generate a dummy matrix...
        MyRes <- data.frame(trait_name=TraitName,
                            total_phenotyped=length(phenotypedIndividuals(object)),
                            total_affected=length(affectedIndividuals(object)),
                            affected_id=NA,
                            family=NA,
                            affected=NA,
                            kinship_sum=NA,
                            freq=NA,
                            pvalue=NA,
                            padj=NA,
                            check.names=FALSE, stringsAsFactors=FALSE)
        warning("No simulation data available! This means that either no simulation was",
                " run yet (using the kinshipClusterTest function or runSimulation method)",
                " or that the simulation returned no result (i.e. too few affected",
                " individuals in the trait).")
        return(MyRes)
    }
    resList <- object@sim
    affIds <- as.character(names(resList$sumKinship))
    ped <- pedigree(object)
    ped <- ped[ped$id %in% affIds, c("id", "family")]
    pedL <- split(ped, f=ped$id)
    fams <- unlist(lapply(pedL, function(z){
        return(paste(unique(z[, "family"]), collapse=", "))
    }))
    fams <- fams[affIds]
    res <- data.frame(
        trait_name=TraitName,
        total_phenotyped=length(phenotypedIndividuals(object)),
        total_affected=length(affectedIndividuals(object)),
        affected_id=affIds,
        family=fams,
        affected=rep(length(resList$affected), length(resList$sumKinship)),
        kinship_sum=resList$sumKinship,
        freq=resList$frequencyKinship,
        pvalue=resList$pvalueKinship,
        padj=p.adjust(resList$pvalueKinship, method=method),
        stringsAsFactors=FALSE)
    rownames(res) <- as.character(res$affected_id)
    ## colnames(res) <- c("trait_name", "total_phenotyped", "total_affected", "affected_id",
    ##                    "kinship_sum", "freq", "pvalue", "padj")
    ## warning("Don't know yet how the Pedigree Size is calculated!!!")
    return(res[order(res$pvalue, -res$kinship_sum), , drop=FALSE])
})



## null hypothesis. H0: sum of kinship values of affected
## with all other affected is random. Test: is the sum of kinship values of
## affected with all other affected larger or equal that we would expect by
## chance (randomly selected individuals being affected)? The result is
## essentially identical to the significant.individuals2 function.
## NOTE: this function is robust against cases in which the observed kinship
## sum was never sampled in the simulation! The significant.individuals.2
## function would return a NA for the p-value and the frequency, while the
## p-value could be non-NA if larger kinship sums were sampled in the simulation!
## strata is supposed to be ordered the same way than pool
## tableSimVals: if TRUE a table of the simulated kinship values is returned.
significant.individuals.althyp = function(ks, affected, pool, nr.sim,
                                          strata=NULL, tableSimVals = FALSE){
    if( length(affected)<=1 ) {
        res = data.frame(id=NA, kc=NA, freq=NA, pval=NA);
        return(res);
    }
    if(!all(affected %in% colnames(ks))){
        stop("Got affected ids for which I don't have a kinship value! This shouldn't happen!")
    }
    obs <- colSums(ks[affected, affected], na.rm=TRUE)
    names(obs) <- affected
    ## return the actual values from the simulation.
    if(is.null(strata)){
        simVals <- bg.ks.distr(ks, affected, pool, nr.sim)
    }else{
        simVals <- bg.ks.distr.strat(ks, affected, pool, nr.sim, strata)
    }
    simCount <- length(simVals)
    estimates <- do.call(rbind,lapply(obs, function(z){
        oFreq <- sum(simVals == z)/simCount * 100
        oPval <- sum(simVals >= z)/simCount
        return(c(freq=oFreq, pval=oPval))
    }))
    sumKin <- obs
    names(sumKin) <- affected
    freqKin <- estimates[, "freq"]
    freqKin[is.na(freqKin)] <- 0
    names(freqKin) <- affected
    pKin <- estimates[, "pval"]
    names(pKin) <- affected
    res <- list(sumKinship=sumKin,
                pvalueKinship=pKin,
                frequencyKinship=freqKin,
                expDensity=density(simVals),
                expHist=hist(simVals, breaks=128, plot=FALSE),
                affected=affected)
    if (tableSimVals) {
        res <- c(res, list(tableSimVals = table(simVals)))
    }
    return(res)
}


## Am just sampling using the index 1:nrow of kinship matrix.
bg.ks.distr = function(ks, affIds, pool, nr.sim){
    ## Assuming that kinship and pool are ordered the same way!
    n <- length(affIds)
    colIdx <- 1:ncol(ks)
    bg <- lapply(1:nr.sim, function(x){
        affidx <- sample(colIdx, n, replace=FALSE)
        return(colSums(ks[affidx, affidx], na.rm=TRUE))
    })
    return(unlist(bg, use.names=FALSE))
}

## same as bg.ks.distr, but for stratified sampling.
## Important is that kinship, pool and strata are in the same ordering!
bg.ks.distr.strat = function(ks, affIds, pool, nr.sim, strata=NULL){
    ksColnames <- colnames(ks)
    names(strata) <- as.character(pool)
    affStrataCounts <- table(strata[affIds])
    affStrataCounts <- affStrataCounts[affStrataCounts > 0]
    bg <- lapply(1:nr.sim, function(x){
        ## subsetting by index is faster...
        affidx <- stratsample(strata, counts=affStrataCounts)
        return(colSums(ks[affidx, affidx], na.rm=TRUE))
    })
    return(unlist(bg, use.names=FALSE))
}



#############
## plotting method...
## plotPed representing the results from the kinship clustering test.
## plotPed for FAKinSumResults does only support plotting for id.
setMethod("plotPed", "FAKinSumResults",
          function(object, id=NULL,
                   family=NULL, filename=NULL,
                   device="plot", only.phenotyped=FALSE, ...){
                       if(!is.null(family))
                           stop("Generating a pedigree for a family is not supported for FAKinshipResult. See help for more information.")
                       ## check the id, whether the id is in the result table.
                       ## res <- result(object)
                       ## if(!any(res$id==id))
                       ##     stop("The submitted id has to be present in the result table (i.e. in column affected_id of result(object))!")
                       callNextMethod(object=object, id=id, family=family, filename=filename,
                                      device=device, proband.id=id,
                                      only.phenotyped=only.phenotyped, ...)
                       ## alternatively use highlight.ids
                   })


setMethod("[", "FAKinSumResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAKinSumResults object is not supported!")
})

## This will be a crazy funky method to plot the simulation results.
setMethod("plotRes", "FAKinSumResults",
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
                  stop("plotRes for FAKinSumResults does only support specifying a single id!")
              ## check if the id is an affected for which we tested...
              if(!any(names(object@sim$pvalueKinship) == id))
                  stop(paste0("No simulation result is available for individual ", id,
                              ". id should be one of result(object)$affected_id."))
              ## OK, now we can really start doing things...
              kinSum <- object@sim$sumKinship[id]
              kinPval <- object@sim$pvalueKinship[id]
              fam <- family(object, id=id)[1, "family"]
              par(xpd=FALSE)
              if(type == "density"){
                  toplot <- object@sim$expDensity
                  XL <- range(c(range(toplot$x), kinSum))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Kinship sum", type="h",
                       lwd=3, col="lightgrey", xlim=XL)
                  points(toplot, col="grey", type="l", lwd=2)
              }
              if(type == "hist"){
                  toplot <- object@sim$expHist
                  XL <- range(c(range(toplot$mids), kinSum))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Kinship sum",
                       col="lightgrey", border="grey", xlim=XL)
              }
              Blue <- "#377EB8"
              abline(v=kinSum, col=Blue)
              if(addLegend){
                  legend("topright", legend=c(paste0("kinship sum: ", format(kinSum, digits=2)),
                                             paste0("p-value     : ", format(kinPval, digits=3))
                                             ## paste0("affect count: ", affCount),
                                             ## paste0("ctrls count : ", ctrlsCount))
                                             ))
              }
})


