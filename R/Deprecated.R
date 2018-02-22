##************************************************
##
##       FAProbResult
##
##       Main class defining/containing results and data.
##       * clique
##************************************************
setClass("FAProbResults",
         contains="FAData",
         slots=c(.cliques="factor",
                 sim="list",
                 nsim="numeric"),
         prototype=list(.cliques=factor(),
                        sim=list(),
                        nsim=0)
         )

setGeneric("probabilityTest", function(object, ...)
    standardGeneric("probabilityTest"))
setGeneric("cliques", function(object, ...)
    standardGeneric("cliques"))
setGeneric("cliques<-", function(object, value)
    standardGeneric("cliques<-"))
setGeneric("cliqueAndTrait", function(object, ...)
    standardGeneric("cliqueAndTrait"))
setGeneric("traitByClique", function(object)
    standardGeneric("traitByClique"))
##********************************************************************
##
##   Data analysis methods
##
##********************************************************************
## probabilistic test.
## based partially on code from Daniel Taliun
## x: FAData
## trait: trait data, need IDs and 0, 1 NA encoding.
## cliques: cliques data: two columns, first is "clique", second is ID.
## nsim: number of simulations.
## traitName: the name of the trait.
.prob_msg <- paste0("Due to problems with the 'gap' package in MS Windows the",
                    " 'probability test' will be removed in Bioconductor ",
                    "version 3.8. We apologize for any inconveniences.")
setMethod("probabilityTest", "FAData",
          function(object, trait, cliques, nsim=50000, traitName, ...){
              if(missing(cliques))
                  stop("cliques is missing!")
              if(missing(trait)){
                  ## check internal trait...
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              OrigTrait <- trait
              OrigCliques <- cliques
              ## building the result data object
              ## object <- as(object, "FATrait")
              object <- as(object, "FAProbResults")
              suppressMessages(
                  trait(object) <- trait
              )
              cliques(object) <- cliques
              if(!missing(traitName))
                  object@traitname <- traitName
              ## run the simulation: calls the runSimulation method for the
              ## FAProbResult class.
              object <- runSimulation(object, nsim=nsim, ...)
              return(object)
          })

##************************************************
##
##       FAProbResults
##
##
##************************************************
setMethod("show", "FAProbResults", function(object){
    ##cat("FAProbResults object with:\n")
    callNextMethod()
    ## print clique info.
    cat("Clique info:\n")
    cat(paste0(" * Total number of cliques: ",
               length(levels(cliques(object, na.rm=TRUE))), ".\n"))
    cat(paste0(" * Sum of individuals in pedigree not part of any clique: ",
        sum(is.na(cliques(object))), ".\n"))
    cat(paste0("Result info:\n"))
    cat(paste0(" * Dimension of result data.frame: ",
               dim(result(object)),".\n"))
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})
## get the clique information.
setMethod("cliques", "FAProbResults", function(object, na.rm=FALSE){
    .Deprecated(msg = .prob_msg)
    cl <- object@.cliques
    if(na.rm){
        cl <- cl[!is.na(cl)]
    }
    cl <- droplevels(cl)
    return(cl)
})
setReplaceMethod("cliques", "FAProbResults", function(object, value){
    .Deprecated(msg = .prob_msg)
    if(!(is.character(value) | is.numeric(value) | is.factor(value))){
        stop("'cliques' should be a numeric or character vector or a factor!")
    }
    if(is.null(names(value)))
        stop("cliques should be a named vector with the names corresponding",
             " to the ids of the individuals in the pedigree!")
    ## Convert to character; for now
    valNames <- names(value)
    value <- as.character(value)
    names(value) <- valNames
    Pedigree <- pedigree(object)
    cliqInPed <- names(value) %in% Pedigree$id
    if(!any(cliqInPed))
        stop("None of the ids in cliques (i.e. names of the input argument",
             " cliques) can be matched to ids in the pedigree")
    value <- value[cliqInPed]
    ## create a cliques vector
    cliques <- rep(NA, nrow(Pedigree))
    names(cliques) <- Pedigree$id
    ## filling with values...
    cliques[match(names(value), names(cliques))] <- value
    cliques <- factor(cliques)
    object@.cliques <- cliques
    ## Reset results, if there are some...
    if(length(object@sim) > 0){
        message("Resetting results.")
        object@sim <- list()
    }
    ## object@traitname <- character()
    return(object)
})

## calling the trait replecement method from FAResult and in addition reset
## the simulation result.
setReplaceMethod("trait", "FAProbResults", function(object, value){
    .Deprecated(msg = .prob_msg)
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})


## get a data.frame with the clique and the data from the trait for each
## individual, rownames are the individual IDs.
setMethod("cliqueAndTrait", "FAProbResults", function(object, na.rm=FALSE){
    .Deprecated(msg = .prob_msg)
    affClique <- data.frame(clique=cliques(object), trait=trait(object))
    if(na.rm){
        affClique <- affClique[!is.na(affClique[, 1]), ]
        affClique <- affClique[!is.na(affClique[, 2]), ]
    }
    droplevels(affClique)
})
## get a matrix with the size of the clique and number of affected, rownames
## correspond to the clique ID.
setMethod("traitByClique", "FAProbResults", function(object){
    .Deprecated(msg = .prob_msg)
    affClique <- cliqueAndTrait(object, na.rm=TRUE)
    CliqueSummary <- split(affClique, f=affClique$clique)
    CliqueSummary <- lapply(CliqueSummary, function(z){
        return(c(nrow(z), sum(z$trait)))
    })
    ## summary for the clique, i.e. size of the clique and number of affected.
    CliqueSummary <- do.call(rbind, CliqueSummary)
    colnames(CliqueSummary) <- c("size" , "affected_count")
    CliqueSummary
})

#######
## the analysis method:
setMethod("runSimulation", "FAProbResults", function(object, nsim=50000){
    .Deprecated(msg = .prob_msg)
    ## only makes sense if we've got trait and cliques
    if(length(trait(object)) == 0)
        stop("No trait information available!")
    if(length(cliques(object)) == 0)
        stop("No cliques information avaliable!")
    ## performing the test
    ## based on trait and clique...
    ## affClique <- getAffectedInClique(x)
    ## no of affected per clique and sum of clique members:
    CliqueSummary <- traitByClique(object)

    ## Dummy info...
    message("Cleaning data set (got in total ", length(object$id),
            " individuals):")
    message(" * not phenotyped individuals...", appendLF=FALSE)
    notPhen <- sum(is.na(trait(object)))
    if(notPhen > 0){
        message(" ", notPhen, " removed.")
    }else{
        message(" none present.")
    }
    message("Done.")

    ## what is the frequency???
    l_freq <- aggregate(rep(1, nrow(CliqueSummary)),
                        by=list(CliqueSummary[, "size"],
                                CliqueSummary[, "affected_count"]),
                        FUN=sum)
    ## ASK X_: what's that frequency???
    ## is it the size of a group, the number of affected and the number of
    ## times that this group size/affected count is occurring?
    colnames(l_freq) <- c("n", "affected", "freq")
    l_freq <- as.matrix(l_freq)

    ## Test if there are any cliques which size exceeds the max size
    ## supported by the gap package
    TooLarge <- l_freq[, "n"] > GAP_MAX_CLIQUE_SIZE
    if(sum(TooLarge) > 0)
        stop(sum(TooLarge), " cliques are larger than supported by the Monte ",
             "Carlo simulation in the gap package! Consider defining smaller ",
             "cliques!")
    if (.Platform$OS.type != "unix")
        warning("On Windows the 'gap' package on which this test is based on ",
                "seems to be buggy. You might encounter errors due to that.")
    fc.sim <- pfc.sim(l_freq, n.sim=nsim)
    object@nsim <- nsim
    object@sim <- fc.sim
    object
})

## results; get the results table.
setMethod("result", "FAProbResults", function(object, method="BH"){
    .Deprecated(msg = .prob_msg)
    if(length(object@sim)==0){
        stop("No simulation performed yet! Please use the probabilityTest ",
             "function or the runSimulation method to start the simulation.")
    }
    method <- match.arg(method, p.adjust.methods)
    ## generate the result table...
    ## get a data.frame with clique
    cltr <- cliqueAndTrait(object, na.rm=TRUE)
    ## get a summary matrix: size and affected per clique.
    ## note: that's the clique summary data for all individuals per clique
    ## for which we do have data in the trait!!! thus it's different from
    ## the table(cliques(Test))
    CliqSum <- traitByClique(object)
    Trait <- trait(object)
    ## subset to those which are part of any of the cliques...
    TraitNClique <- Trait[!is.na(cliques(object))]
    Trait <- Trait[!is.na(Trait)]
    TraitNClique <- TraitNClique[!is.na(TraitNClique)]
    ## Question: what is the Pedigree Size???
    TraitName <- object@traitname
    if(length(TraitName)==0)
        TraitName <- NA
    ## Get the ids and cliques
    Cl <- cliques(object)
    Cl <- Cl[!is.na(Cl)]
    ClIds <- names(Cl)
    names(ClIds) <- as.character(Cl)
    ## Get one representative ID for each clique
    ClIds <- ClIds[rownames(CliqSum)]
    MyRes <- data.frame(trait_name=TraitName,
                        total_phenotyped=length(Trait),
                        total_affected=sum(Trait!=0),
                        phenotyped=length(TraitNClique),
                        affected=sum(TraitNClique!=0),
                        group_id=rownames(CliqSum),
                        family=pedigree(object)[ClIds, "family"],
                        group_phenotyped=CliqSum[, "size"],
                        group_affected=CliqSum[, "affected_count"],
                        pvalue=object@sim$tailpu,
                        padj=p.adjust(object@sim$tailpu, method=method),
                        check.names=FALSE, stringsAsFactors=FALSE)
    rownames(MyRes) <- as.character(MyRes$group_id)
    return(MyRes)
})


#############
## plotting method...
## plotPed representing the results from the probabilistic test.
## plotPed for FAProbResults does only support plotting for id with the id
## being the group id.
## TODO: fix that: should I only plot the pedigree of the clique or of the
## full pedigree highlighing the clique?
setMethod("plotPed", "FAProbResults", function(object, id = NULL,
                                               family = NULL, filename = NULL,
                                               device = "plot", ...) {
    .Deprecated(msg = .prob_msg)
    if (!is.null(family))
        stop("Generating a pedigree for a family is not supported for ",
             "FAProbResults. See help for more information.")
    cliqs <- cliques(object)
    res <- result(object)
    if(!any(res$group_id == id))
        stop("The id should be one of the group ids (clique names; column ",
             "group_id in result(object))!")
    ## get all individuals of the group...
    clique <- names(cliqs)[which(as.character(cliqs) == as.character(id))]
    highlight.ids <- list(`*`=clique)  ## alternatively, use id=clique!
    callNextMethod(object=object, id=id, family=family,
                   filename=filename, proband.id=clique,
                   device=device,...)
})

## this buildPed simply ensures that all phenotyped in the group will be
## included in the pedigree.
## TODO: how to build the pedigree?
setMethod("buildPed", "FAProbResults", function(object, id = NULL,
                                                max.generations.up = 3,
                                                max.generations.down = 16,
                                                prune = FALSE) {
    .Deprecated(msg = .prob_msg)
    if(is.null(id))
        stop("The id of the group has to be speficied!")
    cliqs <- cliques(object)
    res <- result(object)
    if(!any(res$group_id == id))
        stop("The id should be one of the group ids (clique names; column ",
             "group_id in result(object))!")
    ## get all individuals of the group...
    clique <- names(cliqs)[which(as.character(cliqs) == as.character(id))]
    ped <- callNextMethod(object=object, id=clique, prune=prune,
                          max.generations.up=max.generations.up,
                          max.generations.down=max.generations.down)
    return(ped)
})

## this is to get all those that are related with any individuals of the group!
setMethod("shareKinship", "FAProbResults", function(object, id = NULL) {
    if(is.null(id))
        stop("The id of the group has to be speficied!")
    cliqs <- cliques(object)
    res <- result(object)
    if(!any(res$group_id == id))
        stop("The id should be one of the group ids (clique names; column ",
             "group_id in result(object))!")
    ## get all individuals of the group...
    clique <- names(cliqs)[which(as.character(cliqs) == as.character(id))]
    return(doShareKinship(kin=kinship(object), id=clique))
})

setMethod("[", "FAProbResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAProbResults object is not supported!")
})


