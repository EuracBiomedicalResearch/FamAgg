## functions to retrieve matched controls from pedigrees for a list of input
## individuals

##***************************************************************************
##
##   Functions to get ids from the pedigree that could be used as controls
##   for the given list of ids.
##
##   Input: ped, id and any number of optional parameters as input
##   Return:
##
##***************************************************************************

## get all ids from the pedigree matching the generation from the input id(s)
## ped is assumed to have the standard format.
## if id is "spread" across several families, ids are returned for each
## family separately.
## include.ancs and include.off allow to increase the range of generations from
## which individuals are returned (e.g. include.ancs=1 all ids form the ids'
## generation plus one generation of ancestors is returned).
## If a list of pre-calculated generations is submitted the function will run
## (obviously) faster.
##### NOTE:
## If we want to consider also the frequency of a generation being
## picked/present in id: get generationMatched for each id separately,
## concatenate the ids for all runs, make a table (-> get the frequency) and
## use that as probability for the sampling.
doGetGenerationMatched <- function(ped, id=NULL, include.anc=0,
                                   include.off=0, ...){
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("No id has been submitted!")
    id <- as.character(id)
    id <- id[id %in% ped$id]
    if(length(id) == 0)
        stop("None of the ids present in the pedigree!")
    names(id) <- ped[id, "family"]
    ids <- split(id, names(id))
    lapply(ids, function(z){
        gens <- doEstimateGenerationsFor2(ped,
                                          family=names(z)[1])[[1]]
        gens <- gens[!is.na(gens)]
        gensOfIds <- gens[z]
        gensOfIds <- unique(c(gensOfIds, gensOfIds - include.anc,
                              gensOfIds + include.off))
        return(names(gens)[gens %in% gensOfIds])
    })
}

##
doGetGenerationSexMatched <- function(ped, id=NULL, include.anc=0,
                                      include.off=0, ...){
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("No id has been submitted!")
    id <- as.character(id)
    id <- id[id %in% ped$id]
    if(length(id) == 0)
        stop("None of the ids present in the pedigree!")
    names(id) <- ped[id, "family"]
    ids <- split(id, names(id))
    lapply(ids, function(z){
        gens <- doEstimateGenerationsFor2(ped,
                                          family=names(z)[1])[[1]]
        gens <- gens[!is.na(gens)]
        gensOfIds <- gens[z]
        gensOfIds <- unique(c(gensOfIds, gensOfIds - include.anc,
                              gensOfIds + include.off))
        ctrls <- names(gens)[gens %in% gensOfIds]
        ## further subset to sex:
        subPed <- ped[ped$family == names(z)[1], , drop=FALSE]
        rownames(subPed) <- as.character(subPed$id)
        haveSex <- unique(subPed[as.character(subPed$id) %in% z, "sex"])
        haveSex <- haveSex[!is.na(haveSex)]
        if(length(haveSex) == 0)
            return(NA)
        ctrlsSex <- subPed[ctrls, "sex"]
        ctrls <- ctrls[which(ctrlsSex %in% haveSex)]
        return(ctrls)
    })
}


## just get all individuals of a family.
doGetAll <- function(ped, id=NULL, ...){
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("id has to be provided!")
    id <- as.character(id)
    id <- id[id %in% ped$id]
    names(id) <- ped[id, "family"]
    ids <- split(id, names(id))
    lapply(ids, function(z){
        return(as.character(ped[as.character(ped$family) == names(z)[1], "id"]))
    })
}

## just get all individuals of a family matching by sex
doGetSexMatched <- function(ped, id=NULL, ...){
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("id has to be provided!")
    id <- as.character(id)
    id <- id[id %in% ped$id]
    names(id) <- ped[id, "family"]
    ids <- split(id, names(id))
    lapply(ids, function(z){
        ## further subset to sex:
        subPed <- ped[ped$family == names(z)[1], , drop=FALSE]
        rownames(subPed) <- as.character(subPed$id)
        haveSex <- unique(subPed[as.character(subPed$id) %in% z, "sex"])
        haveSex <- haveSex[!is.na(haveSex)]
        if(length(haveSex) == 0)
            return(NA)
        ctrlsSex <- subPed[, "sex"]
        ctrls <- subPed[which(ctrlsSex %in% haveSex), "id"]
        return(as.character(ctrls))
    })
}

doGetExternalMatched <- function(ped, id=NULL, match.using, ...){
    ## first check match.using
    if(missing(match.using))
        stop("Function getExternalMatched requires argument 'match.using'!")
    if(!is.null(dim(match.using)))
        stop("'match.using' should be a one-dimensional vector!")
    if(is.null(names(match.using)))
        stop("'match.using' has to be a named vector with the names ",
             "corresponding to the ids of the individuals in the pedigree!")
    ## check the rest...
    if(is.null(ped))
        stop("No pedigree submitted!")
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("id has to be provided!")
    id <- as.character(id)
    id <- id[id %in% ped$id]
    names(id) <- ped[id, "family"]
    ids <- split(id, names(id))
    lapply(ids, function(z){
        infam <- as.character(ped[ped$family == names(z)[1], "id"])
        if(!all(infam %in% names(match.using)))
            stop("Could not find some of the ids in 'names(match.using)'!")
        ## Now get the value in match.using for the ids, i.e. z
        needVal <- unique(match.using[z])
        needVal <- needVal[!is.na(needVal)]
        if(length(needVal) == 0)
            return(NA)
        ctrlsVal <- match.using[infam]
        ctrls <- infam[which(ctrlsVal %in% needVal)]
        return(as.character(ctrls))
    })
}


