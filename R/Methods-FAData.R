## Methods for the classes
setMethod("pedigree", "missing", function(...){
    return(kinship2::pedigree(...))
})


##************************************************
##
##       Methods for FAData objects.
##
##************************************************
## validate the FAData object. Checks if the pedigree and age fit/match.
validateFAData <- function(object){
    Complaints <- character()
    if(!is.null(object@pedigree)){
        ## check the number of columns and column names.
        ## check the columns...
        value <- object@pedigree
        CN <- .PEDCN
        ##
        if(is.null(noColNames(value)) & ncol(value) == length(CN)){
            colnames(Value) <- CN
        }
        ## check if all expected colnames are there:
        if(all(CN %in% colnames(value))){
            ## all fine, just drop additional columns
            value <- value[, CN]
        }else{
            ## Missing one or more colnames!!!
            missingCN <- CN[!(CN %in% colnames(value))]
            stop(paste0("Required columns ",
                        paste(missingCN, collapse=", "),
                        " not present in pedigree!"))
        }
        object@pedigree$sex <- sanitizeSex(object@pedigree$sex)
    }
    if(length(object@age)>0){
        ## has to be a named numeric vector!
        if(is.null(names(object@age))){
            Complaints <- c(Complaints, "\nAge has to be a named numeric vector!")
        }
    }
    ## check the trait:
    if(length(object@.trait) > 0){
        ## trait has to have the same length than the pedigree has rows.
        if(length(object@.trait) != nrow(object@pedigree))
            Complaints <- c(Complaints,
                            "\ntrait has to have the same length than the pedigree has columns!")
    }
    if(length(Complaints)>0){
        return(Complaints)
    }else{
        return(TRUE)
    }
}
setValidity("FAData", validateFAData)
setMethod("initialize", "FAData", function(.Object,...){
    OK <- validateFAData(.Object)
    if(is.character(OK)){
        stop(OK)
    }
    callNextMethod(.Object, ...)
})
setMethod("show", "FAData", function(object){
    cat(paste0(class(object), " object with:\n"))
    cat(paste0("* Pedigree of length ", nrow(object@pedigree), ".\n"))
    if(!is.null(object@pedigree)){
        if(length(object@pedigree) > 0 & nrow(object@pedigree) > 0){
            UnIds <- unique(object@pedigree[, "id"])
            cat(paste0("* Number of unique individuals: ", length(UnIds), ".\n"))
            if(length(age(object)) > 0){
                ## check for how many individuals we have the age...
                cat(paste0("* Number of individuals with known age: ",
                           sum(names(age(object)) %in% UnIds), ".\n"))
            }
        }
        if(length(object@.trait) > 0){
            if(length(object@traitname)==0){
                cat("Result on an unnamed trait\n")
            }else{
                cat("Result on trait named ", object@traitname, "\n")
            }
            ## print trait info.
            cat("Trait info:\n")
            cat(paste0(" * Number of non-NA values: ", length(trait(object, na.rm=TRUE)),
                       ".\n"))
            cat(paste0(" * Number of non-zero values: ", sum(trait(object, na.rm=TRUE)!=0),
                       ".\n"))
        }
    }
})


########
## getter/setter for pedigree.
setMethod("pedigree", "FAData", function(object, return.type="data.frame"){
    if(length(object@pedigree)==0){
        return(NULL)
    }
    return.type <- match.arg(return.type, c("data.frame", "pedigree"))
    ped <- object@pedigree
    CN <- .PEDCN
    if(ncol(ped)!=length(CN))
        stop("pedigree has more columns than expected!")
    if(noColNames(ped)){
        colnames(ped) <- CN
    }
    if(!is.factor(ped$sex)){
        ped$sex <- sanitizeSex(ped$sex)
    }
    ## should we add the trait information?
    if(length(object@.trait) > 0){
        tr <- trait(object)
        ped <- cbind(ped, affected=tr)
    }
    if(return.type=="pedigree"){
        if(any(colnames(ped) == "affected")){
            ped <- kinship2::pedigree(famid=ped[, "family"],
                                      id=ped[, "id"],
                                      dadid=ped[, "father"],
                                      momid=ped[, "mother"],
                                      sex=ped[, "sex"],
                                      affected=ped[, "affected"])
        }else{
            ped <- kinship2::pedigree(famid=ped[, "family"],
                                      id=ped[, "id"],
                                      dadid=ped[, "father"],
                                      momid=ped[, "mother"],
                                      sex=ped[, "sex"])
        }
    }
    return(ped)
})


## pedigree<- internally checks validity of the submitted value.
setReplaceMethod("pedigree", "FAData", function(object, value){
    if(is.data.frame(value)){
        CN <- .PEDCN
        if(!noColNames(value)){
            ## check the colnames
            if(sum(CN %in% colnames(value)) != length(CN)){
                stop(paste0("pedigree is expected to have column names ",
                            paste(CN, collapse=", "), "."))
            }
            value <- value[, CN]
        }else{
            ## no col names, thus assuming they fit what we expect.
            if(ncol(value)!=length(CN))
                stop(paste0("pedigree is expected to have ", length(CN),
                            " columns but I got ", ncol(value), "!"))
            colnames(value) <- CN

        }
        ## check and sanitize the sex...
        value[, "sex"] <- sanitizeSex(value$sex)
        rownames(value) <- as.character(value$id)
    }else if(is(value, "pedigree") | is(value, "pedigreeList")){
        value <- ped2df(value)
    }else{
        stop("pedigree has to be a data.frame, a pedigree or a pedigreeList object!")
    }

    ## next remove all white spaces!
    for(theCol in c("id", "family", "father", "mother")){
        tmp <- value[, theCol]
        M <- mode(tmp)
        tmp <- gsub(tmp, pattern="\\s", replacement="")
        mode(tmp) <- M
        value[, theCol] <- tmp
        rm(tmp)
    }
    object@pedigree <- value
    validateFAData(object)

    ## save the kinship matrix
    message("Generating the kinship matrix...", appendLF=FALSE)
    object@.kinship <- kinship2::kinship(pedigree(object, return.type="pedigree"))
    message("OK\n")

    return(object)
})


########
## get pedigree size
setMethod("pedigreeSize", "FAData", function(object){
    return(nrow(pedigree(object)))
})


## accessor for columns in pedigree...
setMethod("$", "FAData", function(x, name){
    tn <- x@traitname
    if(length(tn) == 0)
        tn <- ""
    if(name == "trait" | name == tn){
        return(trait(x))
    }
    vals <- eval(substitute(pedigree(x)$NAME_ARG, list(NAME_ARG=name)))
    if(length(vals) > 0){
        ## if we're accessing the id let's transform that to character
        if(name == "id")
            vals <- as.character(vals)
        ids <- pedigree(x)$id
        if(length(vals)==length(ids))
            names(vals) <- ids
    }
    return(vals)
})

## get the ids of all affected individuals in the trait
setMethod("affectedIndividuals", "FAData", function(object){
    if(length(object@.trait) == 0){
        warning("No trait information available")
        return(NULL)
    }
    tr <- trait(object, na.rm=TRUE)
    return(names(tr)[tr!=0])
})
## get the ids of all phenotyped individuals in the trait
setMethod("phenotypedIndividuals", "FAData", function(object){
    if(length(object@.trait) == 0){
        warning("No trait information available")
        return(NULL)
    }
    tr <- trait(object, na.rm=TRUE)
    return(names(tr))
})


########
## getter/setter for age.
setMethod("age", "FAData", function(object){
    ## we do have a pedigree, so we make sure we return ages in the same
    ## order than ids in the pedigree!
    if(length(object@pedigree) > 0){
        ped <- pedigree(object, return.type="data.frame")
        ages <- rep(NA, nrow(ped))
        names(ages) <- ped[, "id"]
        if(length(object@age) > 0){
            ageids <- intersect(names(ages), names(object@age))
            ages[ageids] <- object@age[ageids]
            return(ages)
        }else{
            return(ages)
        }
    }else{
        ## we don't have any pedigree, so just return ages:
        return(object@age)
    }
})
setReplaceMethod("age", "FAData", function(object, value){
    if(!is.numeric(value) | is.null(names(value)))
        stop("value has to be a named numeric vector.")
    object@age <- value
    ## check if the object is (still) valid
    validateFAData(object)
    return(object)
})
## get the trait information.
setMethod("trait", "FAData", function(object, na.rm=FALSE){
    tr <- object@.trait
    if(na.rm){
        return(tr[!is.na(tr)])
    }else{
        return(tr)
    }
})
## trait will be matched to the ids in the pedigree of object.
setReplaceMethod("trait", "FAData", function(object, value){
    if(is.logical(value))
        value <- as.numeric(value)
    if(is.numeric(value)){
        ## check that trait is 0, 1, NA
        if(!all(unique(value) %in% c(0, 1, NA)))
            stop(paste0("trait should be a named logical vector or numeric",
                        " vector with values 0, 1 and NA!"))
    }else{
        stop("trait has to be either a named logical or numerical vector!")
    }
    if(is.null(names(value)))
        stop(paste0("trait has to be a named vector with the names",
                    " corresponding to the IDs used in the pedigree"))
    Pedigree <- pedigree(object)
    traitInPed <- names(value) %in% Pedigree$id
    ## subsetting trait to those...
    if(!any(traitInPed))
        stop(paste0("None of the ids in trait (i.e. names of the input",
                    " argument trait) can be matched to ids in the pedigree!"))
    value <- value[traitInPed]
    message(paste0(sum(traitInPed), " of in total ", length(traitInPed),
                   " trait values can be matched to IDs in the pedigree."))
    trait <- as.numeric(rep(NA, nrow(Pedigree)))
    names(trait) <- Pedigree$id
    ## filling with values...
    trait[match(names(value), names(trait))] <- value
    object@.trait <- trait
    return(object)
})



########
## getter for kinship. if not present in the .kinship slot it will be
## generated and saved.
setMethod("kinship", "FAData", function(id, ...){
    if(nrow(id@.kinship)==0){
        kin <- kinship2::kinship(pedigree(id, return.type="pedigree"))
    }else{
        kin <- id@.kinship
    }
    return(kin[as.character(id$id), as.character(id$id)])
})
#######
## get the family given an id or a family id.
setMethod("family",
          "FAData",
          function(object, id=NULL, family=NULL, return.type="data.frame"){
              if(missing(id))
                  id <- NULL
              if(missing(family))
                  family <- NULL
              return.type <- match.arg(return.type, c("data.frame", "pedigree"))
              if(is.null(id) & is.null(family))
                  stop("Either id or family has to be specified!")
              ## check first if we have id specified...
              ped <- pedigree(object)
              if(!is.null(id)){
                  ## get the ids of all families for these id(s)
                  family <- as.character(
                      ped[as.character(ped$id) %in% as.character(id), "family"]
                  )
              }
              if(!is.null(family) & length(family) > 0){
                  ped <- ped[as.character(ped$family) %in% family, ]
              }else{
                  return(NULL)
              }
              if(return.type=="pedigree")
                  if(any(colnames(ped)=="affected")){
                      ped <- pedigree(id=ped$id, dadid=ped$father,
                                      momid=ped$mother, sex=ped$sex,
                                      famid=ped$family, affected=ped$affected)
                  }else{
                      ped <- pedigree(id=ped$id, dadid=ped$father,
                                      momid=ped$mother, sex=ped$sex,
                                      famid=ped$family)
                  }
              return(ped)
          })

setMethod("buildPed", "FAData",
          function(object, id=NULL, max.generations.up=3, max.generations.down=16,
                   prune=FALSE, ...){
              ped <- pedigree(object)
              if(is.null(id))
                  return(ped)
              id <- as.character(id)
              notThere <- !(id %in% ped$id)
              if(any(notThere)){
                  id <- id[!notThere]
                  warning("Removed ", sum(notThere), " ids since they are not in the pedigree!")
                  if(length(id) == 0)
                      stop("No id available!")
              }
              if((length(id) == 1) & prune){
                  prune <- FALSE
                  warning("Setting prune=FALSE since there is only a single id!")
              }
              if(prune){
                  ## restrict to the smalles pedigree including the ids
                  return(subPedigree(ped, id=id))
              }
              ## else grow the full pedigree
              ## first get all ancestors:
              ancs <- getAncestors(object, id=id, max.generations=max.generations.up)
              ## add eventual missing mates.
              mismate <- doGetMissingMate(ped, id=ancs)
              ## next get all children:
              allids <- unique(c(ancs, id, mismate))
              chlds <- getChildren(object, id=allids,
                                   max.generations=max.generations.down)
              ## get all missing mates:
              allids <- unique(c(allids, chlds))
              allmates <- getMissingMate(object, id=allids)
              allids <- unique(c(id, ancs, chlds, allmates))
              ped <- ped[as.character(allids), , drop=FALSE]
              ## fixfounders:
              ## set all mothers and fathers that are not in id to 0
              ped[!(ped$mother %in% ped$id), "mother"] <- 0
              ped[!(ped$father %in% ped$id), "father"] <- 0
              ped <- doRemoveChildlessFounders(ped)
              return(ped)
          })


## [ subsetting.
.bracketSubset <- function(x, i, j, ..., drop){
    if(!missing(j))
        stop("Subsetting by columns ('j') is not supported")
    haveRows <- nrow(pedigree(x))
    if(length(haveRows) == 0){
        warning("Can not subset an empty object.")
        return(x)
    }
    if(missing(i)){
        i <- 1:haveRows
    }
    ## check i:
    i <- .checkRowindex(i, pedigree(x))
    ## don't use the accessor here!
    pedSub <- x@pedigree[i, , drop=FALSE]
    ageSub <- numeric()
    if(length(x@age) > 0){
        ## subsetting age.
        ageSub <- x@age[names(x@age) %in% rownames(pedSub)]
    }
    ## have to use character colnames here.
    kinSub <- kinship(x)[as.character(pedSub$id), as.character(pedSub$id), drop=FALSE]
    newX <- new("FAData", pedigree=pedSub, age=ageSub, .kinship=kinSub)
    ## subsetting trait...
    if(length(x@.trait) > 0){
        subTrait <- trait(x)[i]
        trait(newX) <- subTrait
        newX@traitname <- x@traitname
    }
    return(newX)
}
setMethod("[", "FAData", .bracketSubset)

## i can be a numeric, character or logical vector which is supposed to be used
## for subsetting data.
## data is a data.frame or matrix with rownames.
## the function checks if i is valid to subset data and returns a numeric vector
## that can be used for subsetting.
.checkRowindex <- function(i, data){
    haveRows <- nrow(data)
    if(!is.logical(i) & !is.numeric(i) & !is.character(i))
        stop("'i' has to be either a logical, numeric or character vector!")
    ## i can be boolean -> has to be the same length than haveRows
    if(is.logical(i)){
        if(length(i) != haveRows)
            stop(paste0("If 'i' is a logical vector its length has to match",
                        " the number individuals in the pedigree!"))
        ## transform to numeric
        i <- which(i)
    }
    ##          character -> have to match rownames of pedigree.
    if(is.character(i)){
        iLen <- length(i)
        i <- match(i, rownames(data))
        if(all(is.na(i))){
            stop(paste0("None of the elements in 'i' matches an id of an",
                        " individual (i.e. rownames of pedigree)!"))
        }
        if(any(is.na(i))){
            warning(paste0(sum(is.na(i)), " elements in 'i' can not be",
                           " matched to ids of individuals and were thus discarded."))
            i <- i[!is.na(i)]
        }
    }
    ##          numeric -> indices have to be within haveRows
    iLen <- length(i)
    i <- i[i >=1 & i <=haveRows]
    if(length(i)==0)
        stop("'i' has to be a numeric between 1 and ", haveRows, "!")
    if(iLen != length(i))
        warning(paste0("Some of the values in 'i' are outside of the",
                       " allowed range [1,",haveRows, "] and were thus discarded"))
    ## have it.
    return(i)
}




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
              trait(object) <- trait
              cliques(object) <- cliques
              if(!missing(traitName))
                  object@traitname <- traitName
              ## run the simulation: calls the runSimulation method for the FAProbResult class.
              object <- runSimulation(object, nsim=nsim, ...)
              return(object)
          })
##
## kinship test.
## based partially on code from Daniel Taliun.
setMethod("kinshipGroupTest", "FAData",
          function(object, trait, nsim=50000, traitName, strata=NULL, ...){
              if(missing(trait)){
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              OrigTrait <- trait
              ## building the result data object
              ## object <- as(object, "FATrait")
              object <- as(object, "FAKinGroupResults")
              trait(object) <- trait
              if(!missing(traitName))
                  object@traitname <- traitName
              ## run the simulation: calls the runSimulation method for the FAKinshipResult class.
              object <- runSimulation(object, nsim=nsim, strata=strata, ...)
              return(object)
          })

## kinship clustering test.
setMethod("kinshipSumTest", "FAData",
          function(object, trait, nsim=50000, traitName, strata=NULL, ...){
              if(missing(trait)){
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              ## object <- as(object, "FAResult")
              object <- as(object, "FAKinSumResults")
              trait(object) <- trait
              if(!missing(traitName))
                  object@traitname <- traitName
              ## run the simulation.
              object <- runSimulation(object, nsim=nsim, strata=strata, ...)
              return(object)
          })

## genealogical index
setMethod("genealogicalIndexTest", "FAData",
          function(object, trait, nsim=50000, traitName,
                   perFamilyTest=FALSE, controlSetMethod="getAll",
                   prune=TRUE, strata=NULL, ...){
              if(missing(trait)){
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              ## building the result data object
              ## object <- as(object, "FAResult")
              object <- as(object, "FAGenIndexResults")
              trait(object) <- trait
              if(!missing(traitName))
                  object@traitname <- traitName
              ## run the simulation: calls the runSimulation method for the FAProbResult class.
              object <- runSimulation(object, nsim=nsim, perFamilyTest=perFamilyTest,
                                      controlSetMethod=controlSetMethod,
                                      prune=prune, strata=strata, ...)
              return(object)
          })


##****************************************************************************
##
##  familial incidence rate as described in Kerber
##
##
setMethod("familialIncidenceRate", "FAData",
          function(object, trait=NULL, timeAtRisk=NULL, prune=TRUE){
              ## .FR is defined in Methods-FAIncidenceRatio.R
              FR <- .FR(ped=pedigree(object), kin=kinship(object), trait=trait,
                        timeAtRisk=timeAtRisk, prune=prune, perFamilyTest=FALSE)
              return(FR[[1]])
          })

##****************************************************************************
##
##  familialIncidenceRateTest method.
##
setMethod("familialIncidenceRateTest", "FAData",
          function(object, trait=NULL, nsim=50000, traitName=NULL, timeAtRisk=NULL,
                   prune=TRUE, strata=NULL, ...){
              if(is.null(trait)){
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              object <- as(object, "FAIncidenceRateResults")
              trait(object) <- trait
              if(!is.null(traitName))
                  object@traitname <- traitName
              ## run the simulation
              object <- runSimulation(object, nsim=nsim, timeAtRisk=timeAtRisk,
                                      strata=strata, prune=prune, ...)
              return(object)
          })


##****************************************************************************
##
##  fsir method.
##  familial standardized incidence rate as described in Kerber
##
setMethod("fsir", "FAData", function(object, trait=NULL, lambda=NULL, timeInStrata=NULL,
                                     prune=TRUE){
    ## First we need to do a lot of checking and testing.
    if(is.null(trait)){
        ## check internal trait...
        if(length(object@.trait) == 0)
            stop("trait is missing!")
        trait <- trait(object)
    }
    ## Check lambda
    if(is.null(lambda))
        stop("lambda missing!")
    if(is.null(names(lambda)))
        stop(paste0("lambda has to be a named vector with the names",
                    "corresponding to the strata names!"))
    ## Check timeInStrata
    if(is.null(timeInStrata))
        stop("timeInStrata missing!")
    if(class(timeInStrata)!="matrix")
        stop("timeInStrata has to be a matrix!")
    if(nrow(timeInStrata) != length(object$id))
        stop(paste0("timeInStrata has to have the same number of rows as",
                    " there are individuals in the pedigree!"))
    if(length(lambda)!=ncol(timeInStrata))
        stop("length of lambda has to match the number of columns of timeInStrata!")
    if(!all(colnames(timeInStrata) %in% names(lambda)))
        stop("Names of lambda does not match the colnames of timeInStrata!")
    lambda <- lambda[colnames(timeInStrata)]
    ## Done.
    trait(object) <- trait
    trait <- trait(object)   # that way we ensure that we have the same ordering.
    kin <- kinship(object)
    ## Just to be on the save side... ensure that the ordering of id/Trait matches the kin
    kin <- kin[names(trait), names(trait)]
    diag(kin) <- 0
    ## Now start subsetting the data:
    ## * NA in trait
    nas <- is.na(trait)
    if(any(nas)){
        warning(paste0("Excluding ", sum(nas),
                       " individuals because of a missing value in the trait."))
        trait <- trait[!nas]
        kin <- kin[!nas, !nas]
        timeInStrata <- timeInStrata[!nas, , drop=FALSE]
    }
    ## * NA in timeInStrata.
    nas <- apply(timeInStrata, MARGIN=1, function(z){
        any(is.na(z))
    })
    if(any(nas)){
        warning(paste0("Excluding ", sum(nas),
                       " individuals because of a missing value in timeInStrata."))
        trait <- trait[!nas]
        kin <- kin[!nas, !nas]
        timeInStrata <- timeInStrata[!nas, , drop=FALSE]
    }
    if(prune){
        ## * Not related, i.e. individuals with a kinship sum of 0
        nas <- colSums(kin) == 0
        if(any(nas)){
            warning(paste0("Excluding ", sum(nas),
                           " individuals because they do not share kinship",
                           " with any individual in the pedigree."))
            trait <- trait[!nas]
            kin <- kin[!nas, !nas]
            timeInStrata <- timeInStrata[!nas, , drop=FALSE]
        }
    }
    ## Well done. Now let's do the test:
    fsirs <- doFsir(affected=trait, kin=kin, lambda=lambda, timeInStrata=timeInStrata)
    ## Prepare the results vector.
    allIds <- object$id
    allFsirs <- rep(NA, length(allIds))
    names(allFsirs) <- allIds
    allFsirs[names(trait)] <- as.numeric(fsirs)
    return(allFsirs)
})

##****************************************************************************
##
##  fsirTest method.
##
setMethod("fsirTest", "FAData",
          function(object, trait=NULL, nsim=50000, traitName=NULL, lambda=NULL,
                   timeInStrata=NULL, prune=TRUE, strata=NULL, ...){
              if(is.null(trait)){
                  if(length(object@.trait) == 0)
                      stop("trait is missing!")
                  trait <- trait(object)
              }
              object <- as(object, "FAStdIncidenceRateResults")
              trait(object) <- trait
              if(!is.null(traitName))
                  object@traitname <- traitName
              ## run the simulation
              object <- runSimulation(object, nsim=nsim, lambda=lambda,
                                      timeInStrata=timeInStrata,
                                      strata=strata, prune=prune, ...)
              return(object)
          })



.checkLabels <- function(label, fam){
    if(is.null(label))
        return(NULL)
    if(length(label) == nrow(fam)){
        names(label) <- fam$id
    }else{
        ## do we have names?
        if(is.null(names(label))){
            warning("Label vector does not have names and was thus discarded!")
            return(NULL)
        }
    }
    ## check and match by names.
    subset <- names(label) %in% fam$id
    if(sum(subset) == 0){
        warning(paste0("None of the names of label match to ids in the",
                       " pedigree! Argument label was thus discarded!"))
        return(NULL)
    }
    label <- label[subset]
    NAs <- is.na(label)
    if(is.numeric(label)){
        label <- format(label, digits=2)
        label[NAs] <- ""
    }
    realLabel <- rep(NA, nrow(fam))
    names(realLabel) <- fam$id
    realLabel[names(label)] <- label
    return(realLabel)
}

##********************************************************************
##
##   Plotting
##
##********************************************************************
## plotting method...
## if family not NULL: plot the pedigree for the FULL family.
## if id is not NULL: build the pedigree around that id and plot that.
## highlight.ids: character vector of ids or named list with character vector(s)
setMethod("plotPed", "FAData",
          function(object, id=NULL, family=NULL, filename=NULL,
                   device="plot", symbol.related=NA, proband.id=NULL, highlight.ids=NULL,
                   only.phenotyped=FALSE,
                   label1=age(object), label2=NULL, label3=NULL, ...){
              ## if id was defined build the pedigree for that individual,
              ## otherwise plot the full family.
              Args <- list(...)
              if(!is.null(id)){
                  ## only extract arguments we might need for buildPed:
                  max.generations.up <- 3
                  max.generations.down <- 16
                  prune <- FALSE
                  if(any(names(Args) == "max.generations.up"))
                      max.generations.up <- Args$max.generations.up
                  if(any(names(Args) == "max.generations.down"))
                      max.generations.down <- Args$max.generations.down
                  if(any(names(Args) == "prune"))
                      prune <- Args$prune
                  fam <- buildPed(object, id=id,
                                  max.generations.up=max.generations.up,
                                  max.generations.down=max.generations.down,
                                  prune=prune)
              }else{
                  fam <- family(object, id=id, family=family)
              }
              if(nrow(fam)==0)
                  stop("No data left for plotting after sub-setting.")
              ## check if we've got affected
              if(any(colnames(fam)=="affected")){
                  affected <- fam[, "affected"]
                  haveAffected <- TRUE
              }else{
                  affected <- rep(NA, nrow(fam))
                  names(affected) <- as.character(fam$id)
                  haveAffected <- FALSE
              }
              if(is.null(symbol.related))
                  symbol.related <- NA
              ## checking labels.
              label1 <- .checkLabels(label1, fam)
              label2 <- .checkLabels(label2, fam)
              label3 <- .checkLabels(label3, fam)
              ## ## check ages...
              ## ages <- age(object)[as.character(fam$id)]
              ## ages[!is.na(ages)] <- paste0(format(ages[!is.na(ages)], digits=2), "y")
              ## ages[is.na(ages)] <- ""
              ## do not have affected... obviously...
              is.proband <- rep(FALSE, nrow(fam))
              names(is.proband) <- as.character(fam$id)
              if(!missing(proband.id)){
                  proband.id <- as.character(proband.id)
                  is.proband[names(is.proband) %in% proband.id] <- TRUE
                  if(sum(is.proband)!=length(proband.id))
                      warning("Not all probands specified in proband.id are in the pedigree!")
              }
              text.inside.symbol <- rep("", nrow(fam))
              names(text.inside.symbol) <- as.character(fam$id)
              if(!is.null(id)){
                  ## is.proband[as.character(id)] <- TRUE
                  ## get individuals that share kinship
                  ##related <- doShareKinship(ped=fam, id=id)
                  related <- shareKinship(object, id=id)
                  related <- related[related %in% fam$id]
                  ## could also use shareKinship(object, id), but that's not that efficient!
                  related <- related[related %in% names(text.inside.symbol)]
                  text.inside.symbol[related] <- symbol.related
              }
              text2.below.symbol=NULL
              text3.below.symbol=NULL
              text4.below.symbol=NULL
              if(!is.null(highlight.ids)){
                  if(is.character(highlight.ids)){
                      highlight.ids <- list(`*`=highlight.ids)
                  }
                  if(is.list(highlight.ids)){
                      ## support at max 3 highlight lines.
                      for(i in 1:min(c(3, length(highlight.ids)))){
                          if(is.null(names(highlight.ids)[i])){
                              symb <- "*"
                          }else{
                              symb <- names(highlight.ids)[i]
                          }
                          texts <- rep("", nrow(fam))
                          names(texts) <- as.character(fam$id)
                          texts[names(texts) %in% highlight.ids[[i]]] <- symb
                          if(i==1)
                              text2.below.symbol <- texts
                          if(i==2)
                              text3.below.symbol <- texts
                          if(i==3)
                              text4.below.symbol <- texts
                      }
                  }else{
                      warning(paste0("Discarding argument highlight.ids. It",
                                     " should be a character vector of ids or a",
                                     " (named) list of character vectors with ids."))
                  }
              }
              ## Note that label1 to label3 has precedence to any other argument!
              if(!is.null(label2))
                  text2.below.symbol <- label2
              if(!is.null(label3))
                  text3.below.symbol <- label3
              ## If we want to plot only phenotyped individuals...
              if(only.phenotyped & haveAffected){
                  ## use the buildPed...
                  keepIds <- fam[!is.na(fam[, "affected"]), "id"]
                  subped <- subPedigree(fam, id=as.character(keepIds), all=FALSE)
                  plotMe <- as.character(fam$id) %in% as.character(subped$id)
              }else{
                  plotMe <- rep(TRUE, nrow(fam))
              }
              ## OK, now plot!!!
              res <- doPlotPed(family=fam$family[plotMe], individual=fam$id[plotMe],
                               father=fam$father[plotMe], mother=fam$mother[plotMe],
                               gender=fam$sex[plotMe], is.proband=is.proband[plotMe],
                               text1.below.symbol=label1[plotMe],
                               text2.below.symbol=text2.below.symbol[plotMe],
                               text3.below.symbol=text3.below.symbol[plotMe],
                               text4.below.symbol=text4.below.symbol[plotMe],
                               text.inside.symbol=text.inside.symbol[plotMe],
                               affected=affected[plotMe],
                               filename=filename, device=device, ...)
              invisible(res)
          })





##********************************************************************
##
##   Pedigree utilities
##
##********************************************************************
setMethod("countGenerations", "FAData",
          function(object, id=NULL, direction="down", ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doCountGenerations(pedigree(object), id=id, direction=direction))
          })
setMethod("estimateGenerations", "FAData",
          function(object, family=NULL, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doEstimateGenerationsFor2(pedigree(object), family=family))
          })
setMethod("findFounders", "FAData",
          function(object, family=NULL, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doFindFounders(pedigree(object), family=family))
          })
setMethod("generationsFrom", "FAData",
          function(object, id=NULL, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doGetGenerationFrom2(pedigree(object), id=id, ...))
          })
setMethod("getAncestors", "FAData",
          function(object, id=NULL, max.generations=3, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doGetAncestors(pedigree(object), id=id, maxlevel=max.generations, ...))
          })
setMethod("getChildren", "FAData",
          function(object, id=NULL, max.generations=16, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doGetChildren(pedigree(object), id=id, maxlevel=max.generations, ...))
          })
setMethod("getMissingMate", "FAData",
          function(object, id=NULL, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doGetMissingMate(pedigree(object), id=id, ...))
          })
setMethod("getSiblings", "FAData",
          function(object, id=NULL, ...){
              if(nrow(pedigree(object)) == 0)
                  stop("No pedigree available!")
              return(doGetSiblings(pedigree(object), id=id, ...))
          })

setMethod("shareKinship", "FAData",
          function(object, id=NULL){
              if(is.null(id))
                  stop("id has to be specified!")
              return(doShareKinship(kin=kinship(object), id=id))
          })

setMethod("getCommonAncestor", "FAData",
          function(object, id, method="min.dist"){
              return(doGetCommonAncestor(pedigree(object), id=id,
                                         method=method))
          })



##***************************************************************************
##
##   Methods to get ids from the pedigree that could be used as controls
##   for the given list of ids.
##
##***************************************************************************
## getAll
setMethod("getAll", "FAData",
          function(object, id=NULL, ...){
              return(doGetAll(pedigree(object), id=id, ...))
          })
setMethod("getExternalMatched", "FAData",
          function(object, id=NULL, match.using, ...){
              return(doGetExternalMatched(pedigree(object), id=id, match.using, ...))
          })
setMethod("getGenerationMatched", "FAData",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              return(doGetGenerationMatched(pedigree(object), id=id,
                                            include.anc=include.anc,
                                            include.off=include.off, ...))
          })
setMethod("getGenerationSexMatched", "FAData",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              return(doGetGenerationSexMatched(pedigree(object), id=id,
                                               include.anc=include.anc,
                                               include.off=include.off, ...))
          })
setMethod("getSexMatched", "FAData",
          function(object, id=NULL, ...){
              return(getSexMatched(pedigree(object), id=id, ...))
          })


##********************************************************************
##
##   Utility functions
##
##********************************************************************
sanitizeSex <- function(x){
    if(is.numeric(x)){
        ## require values 1, 2, NA
        numrange <- unique(x)
        if(!all(c(1, 2) %in% numrange))
            stop(paste0("If column sex is numeric it has to contain 1 (=male) ",
                        "and 2 (=female) values! Unknown or missing can be encoded",
                        " by NA or any other number than 1 or 2."))
        news <- rep(NA, length(x))
        news[which(x == 1)] <- "M"
        news[which(x == 2)] <- "F"
        return(factor(news, levels=c("M", "F")))
    }
    if(is.factor(x))
        x <- as.character(x)
    if(is.character(x)){
        ## remove all white spaces:
        x <- gsub(x, pattern="\\s", replacement="")
        ## that's tricky...kind of...
        ## just get the first character; want to have M and F...
        x <- substr(toupper(x), 1, 1)
        chars <- unique(x)
        if(length(chars) == 2){
            if(!all(c("M", "F") %in% chars))
                stop(paste0("Male and female individuals have to be represented",
                            " by 'M' and 'F', respectively!"))
        }
        if(length(chars) == 1){
            if(!any(c("M", "F") %in% chars))
                stop(paste0("Male and female individuals have to be represented",
                            " by 'M' and 'F', respectively!"))
        }
        if(length(chars[!is.na(chars)]) > 2)
            warning(paste0("All characters in column sex other than M and F",
                           " (or male, female etc) are set to NA!"))
        news <- rep(NA, length(x))
        news[which(x == "M")] <- "M"
        news[which(x == "F")] <- "F"
        return(factor(news, levels=c("M", "F")))
    }
}




