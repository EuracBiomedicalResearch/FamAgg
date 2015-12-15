## Generic method definitions.
setGeneric("affectedIndividuals", function(object)
    standardGeneric("affectedIndividuals"))

if(!isGeneric("affectedKinshipGroups"))
    setGeneric("affectedKinshipGroups", function(object)
        standardGeneric("affectedKinshipGroups"))
##if(!isGeneric("age"))
setGeneric("age", function(object)
    standardGeneric("age"))

##if(!isGeneric("age<-"))
setGeneric("age<-", function(object, value)
    standardGeneric("age<-"))

setGeneric("buildPed", function(object, id=NULL, ...)
    standardGeneric("buildPed"))

setGeneric("cliques", function(object, ...)
    standardGeneric("cliques"))

## setGeneric("cliques", function(object, ...)
##     standardGeneric("cliques"))
setGeneric("cliques<-", function(object, value)
    standardGeneric("cliques<-"))

setGeneric("cliqueAndTrait", function(object, ...)
    standardGeneric("cliqueAndTrait"))

setGeneric("countGenerations", function(object, id=NULL, ...)
    standardGeneric("countGenerations"))

setGeneric("estimateGenerations", function(object, family=NULL, ...)
    standardGeneric("estimateGenerations"))

if(!isGeneric("export"))
    setGeneric("export", function(object, con, format, ...)
               standardGeneric("export"))

setGeneric("familialIncidenceRate", function(object, trait=NULL, timeAtRisk=NULL,
                                             prune=TRUE, ...)
    standardGeneric("familialIncidenceRate"))

setGeneric("familialIncidenceRateTest", function(object, ...)
    standardGeneric("familialIncidenceRateTest"))

setGeneric("fsir", function(object, trait=NULL, lambda=NULL, timeInStrata=NULL,
                            prune=TRUE, ...)
    standardGeneric("fsir"))

setGeneric("fsirTest", function(object, ...)
    standardGeneric("fsirTest"))

setGeneric("genealogicalIndexTest", function(object, trait, nsim=50000, traitName,
                                             perFamilyTest=FALSE, controlSetMethod="getAll",
                                             prune=TRUE, strata=NULL, ...)
    standardGeneric("genealogicalIndexTest"))

setGeneric("findFounders", function(object, ...)
    standardGeneric("findFounders"))

setGeneric("generationsFrom", function(object, id=NULL, ...)
    standardGeneric("generationsFrom"))

setGeneric("getAncestors", function(object, id=NULL, ...)
    standardGeneric("getAncestors"))

setGeneric("getChildren", function(object, id=NULL, ...)
    standardGeneric("getChildren"))

setGeneric("getCommonAncestor", function(object, id=NULL, ...)
    standardGeneric("getCommonAncestor"))

setGeneric("getMissingMate", function(object, ...)
    standardGeneric("getMissingMate"))

setGeneric("getSiblings", function(object, id=NULL, ...)
    standardGeneric("getSiblings"))

setGeneric("kinship", function(id, ...)
    standardGeneric("kinship"))

if(!isGeneric("independentGroupCount"))
    setGeneric("independentGroupCount", function(object, ...)
        standardGeneric("independentGroupCount"))

setGeneric("kinshipGroupTest", function(object, ...)
    standardGeneric("kinshipGroupTest"))

setGeneric("kinshipSumTest", function(object, ...)
    standardGeneric("kinshipSumTest"))

setGeneric("lambda", function(object, ...)
    standardGeneric("lambda"))

setGeneric("oldResult", function(object, ...)
    standardGeneric("oldResult"))

##if(!isGeneric("pedigree"))
setGeneric("pedigree", function(object,...)
    standardGeneric("pedigree"))

##if(!isGeneric("pedigree<-"))
setGeneric("pedigree<-", function(object, value)
    standardGeneric("pedigree<-"))

setGeneric("pedigreeSize", function(object)
    standardGeneric("pedigreeSize"))

setGeneric("phenotypedIndividuals", function(object)
    standardGeneric("phenotypedIndividuals"))

setGeneric("plotPed", function(object, id=NULL, family=NULL, filename=NULL,
                               device="pdf", ...)
    standardGeneric("plotPed"))

setGeneric("plotRes", function(object, id=NULL, family=NULL, ...)
    standardGeneric("plotRes"))

setGeneric("probabilityTest", function(object, ...)
    standardGeneric("probabilityTest"))

setGeneric("shareKinship", function(object, ...)
    standardGeneric("shareKinship"))

setGeneric("result", function(object, ...)
    standardGeneric("result"))

setGeneric("resultForId", function(object, id=NULL)
    standardGeneric("resultForId"))

setGeneric("runSimulation", function(object, nsim, ...)
    standardGeneric("runSimulation"))

## setGeneric("subPedigree", function(ped, id=NULL, all=TRUE)
##     standardGeneric("subPedigree"))
setGeneric("timeAtRisk", function(object, ...)
    standardGeneric("timeAtRisk"))

setGeneric("timeAtRisk<-", function(object, value)
    standardGeneric("timeAtRisk<-"))

setGeneric("timeInStrata", function(object, ...)
    standardGeneric("timeInStrata"))

setGeneric("trait", function(object, ...)
    standardGeneric("trait"))

setGeneric("trait<-", function(object, value)
    standardGeneric("trait<-"))

setGeneric("traitByClique", function(object)
    standardGeneric("traitByClique"))

##***************************************************************************
##
## Methods to get matched controls.
##
##***************************************************************************
setGeneric("getAll", function(object, id=NULL, ...)
    standardGeneric("getAll"))
setGeneric("getExternalMatched", function(object, id=NULL, match.using, ...)
    standardGeneric("getExternalMatched"))
setGeneric("getGenerationMatched", function(object, id=NULL, include.anc=0,
                                            include.off=0, ...)
    standardGeneric("getGenerationMatched"))
setGeneric("getGenerationSexMatched", function(object, id=NULL, include.anc=0,
                                               include.off=0, ...)
    standardGeneric("getGenerationSexMatched"))
setGeneric("getSexMatched", function(object, id=NULL, ...)
    standardGeneric("getSexMatched"))


