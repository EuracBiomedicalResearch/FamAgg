## Generic method definitions.
## A
setGeneric("affectedIndividuals", function(object)
    standardGeneric("affectedIndividuals"))
if(!isGeneric("affectedKinshipGroups"))
    setGeneric("affectedKinshipGroups", function(object)
        standardGeneric("affectedKinshipGroups"))
setGeneric("age", function(object)
    standardGeneric("age"))
setGeneric("age<-", function(object, value)
    standardGeneric("age<-"))

## B
setGeneric("buildPed", function(object, ...)
    standardGeneric("buildPed"))

## C
setGeneric("cliques", function(object, ...)
    standardGeneric("cliques"))
setGeneric("cliques<-", function(object, value)
    standardGeneric("cliques<-"))
setGeneric("cliqueAndTrait", function(object, ...)
    standardGeneric("cliqueAndTrait"))
setGeneric("countGenerations", function(object, id=NULL, ...)
    standardGeneric("countGenerations"))

## D

## E
setGeneric("estimateGenerations", function(object, family=NULL, ...)
    standardGeneric("estimateGenerations"))
if(!isGeneric("export"))
    setGeneric("export", function(object, con, format, ...)
               standardGeneric("export"))

## F
setGeneric("familialIncidenceRate", function(object, trait=NULL,
                                             timeAtRisk=NULL, ...)
    standardGeneric("familialIncidenceRate"))
setGeneric("familialIncidenceRateTest", function(object, ...)
    standardGeneric("familialIncidenceRateTest"))
setGeneric("findFounders", function(object, ...)
    standardGeneric("findFounders"))
setGeneric("fsir", function(object, trait=NULL, lambda=NULL,
                            timeInStrata=NULL, ...)
    standardGeneric("fsir"))
setGeneric("fsirTest", function(object, ...)
    standardGeneric("fsirTest"))

## G
setGeneric("genealogicalIndexTest", function(object, trait, nsim=50000, traitName,
                                             perFamilyTest=FALSE,
                                             controlSetMethod="getAll",
                                             rm.singletons=TRUE, strata=NULL, ...)
    standardGeneric("genealogicalIndexTest"))
setGeneric("generationsFrom", function(object, id=NULL, ...)
    standardGeneric("generationsFrom"))
setGeneric("getAncestors", function(object, id=NULL, ...)
    standardGeneric("getAncestors"))
setGeneric("getChildren", function(object, id=NULL, ...)
    standardGeneric("getChildren"))
setGeneric("getCommonAncestor", function(object, id=NULL, ...)
    standardGeneric("getCommonAncestor"))
setGeneric("getFounders", function(object, ...)
    standardGeneric("getFounders"))
setGeneric("getMissingMate", function(object, ...)
    standardGeneric("getMissingMate"))
setGeneric("getSiblings", function(object, id=NULL, ...)
    standardGeneric("getSiblings"))
setGeneric("getSingletons", function(object, ...)
    standardGeneric("getSingletons"))
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

## H

## I
setGeneric("independentGroupCount", function(object, ...)
    standardGeneric("independentGroupCount"))

## J

## K
setGeneric("kinship", function(id, ...)
    standardGeneric("kinship"))
setGeneric("kinshipGroupTest", function(object, ...)
    standardGeneric("kinshipGroupTest"))
setGeneric("kinshipSumTest", function(object, ...)
    standardGeneric("kinshipSumTest"))

## L
setGeneric("lambda", function(object, ...)
    standardGeneric("lambda"))

## M

## N

## O
setGeneric("oldResult", function(object, ...)
    standardGeneric("oldResult"))

## P
setGeneric("pedigree", function(object,...)
    standardGeneric("pedigree"))
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

## Q

## R
setGeneric("removeSingletons", function(object, ...)
           standardGeneric("removeSingletons"))
setGeneric("result", function(object, ...)
    standardGeneric("result"))
setGeneric("resultForId", function(object, id=NULL)
    standardGeneric("resultForId"))
setGeneric("runSimulation", function(object, nsim, ...)
    standardGeneric("runSimulation"))

## S
setGeneric("shareKinship", function(object, ...)
    standardGeneric("shareKinship"))

## T
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

## U
## V
## W
## X
## Y
## Z


