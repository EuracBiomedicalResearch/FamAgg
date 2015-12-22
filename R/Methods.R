## some other methods and methods for data.frame, pedigree or pedigreeList
## countGenerations
setMethod("countGenerations", "data.frame",
          function(object, id=NULL, direction="down"){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doCountGenerations(object, id=id, direction=direction))
          })
setMethod("countGenerations", "pedigreeList",
          function(object, id=NULL, direction="down"){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doCountGenerations(object, id=id, direction=direction))
          })
setMethod("countGenerations", "pedigree",
          function(object, id=NULL, direction="down"){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doCountGenerations(object, id=id, direction=direction))
          })

## estimateGenerations
setMethod("estimateGenerations", "data.frame",
          function(object, family=NULL, ...){
              object <- sanitizePed(object)
              object <- checkPedCol(object)
              return(doEstimateGenerationsFor2(object, family=family))
          })
setMethod("estimateGenerations", "pedigree",
          function(object, family=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doEstimateGenerationsFor2(object, family=family))
          })
setMethod("estimateGenerations", "pedigreeList",
          function(object, family=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doEstimateGenerationsFor2(object, family=family))
          })


## findFounders
setMethod("findFounders", "data.frame",
          function(object, family=NULL, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doFindFounders(object, family=family))
          })
setMethod("findFounders", "pedigreeList",
          function(object, family=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doFindFounders(object, family=family))
          })
setMethod("findFounders", "pedigreeList",
          function(object, family=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doFindFounders(object, family=family))
          })

## generationsFrom
setMethod("generationsFrom", "data.frame",
          function(object, id=NULL, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationFrom2(object, id=id, ...))
          })
setMethod("generationsFrom", "pedigree",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationFrom2(object, id=id, ...))
          })
setMethod("generationsFrom", "pedigreeList",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationFrom2(object, id=id, ...))
          })


## getAncestors
setMethod("getAncestors", "data.frame",
          function(object, id=NULL, max.generations=3, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetAncestors(ped=object, id=id, maxlevel=max.generations, ...))
          })
setMethod("getAncestors", "pedigreeList",
          function(object, id=NULL, max.generations=3, ...){
              if(!is(object, "pedigree") & !is(object, "pedigreeList"))
                  stop("object should be either a 'pedigree' or 'pedigreeList' object!")
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetAncestors(ped=object, id=id, maxlevel=max.generations, ...))
          })
setMethod("getAncestors", "pedigree",
          function(object, id=NULL, max.generations=3, ...){
              if(!is(object, "pedigree") & !is(object, "pedigreeList"))
                  stop("object should be either a 'pedigree' or 'pedigreeList' object!")
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetAncestors(ped=object, id=id, maxlevel=max.generations, ...))
          })

## getChildren
setMethod("getChildren", "data.frame",
          function(object, id=NULL, max.generations=3, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetChildren(ped=object, id=id, maxlevel=max.generations, ...))
          })
setMethod("getChildren", "pedigreeList",
          function(object, id=NULL, max.generations=3, ...){
              if(!is(object, "pedigree") & !is(object, "pedigreeList"))
                  stop("object should be either a 'pedigree' or 'pedigreeList' object!")
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetChildren(ped=object, id=id, maxlevel=max.generations, ...))
          })
setMethod("getChildren", "pedigree",
          function(object, id=NULL, max.generations=3, ...){
              if(!is(object, "pedigree") & !is(object, "pedigreeList"))
                  stop("object should be either a 'pedigree' or 'pedigreeList' object!")
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetChildren(ped=object, id=id, maxlevel=max.generations, ...))
          })

## getCommonAncestor
setMethod("getCommonAncestor", "data.frame",
          function(object, id=NULL, method="min.dist"){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetCommonAncestor(object, id=id, method=method))
          })
setMethod("getCommonAncestor", "pedigree",
          function(object, id=NULL, method="min.dist"){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetCommonAncestor(object, id=id, method=method))
          })
setMethod("getCommonAncestor", "pedigreeList",
          function(object, id=NULL, method="min.dist"){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetCommonAncestor(object, id=id, method=method))
          })

## getMissingMate
setMethod("getMissingMate", "data.frame",
          function(object, id=NULL, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetMissingMate(object, id=id))
          })
setMethod("getMissingMate", "pedigreeList",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetMissingMate(object, id=id))
          })
setMethod("getMissingMate", "pedigree",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetMissingMate(object, id=id))
          })

## getSiblings
setMethod("getSiblings", "data.frame",
          function(object, id=NULL, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSiblings(object, id=id))
          })
setMethod("getSiblings", "pedigreeList",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSiblings(object, id=id))
          })
setMethod("getSiblings", "pedigree",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSiblings(object, id=id))
          })




##***************************************************************************
##
##   Methods to get ids from the pedigree that could be used as controls
##   for the given list of ids.
##
##   Input: ped, id and any number of optional parameters as input
##   Return: a list of ids, names of the list are the family ids.
##
##***************************************************************************
## getAll
setMethod("getAll", "data.frame",
          function(object, id=NULL, ...){
              object <- checkPedCol(object)
              return(doGetAll(ped=object, id=id, ...))
          })
setMethod("getAll", "pedigreeList",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              return(doGetAll(ped=object, id=id, ...))
          })
setMethod("getAll", "pedigree",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              return(doGetAll(ped=object, id=id, ...))
          })

## getExternalMatched
setMethod("getExternalMatched", "data.frame",
          function(object, id=NULL, match.using, ...){
              object <- checkPedCol(object)
              return(doGetExternalMatched(ped=object, id=id, match.using, ...))
          })
setMethod("getExternalMatched", "pedigreeList",
          function(object, id=NULL, match.using, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              return(doGetExternalMatched(ped=object, id=id, match.using, ...))
          })
setMethod("getExternalMatched", "pedigree",
          function(object, id=NULL, match.using, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              return(doGetExternalMatched(ped=object, id=id, match.using, ...))
          })

## getGenerationMatched
setMethod("getGenerationMatched", "data.frame",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationMatched(ped=object, id=id, include.anc=include.anc,
                                            include.off=include.off, ...))
          })
setMethod("getGenerationMatched", "pedigreeList",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationMatched(ped=object, id=id, include.anc=include.anc,
                                            include.off=include.off, ...))
          })
setMethod("getGenerationMatched", "pedigree",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationMatched(ped=object, id=id, include.anc=include.anc,
                                            include.off=include.off, ...))
          })

## getGenerationSexMatched
setMethod("getGenerationSexMatched", "data.frame",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationSexMatched(ped=object, id=id, include.anc=include.anc,
                                               include.off=include.off, ...))
          })
setMethod("getGenerationSexMatched", "pedigreeList",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationSexMatched(ped=object, id=id, include.anc=include.anc,
                                               include.off=include.off, ...))
          })
setMethod("getGenerationSexMatched", "pedigree",
          function(object, id=NULL, include.anc=0, include.off=0, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetGenerationSexMatched(ped=object, id=id, include.anc=include.anc,
                                               include.off=include.off, ...))
          })

## getSexMatched
setMethod("getSexMatched", "data.frame",
          function(object, id=NULL, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSexMatched(ped=object, id=id, ...))
          })
setMethod("getSexMatched", "pedigreeList",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSexMatched(ped=object, id=id, ...))
          })
setMethod("getSexMatched", "pedigree",
          function(object, id=NULL, ...){
              object <- ped2df(object)
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSexMatched(ped=object, id=id, ...))
          })

## getFounders
setMethod("getFounders", "data.frame",
          function(object, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetFounders(object, ...))
          })

## getSingletons
setMethod("getSingletons", "data.frame",
          function(object, ...){
              object <- checkPedCol(object)
              object <- sanitizePed(object)
              return(doGetSingletons(object, ...))
          })

