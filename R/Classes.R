## Definition of all classes.

## the column names for the pedigree
.PEDCN <- c("family", "id", "father", "mother", "sex")

##************************************************
##
##       FAData
##
##       Main class defining/containing all required data.
##       *
##************************************************
setClass("FAData",
         slots=c(
             pedigree="data.frame",
             age="numeric",
             .kinship="dsCMatrix",
             .trait="numeric",
             traitname="character"
         ),
         prototype=list(
             pedigree = data.frame(family = character(),
                                   id = character(),
                                   father = character(),
                                   mother = character(),
                                   sex = factor(levels = c("F", "M"))),
             age=numeric(),
             .kinship=new("dsCMatrix"),
             .trait=numeric(),
             traitname=character()
         )
         )

##************************************************
##
##       FAResult
##
##       main result class...
## setClass("FATrait",
##          contains="FAData",
##          slots=c(
##              ##             "VIRTUAL",
##              .trait="numeric",
##              traitname="character"
##          ),
##          prototype=list(.trait=numeric(),
##                         traitname=character()
##                         )
##          )

##************************************************
##
##       FAKinGroupResult
##
##       sim: the results from the simulation.
##       nsim: the number of simulations.
##       indGroupCount: number of independent affected kinship groups.
##       affectedKinGroups: list of kinship groups of each affected individual.
##************************************************
setClass("FAKinGroupResults",
         contains="FAData",
         slots=c(
             sim="list",
             nsim="numeric",
             affectedKinshipGroups="list"
         ),
         prototype=list(
             sim=list(),
             nsim=0,
             affectedKinshipGroups=list()
         )
         )

##************************************************
##
##       FAKinSumResult
##
##       sim: the results from the simulation.
##       nsim: the number of simulations.
##************************************************
setClass("FAKinSumResults",
         contains="FAData",
         slots=c(
             sim="list",
             nsim="numeric"
         ),
         prototype=list(
             sim=list(),
             nsim=0
         )
         )


##************************************************
##
##       FAGenIndexResults
##
##       sim: the results from the simulation.
##       nsim: the number of simulations.
##************************************************
setClass("FAGenIndexResults",
         contains="FAData",
         slots=c(
             sim="list",
             nsim="numeric",
             controlSetMethod="character",
             perFamilyTest="logical"
         ),
         prototype=list(
             sim=list(),
             nsim=0,
             controlSetMethod="getAll",
             perFamilyTest=FALSE
         )
         )


##************************************************
##
##       FAIncidenceRateResults
##
##       sim: the results from the simulation.
##       nsim: the number of simulations.
##************************************************
setClass("FAIncidenceRateResults",
         contains="FAData",
         slots=c(
             sim="list",
             nsim="numeric",
             timeAtRisk="numeric"
         ),
         prototype=list(
             sim=list(),
             nsim=0,
             timeAtRisk=numeric()
         )
         )

##************************************************
##
##       FAStdIncidenceRateResults
##
##       sim: the results from the simulation.
##       nsim: the number of simulations.
##************************************************
setClass("FAStdIncidenceRateResults",
         contains="FAData",
         slots=c(
             sim="list",
             nsim="numeric",
             timeInStrata="matrix",
             lambda="numeric"
         ),
         prototype=list(
             sim=list(),
             nsim=0,
             timeInStrata=matrix(),
             lambda=numeric()
         )
         )

setClass("FABinTestResults",
         contains = "FAData",
         slots = c(result = "data.frame")
        ,
         prototype = list(
             result = data.frame(total_phenotyped = integer(),
                                 total_affected = integer(),
                                 family = character(),
                                 phenotyped = integer(),
                                 affected = integer(),
                                 pvalue = numeric(),
                                 prob = numeric(),
                                 check.names = FALSE, stringsAsFactors = FALSE)
         )
         )

