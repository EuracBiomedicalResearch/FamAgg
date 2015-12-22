## incidence ratio test as defined by Kerber 1995
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- FamAgg:::.PEDCN
PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
## generate the FAData.
fad <- FAData(pedigree=PedDf)
## specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id
## the time at risk...
tar <- mbsub$endage


test_incidence_ratio <- function(){

    fr <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                       timeAtRisk=tar)
    ## the same but per family
    frFam <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                          timeAtRisk=tar, perFamilyTest=TRUE)
    fr3 <- familialIncidenceRate(fad, trait=tcancer, timeAtRisk=tar)
    checkEquals(fr[[1]], fr3)
    ## The argument rm.singletons has been removed, so no need for the test below.
    ## Using the method: check if rm.singletons has an effect.
    ## fr <- familialIncidenceRate(fad, trait=tcancer, timeAtRisk=tar, rm.singletons=FALSE)
    ## fr2 <- familialIncidenceRate(fad, trait=tcancer, timeAtRisk=tar, rm.singletons=TRUE)
    ## checkEquals(sum(fr, na.rm=TRUE)==sum(fr2, na.rm=TRUE), TRUE)
    ## Note: we're getting the same results here as we are summing kinship coefs. singletons have
    ## a kinship of 0 with any other, thus they don't add to the value anyway.

    ## So, it's not the same if we run this within family or for the whole pedigree...
    ## frVals <- unlist(fr, use.names=FALSE)
    ## frFamVals <- unlist(frFam, use.names=FALSE)
    ## checkEquals(frVals[!is.na(frVals)], frFamVals[!is.na(frFamVals)])
    ## checkEquals(unlist(fr, use.names=FALSE), unlist(frFam, use.names=FALSE))
    ## fr[[1]][1:10]
    ## frFam[["4"]]
}

test_estimate_time_at_risk <- function(){
    sdates <- c("2012-04-17", "2014-05-29", "1999-12-31", "2002-10-10")
    edates <- c("2015-09-15", "2015-09-15", "2005-09-15", "2015-09-15")
    idates <- c(NA, NA, "2007-07-13", "2013-12-23")
    ddates <- c(NA, NA, NA, "2014-03-14")

    checkException(estimateTimeAtRisk(sdates))
    ests <- estimateTimeAtRisk(startDate=sdates, endDate=edates)
    ests2 <- estimateTimeAtRisk(startDate=sdates, endDate=edates,
                                incidenceDate=idates)
    ests3 <- estimateTimeAtRisk(startDate=sdates, endDate=edates,
                                incidenceDate=idates, deathDate=ddates)
}

test_fr_simulation <- function(){
    doPlot <- FALSE
    fr <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                       timeAtRisk=tar)
    fr <- fr[[1]]
    fr2 <- familialIncidenceRate(fad, trait=tcancer,
                                 timeAtRisk=tar)
    checkEquals(fr, fr2[names(fr)])
    ## Use the familialIncidenceRateTest
    set.seed(18011977)
    Test <- FamAgg:::.FRSimulation(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                                   timeAtRisk=tar, prune=TRUE, nsim=1000)
    set.seed(18011977)
    frRes <- familialIncidenceRateTest(fad, trait=tcancer,
                                       timeAtRisk=tar, nsim=1000)
    checkEquals(frRes@sim$fir, fr2)
    ## Do the simulation using dummy strata.
    set.seed(18011977)
    frStrat <- familialIncidenceRateTest(fad, trait=tcancer,
                                                  timeAtRisk=tar, nsim=1000,
                                                  strata=rep(1, length(fad$id)))
    checkEquals(frRes@sim$fir, frStrat@sim$fir)
    checkEquals(frRes@sim$pvalue, frStrat@sim$pvalue)
    ## Repeat using the low mem version.
    set.seed(18011977)
    frLM <- familialIncidenceRateTest(fad, trait=tcancer,
                                               timeAtRisk=tar, nsim=1000,
                                               lowMem=TRUE)
    checkEquals(frRes@sim$fir, frLM@sim$fir)
    checkEquals(frRes@sim$pvalue, frLM@sim$pvalue)
    ##
    ## plotting...
    if(doPlot){
        res <- result(frRes)
        plotPed(frRes, id="4")
        plotPed(frRes, family=19)
        plotPed(frRes, id=res[1, "id"])
        ## plotRes.
        plotRes(frRes, id="4")
        plotRes(frRes, id=res[1, "id"])
    }

    ## Testing the $ accessors.
    checkEquals(frRes$fir, frRes@sim$fir)
    checkEquals(frRes$tar, FamAgg:::timeAtRisk(frRes))

    ## Testing with and without rm.singletons.
    fr <- familialIncidenceRateTest(fad, trait=tcancer, nsim=400,
                                    timeAtRisk=tar, rm.singletons=TRUE)
    fr2 <- familialIncidenceRateTest(fad, trait=tcancer, nsim=400,
                                     timeAtRisk=tar, rm.singletons=FALSE)
    checkEquals(fr$fir, fr2$fir)
}


