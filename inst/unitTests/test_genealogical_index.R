## load the data we're going to analyze
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid %in% c(4, 5, 6, 8, 14, 16), ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- FamAgg:::.PEDCN
PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
## generate the FAData.
fad <- FAData(pedigree=PedDf)
## specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## per family tests: that way I can better control the control methods.
test_gen_index_per_family <- function(){
    do.plot <- FALSE

    ## perform the test:
    set.seed(18011977)
    ## family 14 looks promising!
    trait(fad) <- tcancer
    ##plotPed(far, family="14", cex=0.5)

    gi <- as(fad, "FAGenIndexResults")
    Test <- runSimulation(gi, controlSetMethod="getAll", nsim=100, perFamilyTest=TRUE)
    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSetMethod="getAll", nsim=100,
    ##                            perFamilyTest=TRUE)
    All <- Test@sim[["14"]]
    All8 <- Test@sim[["8"]]
    expectAff <- as.character(c(444, 447, 450, 451, 452))
    checkEquals(sort(All$affected), sort(expectAff))

    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSetMethod="getSexMatched", nsim=100,
    ##                            perFamilyTest=TRUE)
    ## SexMatch <- Test[["14"]]
    ## SexMatch8 <- Test[["8"]]
    Test <- runSimulation(gi, controlSetMethod="getSexMatched", nsim=100,
                          perFamilyTest=TRUE)
    SexMatch <- Test@sim[["14"]]
    SexMatch8 <- Test@sim[["8"]]
    expectCtrls <- as.character(c(441, 444, 446, 447, 442, 443, 462, 450, 451, 452,
                                  453, 454, 468, 469))
    checkEquals(sort(SexMatch$ctrls), sort(expectCtrls))

    ## Test generation matched.
    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSetMethod="getGenerationMatched", nsim=100,
    ##                            perFamilyTest=TRUE)
    ## GenMatch <- Test[["14"]]
    Test <- runSimulation(gi, controlSetMethod="getGenerationMatched", nsim=100,
                          perFamilyTest=TRUE)
    GenMatch <- Test@sim[["14"]]
    expectCtrls <- as.character(c(459, 444, 445, 460, 446, 461, 447, 442, 458, 443,
                                  448, 449, 462, 450, 464, 451, 465, 452, 466, 453,
                                  467, 454))
    checkEquals(sort(GenMatch$ctrls), sort(expectCtrls))

    ## Test generation & sex matched. Note that the result is no longer significant!
    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSetMethod="getGenerationSexMatched", nsim=100,
    ##                            perFamilyTest=TRUE)
    ## GenSexMatch <- Test[["14"]]
    Test <- runSimulation(gi, controlSetMethod="getGenerationSexMatched", nsim=100,
                          perFamilyTest=TRUE)
    GenSexMatch <- Test@sim[["14"]]
    expectCtrls <- as.character(c(444, 446, 447, 442, 443, 462, 450, 451, 452, 453,
                                  454))
    checkEquals(sort(GenSexMatch$ctrls), sort(expectCtrls))

    if(do.plot){
        par(mfrow=c(4, 1))
        plot(All$expDensity, main="getAll", type="h", lwd=2, col="lightgrey")
        points(All$expDensity, type="l", col="grey")
        abline(v=All$meanKinship)
        plot(SexMatch$expDensity, "getSexMatched", type="h", lwd=2, col="lightgrey")
        points(SexMatch$expDensity, type="l", col="grey")
        abline(v=SexMatch$meanKinship)
        plot(GenMatch$expDensity, "getGenerationMatched", type="h", lwd=2, col="lightgrey")
        points(GenMatch$expDensity, type="l", col="grey")
        abline(v=GenMatch$meanKinship)
        plot(GenSexMatch$expDensity, "getGenerationSexMatched", type="h", lwd=2, col="lightgrey")
        points(GenSexMatch$expDensity, type="l", col="grey")
        abline(v=GenSexMatch$meanKinship)
    }

    ## families with male cases: 4, 8, 9, 10, 19
    ## other interesting family: 8
    ## LLLL Test with or without strata sampling for family 8.
    plotPed(fad, family="8", cex=0.8)
    ## sex matching and getting all should not make any difference here.
    checkEquals(sort(All8$ctrls), sort(SexMatch8$ctrls))

    ## do the test again using sex as strata:
    set.seed(18011977)
    Test <- runSimulation(gi, controlSetMethod="getAll", nsim=100, perFamilyTest=TRUE,
                          strata=fad$sex)
    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSeteMethod="getAll", nsim=100, perFamilyTest=TRUE,
    ##                            strata=fad$sex)
    ## All8Strata <- Test[["8"]]
    All8Strata <- Test@sim[["8"]]
    checkTrue(All8Strata$pvalueKinship != All8$pvalueKinship)

    ## the results for family 14 have to be identical:
    Test <- runSimulation(gi, controlSetMethod="getSexMatched", nsim=100, perFamilyTest=TRUE,
                          strata=gi$sex)
    ## Test <- FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
    ##                            controlSeteMethod="getSexMatched", nsim=100, perFamilyTest=TRUE,
    ##                            strata=fad$sex)
    checkEquals(SexMatch$pvalueKinship, Test@sim[["14"]]$pvalueKinship)
}

## The same tests, but once for the full pedigree.
test_gen_index_pedigree <- function(){
    set.seed(18011977)
    ## use getAll.
    giAll <- genealogicalIndexTest(fad, trait=tcancer, traitName="cancer", nsim=100,
                                   controlSetMethod="getAll")
    checkTrue(result(giAll)[1, "pvalue"] < 0.05)

    ## use getSexMatched
    giSexM <- genealogicalIndexTest(fad, trait=tcancer, traitName="cancer", nsim=100,
                                    controlSetMethod="getSexMatched")
    ## it is NOT the same as getAll, since some individuals with unknown sex were
    ## removed.

    ## use getSexMatched, strata.
    giSexMStrata <- genealogicalIndexTest(fad, trait=tcancer, traitName="cancer", nsim=100,
                                          controlSetMethod="getSexMatched", strata=fad$sex)

    ## checkError for getGenerationMatched
    checkException(FamAgg:::.genIndex(pedigree(fad), kin=kinship(fad), trait=tcancer,
                                      controlSetMethod="getGenerationMatched", nsim=100))
}


test_gen_index_with_gap <- function(){
    library(gap)
    gi <- genealogicalIndexTest(fad, trait=tcancer, nsim=100, rm.singletons=FALSE)

    ## Now use gap:
    gapPed <- pedigree(fad)
    gapPed[is.na(gapPed$mother), "mother"] <- 0
    gapPed[is.na(gapPed$father), "father"] <- 0
    ## get the ids of the affected
    affIds <- names(tcancer[which(tcancer == 1)])
    giGap <- gif(gapPed[, c("id", "father", "mother")], as.numeric(affIds))
    checkEquals(result(gi)$genealogical_index, giGap[[1]])
    ## Doing the same with removal of singletons:
    gi <- genealogicalIndexTest(fad, trait=tcancer, nsim=100, rm.singletons=TRUE)
    gapPed <- pedigree(fad)
    gapPed <- removeSingletons(gapPed)
    gapPed[is.na(gapPed$mother), "mother"] <- 0
    gapPed[is.na(gapPed$father), "father"] <- 0
    affIds <- affIds[affIds %in% gapPed$id]
    giGap <- gif(gapPed[, c("id", "father", "mother")], as.numeric(affIds))
    checkEquals(result(gi)$genealogical_index, giGap[[1]])
}

