data(minnbreast)
mbsub <- minnbreast[minnbreast$famid %in% c(4, 42, 165, 432),]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- FamAgg:::.PEDCN
PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
## data(minnbreast)
## mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
## colnames(PedDf) <- FamAgg:::.PEDCN
## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
## generate the FAData.
fad <- FAData(pedigree=PedDf)
## specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## OK, 0.0.6
## evaluating the kinship analysis method and for FAKinshipResult objects.
test_kinship_group <- function(){
    ##
    ## run the analysis.
    set.seed(18011977)
    far <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer",
                            nsim=500)
    ## repeating the test and evaluate whether we get the same results.
    set.seed(18011977)
    far2 <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer",
                        nsim=500)
    checkEquals(far, far2)

    ## replacing the trait with a random trait.
    tpreg <- mbsub$everpreg
    names(tpreg) <- mbsub$id
    ## randomizing to see if re-ordering works...
    tpregrand <- sample(tpreg, length(tpreg))
    trait(far) <- tpregrand
    checkEquals(trait(far), tpreg)

    ## testing to subset the object... which is not supported
    checkException(far[1:10, ])
    return(TRUE)
}

## question if for the plotting: how to build the pedigree: just the pheno and aff
## from the affectedKinshipGroups?
## plot the kinshipGroup.
test_plot_kinship <- function(){
    do.plot <- FALSE

    set.seed(18011977)
    far <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer", nsim=500)
    res <- result(far)
    res <- res[order(res$kinship_pvalue), ]

    if(do.plot){
        switchPlotfun("ks2paint")
        plotPed(far, id="17517", prune=TRUE, device="plot", cex=0.4)
        plotPed(far, id="17517", prune=FALSE, device="plot")
    }
    ## switchPlotfun()
}


test_buildped_kinship <- function(){
    ## data(minnbreast)
    ## mbsub <- minnbreast[minnbreast$famid %in% c(4, 42, 165, 432),]
    ## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(PedDf) <- FamAgg:::.PEDCN
    ## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    ## ## generate the FAData.
    ## fad <- FAData(pedigree=PedDf)
    ## ## specify the trait.
    ## tcancer <- mbsub$cancer
    ## names(tcancer) <- mbsub$id

    set.seed(18011977)
    far <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer", nsim=500)
    res <- result(far)
    res <- res[order(res$kinship_pvalue), ]

    FullPed <- buildPed(far, id=res[1, "group_id"])
    affGroup <- affectedKinshipGroups(far)[[res[1, "group_id"]]]
    allids <- unique(c(affGroup$aff, affGroup$pheno))
    ## check if all ids are in the pedigree.
    checkEquals(length(allids), sum(allids %in% FullPed$id))
    ## check if the same works for the pruned pedigree:
    SmallPed <- buildPed(far, id=res[1, "group_id"], prune=TRUE)
    checkEquals(length(allids), sum(allids %in% SmallPed$id))
}

notrun_test_compare_strat_normal <- function(){
    ## test the strata test and compare that to the "normal" test
    set.seed(18011977)
    fag <- kinshipGroupTest(fad, trait=tcancer, nsim=500)
    fag2 <- fag
    set.seed(18011977)
    fag2 <- runSimulation(fag2, nsim=500)
    ## Compare
    checkEquals(fag@sim$pvalueRatio, fag2@sim$pvalueRatio)
    checkEquals(fag@sim$pvalueKinship, fag2@sim$pvalueKinship)
    ## now do it with strata = 1.
    strat <- rep("A", length(fad$id))
    set.seed(18011977)
    fag2 <- FamAgg:::runSimulationStrata(fag, nsim=500, strata=strat)
    checkEquals(fag@sim$pvalueRatio, fag2@sim$pvalueRatio)
    checkEquals(fag@sim$pvalueKinship, fag2@sim$pvalueKinship)
    ## Test if the NAs are correctly removed...
    strat[!is.na(tcancer)][1:10] <- NA
    fag2 <- FamAgg:::runSimulationStrata(fag, nsim=500, strata=strat)
    ## Have to get the message that 10 phenotyped are removed.
    ## OK, that's it. Now do the same but with sex strata.
    fasex <- FamAgg:::runSimulationStrata(fag, nsim=500, strata=fad$sex)
    Col <- rep("blue", length(fad$id))
    Col[which(fad$sex == "F")] <- "red"
    plot(fasex@sim$pvalueKinship, fag@sim$pvalueKinship, xlab="Sex stratified",
         ylab="No stratification", col=Col)
    abline(0, 1)
}


