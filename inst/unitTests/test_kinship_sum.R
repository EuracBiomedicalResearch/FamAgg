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

## OK, 0.0.6
## testing the kinship sum test.
test_kinship_sum <- function(){
    do.plot <- FALSE

    ## perform the test:
    set.seed(18011977)
    far <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                          nsim=1000)
    Result <- result(far)
    dim(Result)
    ## check the pedigree/family of this individual
    rownames(PedDf) <- as.character(PedDf$id)
    ## use the family method to extract the pedigree of that family.
    GotFam <- family(far, id=Result[1, "affected_id"])
    AffFam <- PedDf[PedDf$family == PedDf[Result[1, "affected_id"], "family"], ]
    AffFam <- FamAgg:::sanitizePed(AffFam)
    AffFam <- cbind(AffFam, affected=tcancer[AffFam$id])
    checkEquals(GotFam, AffFam)

    if(do.plot){
        plotPed(far, id=Result[1, "affected_id"])
    }

    ## now replacing the trait and re-running.
    tpreg <- mbsub$everpreg
    names(tpreg) <- mbsub$id
    ## randomizing to see if re-ordering works...
    tpregrand <- sample(tpreg, length(tpreg))
    trait(far) <- tpregrand
    checkEquals(trait(far), tpreg)
    ## running the simulation.
    far <- runSimulation(far, nsim=1000)

    ## testing to subset the object... which is not supported
    checkException(far[1:10, ])

    ##
    return(TRUE)
}

test_kinship_sum_strata <- function(){
    set.seed(18011977)
    far <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                          nsim=500)
    ## do the same using a fake strata
    set.seed(18011977)
    farStrat <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                               nsim=500, strata=rep("OK", length(tcancer)))
    checkEquals(far@sim$sumKinship, farStrat@sim$sumKinship)
    checkEquals(far@sim$pvalueKinship, farStrat@sim$pvalueKinship)
    ## check again, subsetting to male only.
    fadM <- fad[which(fad$sex == "M"), ]
    set.seed(18011977)
    farM <- kinshipSumTest(fadM, trait=tcancer, nsim=500)
    sexStrat <- fad$sex
    sexStrat[sexStrat == "F"] <- NA
    set.seed(18011977)
    farStrat <- kinshipSumTest(fad, trait=tcancer, nsim=500, strata=sexStrat)
    ## compare that
    checkEquals(farM@sim$sumKinship, farStrat@sim$sumKinship)
    checkEquals(farM@sim$pvalueKinship, farStrat@sim$pvalueKinship)
}

test_plot_kinclust <- function(){
    ## data(minnbreast)
    ## mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
    ## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(PedDf) <- FamAgg:::.PEDCN
    ## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    ## ## generate the FAData.
    ## fad <- FAData(pedigree=PedDf)
    ## ## specify the trait.
    ## tcancer <- mbsub$cancer
    ## names(tcancer) <- mbsub$id

    ## perform the test:
    set.seed(18011977)
    far <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                          nsim=500)
    res <- result(far)

    id <- "410"

    cols <- rep(1, 32)
    cols[1] <- 2

    ## kinship2 plotting
    switchPlotfun("ks2paint")
    plotPed(far, id=id, device="plot", col=cols)
    ##switchPlotfun()
}


