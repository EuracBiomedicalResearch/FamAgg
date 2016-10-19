.getMinSub <- function(){
    data(minnbreast)
    mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
    return(mbsub)
}

mbsub <- .getMinSub()

test_show <- function() {
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN
    fad <- FAData(pedigree=mbped)
    res <- capture.output(show(fad))
    checkEquals(res[1], "FAData object with:")
}

test_construct <- function(){
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN
    fad <- FAData(pedigree=mbped)
    tcancer <- mbsub$cancer
    names(tcancer) <- mbsub$id
    trait(fad) <- tcancer
    ## constructor with the trait
    fad <- FAData(pedigree=mbped, trait=mbsub$cancer)

    tcancershuffle <- sample(tcancer, length(tcancer))
    trait(fad) <- tcancershuffle
    ## the function internally has to rematch the ids
    checkEquals(trait(fad), tcancer)
    ## testing the $ accessor
    checkEquals(fad$trait, tcancer)
    fad@traitname <- "cancer"
    checkEquals(fad$cancer, tcancer)

    ## check if the pedigree returns the correct extended pedigree
    ## with the additional column affected
    ped <- pedigree(fad)
    fromPed <- ped$affected
    names(fromPed) <- ped$id
    checkEquals(fromPed, tcancer)

    ## test affected and phenotyped individuals
    affIds <- affectedIndividuals(fad)
    checkEquals(affIds, as.character(ped$id[which(ped$affected == 1)]))
    phenoIds <- phenotypedIndividuals(fad)
    checkEquals(phenoIds, as.character(ped$id[!is.na(ped$affected)]))

}

test_FAData_from_file <- function(){
    ## create from file
    pedfile <- system.file("txt/minnbreastsub.txt", package="FamAgg")
    checkException(fad2 <- FAData(pedfile, id.col="ID" ))
    fad2 <- FAData(pedfile, id.col="id", family.col="famid", father.col="fatherid",
                   mother.col="motherid", sex.col="sex", header=TRUE)
    ## Now read it from ped
    pedfile <- system.file("txt/minnbreastsub.ped.gz", package="FamAgg")
    fadPed <- FAData(pedfile)
    pedfile <- system.file("txt/minnbreastsub-bad.ped.gz", package="FamAgg")
    checkException(FAData(pedfile))
}


## Testing the validation methods
test_validate <- function(){
    tmpPed <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(tmpPed) <- FamAgg:::.PEDCN
    ## Don't allow NAs in id, father or mother
    badPed <- tmpPed
    badPed[13, "id"] <- NA
    checkException(FAData(badPed))
    badPed <- tmpPed
    badPed[13, "father"] <- NA
    checkException(FAData(badPed))
    badPed <- tmpPed
    badPed[13, "mother"] <- NA
    checkException(FAData(badPed))
    ## Don't allow father or mother ids other than 0 or %in% $id
    badPed <- tmpPed
    badPed[13, "mother"] <- "123345"
    checkException(FAData(badPed))
    badPed <- tmpPed
    badPed[13, "father"] <- "123345"
    checkException(FAData(badPed))
}


## OK, 0.0.6
test_age <- function(){
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN

    tmp <- new("FAData")
    ages <- mbsub$endage
    names(ages) <- mbped$id
    age(tmp) <- ages
    checkEqualsNumeric(length(age(tmp)), length(ages))
    ## pedigreee
    pedigree(tmp) <- mbped
    ## once the pedigree is set, the age accessor should return the
    ## ages in the same order than ids in the pedigree!
    checkEqualsNumeric(length(age(tmp)), nrow(mbped))
    checkEquals(names(age(tmp)), as.character(mbped[, 2]))
}

test_pedigree <- function(){
    tmp <- new("FAData")

    ## mbsub <- .getMinSub()

    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    pedigree(tmp) <- PedDf
    checkTrue(is(pedigree(tmp, return.type="data.frame"), "data.frame"))
    tmp2 <- pedigree(tmp, return.type="data.frame")
    ## have to transform sex into a factor...
    PedDf$sex <- factor(PedDf$sex, levels=c("M", "F"))
    rownames(PedDf) <- as.character(PedDf$id)

    checkEquals(tmp2, FamAgg:::sanitizePed(PedDf))

    ## submit a pedigree(List) object...
    Ped <- pedigree(id=PedDf[, 2], famid=PedDf[, 1], dadid=PedDf[, 3],
                    momid=PedDf[, 4], sex=PedDf[, 5])
    pedigree(tmp) <- Ped
    checkTrue(is(pedigree(tmp, return.type="pedigree"), "pedigreeList"))
    ## check the Ped with the pedigree(tmp)
    tmp1 <- do.call(cbind, Ped)
    tmp2 <- do.call(cbind, pedigree(tmp, return.type="pedigree"))
    checkEquals(tmp1, tmp2)
}


test_subset <- function(){
    do.plot <- FALSE
    ## mbsub <- .getMinSub()

    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN

    ## using the constructor
    fad <- FAData(pedigree=mbped)
    ##----------------------------------
    ## subset using boolean
    fadSub <- fad[fad$family==4, ]
    pedSub <- mbped[mbped$family==4, ]
    pedSub$sex <- FamAgg:::sanitizeSex(pedSub$sex)
    pedSub <- FamAgg:::sanitizePed(pedSub)
    checkEquals(pedigree(fadSub), pedSub)
    ## checking kinship
    kinSub <- kinship(pedigree(id=pedSub$id, dadid=pedSub$father, momid=pedSub$mother,
                               famid=pedSub$family, sex=as.numeric(pedSub$sex)))
    checkEquals(kinship(fadSub), kinSub)
    ##-----------------------------------
    ## subset using character
    ids <- c("524", "471", "17", "5", "4", "6")
    pedSub <- mbped[ids, ]
    pedSub$sex <- FamAgg:::sanitizeSex(pedSub$sex)
    fadSub <- fad[ids, ]
    ## Just compare columns other than father and mother, as they have been changed!
    checkEquals(pedigree(fadSub)[, c("family", "id", "sex")], pedSub[, c("family", "id", "sex")])
    ## kinship checking makes no sense... pedigree function does not work!
    ##-----------------------------------
    ## subset using numeric
    idx <- c(4, 2, 6, 112, 12, 14)
    pedSub <- mbped[idx, ]
    pedSub$sex <- FamAgg:::sanitizeSex(pedSub$sex)
    fadSub <- fad[idx, ]
    checkEquals(pedigree(fadSub)[, c("family", "id", "sex")], pedSub[, c("family", "id", "sex")])
    ##-----------------------------------
    ## error testing
    ids <- c("a", "b", "c")
    checkException(fad[ids, ])

    ## Subset to male only:
    fam4 <- fad[fad$family == "4", ]
    fam4M <- fam4[fam4$sex == "M", ]
    checkEquals(all(is.na(fam4M$mother)), TRUE)
    if(do.plot){
        ## That should work.
        plotPed(fam4, family="4")
        ## Now, how is this going to work if we've got females?
        ## plotPed(fam4M, family="4")
    }


    ## manually defining the object.
    fad <- new("FAData")
    pedigree(fad) <- mbped
    fadSub <- fad[fad$family==4, ]
    pedSub <- mbped[mbped$family==4, ]
    pedSub$sex <- FamAgg:::sanitizeSex(pedSub$sex)
    pedSub <- FamAgg:::sanitizePed(pedSub)
    checkEquals(pedigree(fadSub), pedSub)
    ## checking kinship
    kinSub <- kinship(pedigree(id=pedSub$id, dadid=pedSub$father, momid=pedSub$mother,
                               famid=pedSub$family, sex=as.numeric(pedSub$sex)))
    checkEquals(kinship(fadSub), kinSub)

    ## with age.
    ages <- mbsub$endage
    names(ages) <- mbped$id
    fad <- FAData(pedigree=mbped, age=ages)
    ##-----------------------------------
    ## subset using boolean
    fadSub <- fad[fad$family==4, ]
    ageSub <- ages[mbped$family==4]
    checkEquals(age(fadSub), ageSub)
    ##-----------------------------------
    ## subset using character
    ids <- c("524", "471", "17", "5", "4", "6")
    fadSub <- fad[ids, ]
    ageSub <- ages[ids]
    checkEquals(age(fadSub), ageSub)
    ##-----------------------------------
    ## subset using numeric
    idx <- c(4, 2, 6, 112, 12, 14)
    fadSub <- fad[idx, ]
    ageSub <- ages[idx]
    checkEquals(age(fadSub), ageSub)

    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    tcancer <- mbsub$cancer
    names(tcancer) <- mbsub$id
    trait(fad) <- tcancer
    ## subset by boolean.
    pedSub <- PedDf[PedDf$family == 5, ]
    tcancerSub <- tcancer[PedDf$family == 5]
    farSub <- fad[fad$family == 5, ]
    checkEquals(tcancerSub, trait(farSub))
    checkEquals(tcancerSub, farSub$trait)

    ## subset by id.
    ids <- sample(rownames(pedigree(fad)), size=100)
    pedSub <- PedDf[ids, ]
    tcancerSub <- tcancer[ids]
    farSub <- fad[ids, ]
    checkEquals(tcancerSub, trait(farSub))
    checkEquals(tcancerSub, farSub$trait)

    ## subset by idx.
    idx <- sample(1:nrow(pedigree(fad)), size=100)
    pedSub <- PedDf[idx, ]
    tcancerSub <- tcancer[idx]
    farSub <- fad[idx, ]
    checkEquals(tcancerSub, trait(farSub))
    checkEquals(tcancerSub, farSub$trait)

}


## testing the "sex sanitation"
test_sex_columns <- function(){
    sexnum <- c(1, 1, 2, 2, 3, 1, 1, 2, 2, NA)
    expect <- factor(c("M", "M", "F", "F", NA, "M", "M", "F", "F", NA),
                     levels=c("M", "F"))
    got <- FamAgg:::sanitizeSex(sexnum)
    checkEquals(got, expect)
    sexchar <- c("male", "male", "female", "female", "other", "male", "Man", "female", "female", NA)
    got <- FamAgg:::sanitizeSex(sexchar)
    checkEquals(got, expect)
    got <- FamAgg:::sanitizeSex(expect)
    checkEquals(got, expect)
    ## now checking some errors:
    checkException(FamAgg:::sanitizeSex(c(3, 2, 5, 6)))
    checkException(FamAgg:::sanitizeSex(c("woman", "man", "woman", "man")))
}

## test the kinship...
test_kinship <- function(){
    ## mbsub <- .getMinSub()
    ped <- pedigree(id=mbsub$id, dadid=mbsub$fatherid, momid=mbsub$motherid, sex=mbsub$sex,
                    famid=mbsub$famid)
    kin <- kinship(ped)
    fad <- FAData(pedigree=ped)
    return(checkEquals(kin[fad$id, fad$id], kinship(fad)))
}

## testing the constructor.
test_FAData <- function(){
    pedFile <- system.file("txt/minnbreastsub.txt", package="FamAgg")
    checkException(fad <- FAData(pedigree=pedFile, header=TRUE))
    fad <- FAData(pedigree=pedFile, header=TRUE, family.col="famid", id.col="id",
                  father.col="fatherid", mother.col="motherid")
    mbsub <- .getMinSub()
    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    rownames(PedDf) <- as.character(PedDf$id)
    PedDf <- FamAgg:::sanitizePed(PedDf)
    checkEquals(PedDf, pedigree(fad))

    ## adding also age...
    Ages <- as.numeric(mbsub$endage)
    names(Ages) <- mbsub$id
    fad <- FAData(pedigree=pedFile, header=TRUE, family.col="famid", id.col="id",
                  father.col="fatherid", mother.col="motherid", age=Ages)
    checkEquals(Ages, age(fad))
    Ages[1] <- 120
    age(fad) <- Ages
    checkEquals(Ages, age(fad))

    ## construct by submitting a data.frame
    fad <- FAData(PedDf)
    checkEquals(pedigree(fad), PedDf)
}

test_family <- function(){
    ## mbsub <- .getMinSub()
    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    rownames(PedDf) <- as.character(PedDf$id)
    fad <- FAData(PedDf)

    PedDf <- FamAgg:::sanitizePed(PedDf)
    Test <- family(fad, id=7)
    checkEquals(Test, PedDf[PedDf$family==4, ])
}


test_getCommonAncestor <- function(){
    do.plot <- FALSE
    ## mbsub <- .getMinSub()
    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    rownames(PedDf) <- NULL
    fad <- FAData(PedDf)
    ## get common ancestors for some example families
    ancs <- c("759", "760")
    checkEquals(sort(getCommonAncestor(fad, id=c("761", "763", "772"))), sort(ancs))
    if(do.plot)
        plotPed(fad, id=c("761", "763", "772"), device="plot", prune=TRUE)
    ##
    checkException(getCommonAncestor(fad, id="761"))
    checkException(getCommonAncestor(fad, id=c("761", "bla")))
}

test_buildPed <- function(){
    do.plot <- FALSE
    ## mbsub <- .getMinSub()
    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- FamAgg:::.PEDCN
    PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    rownames(PedDf) <- NULL
    fad <- FAData(PedDf)
    if(do.plot){
        switchPlotfun("ks2paint")
        plotPed(fad, id="4", device="plot")
        switchPlotfun()
    }
    expect <- c(5, 26, 11, 12, 13)
    tmp <- buildPed(fad, id=c(26, 13, 12, 11), prune=TRUE)
    checkEquals(sort(tmp$id), sort(expect))

    expect <- c(23, 3, 24, 4, 25, 5, 26, 11, 1, 2)
    tmp <- buildPed(fad, id=c(23, 11), prune=TRUE)
    checkEquals(sort(tmp$id), sort(expect))

    ## build ped for individuals of family 14:
    ## pedigree for individual 447 should include individuals of a second
    ## branch, e.g. individual 26062
    checkTrue(any(buildPed(fad, id="447")$id == "26062"))
    ## for individual 440 (father of 447) this second branch should be
    ## missing
    checkTrue(!any(buildPed(fad, id="440")$id == "26062"))

    checkException(buildPed(fad, family = 2))
    res <- buildPed(fad, family = "5")
    res_2 <- family(fad, family = "5")
    checkTrue(nrow(res) < nrow(res_2))
    res_2 <- removeSingletons(res_2)
    checkEquals(res, res_2)
}


test_dollar <- function(){
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN
    far <- FAData(pedigree=mbped)
    tcancer <- mbsub$cancer
    names(tcancer) <- mbsub$id
    trait(far) <- tcancer

    checkEquals(far$affected, tcancer)
    ## basic stuff from pedigree...
    ids <- as.character(pedigree(far)$id)
    names(ids) <- as.character(ids)
    checkEquals(far$id, ids)
    ## family
    tmp <- pedigree(far)$family
    names(tmp) <- as.character(pedigree(far)$id)
    checkEquals(far$family, tmp)
    ## father
    tmp <- pedigree(far)$father
    names(tmp) <- as.character(pedigree(far)$id)
    checkEquals(far$father, tmp)
    ## mother
    tmp <- pedigree(far)$mother
    names(tmp) <- as.character(pedigree(far)$id)
    checkEquals(far$mother, tmp)

}

## Here we just want to ensure that it's OK to have either 0 or NA for founders in
## the pedigree.
test_kinship2_kinship <- function(){
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN
    pediNumeric <- pedigree(id=mbped$id, famid=mbped$family, dadid=mbped$father,
                            momid=mbped$mother, sex=mbped$sex)
    kinNumeric <- kinship(pediNumeric)
    ## Replacing 0s with NAs:
    mbped$father[mbped$father == 0] <- NA
    mbped$mother[mbped$mother == 0] <- NA
    pediNumericNA <- pedigree(id=mbped$id, famid=mbped$family, dadid=mbped$father,
                              momid=mbped$mother, sex=mbped$sex)
    kinNumericNA <- kinship(pediNumeric)
    checkEquals(kinNumeric, kinNumericNA)
    ## And now character:
    mbped$id <- as.character(mbped$id)
    mbped$father <- as.character(mbped$father)
    mbped$mother <- as.character(mbped$mother)
    pediCharNA <- pedigree(id=mbped$id, famid=mbped$family, dadid=mbped$father,
                           momid=mbped$mother, sex=mbped$sex)
    kinCharNA <- kinship(pediCharNA)
    checkEquals(kinNumeric, kinCharNA)
    ## The same but with "0" for founders.
    mbped$father[is.na(mbped$father)] <- 0
    mbped$mother[is.na(mbped$mother)] <- 0
    checkException(pedigree(id=mbped$id, famid=mbped$family, dadid=mbped$father,
                            momid=mbped$mother, sex=mbped$sex))

    ## "" for "0"
    mbped$father[mbped$father == "0"] <- ""
    mbped$mother[mbped$mother == "0"] <- ""
    pediCharZ <- pedigree(id=mbped$id, famid=mbped$family, dadid=mbped$father,
                          momid=mbped$mother, sex=mbped$sex)
    kinCharZ <- kinship(pediCharZ)
    checkEquals(kinNumeric, kinCharZ)
}

test_export <- function(){
    mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(mbped) <- FamAgg:::.PEDCN
    far <- FAData(pedigree=mbped)
    fam4 <- far[far$family == 4, ]
    tmpF <- paste0(tempfile(), ".ped")
    export(fam4, tmpF, format="ped")
    fam4imported <- FAData(tmpF)
    checkEquals(as.character(pedigree(fam4)$family), pedigree(fam4imported)$family)
    checkEquals(as.character(pedigree(fam4)$id), pedigree(fam4imported)$id)
    checkEquals(as.character(pedigree(fam4)$father), pedigree(fam4imported)$father)
    checkEquals(as.character(pedigree(fam4)$mother), pedigree(fam4imported)$mother)
    checkEquals(as.character(pedigree(fam4)$sex), as.character(pedigree(fam4imported)$sex))
}


