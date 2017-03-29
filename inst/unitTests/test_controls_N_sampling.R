## if this setup of a global data set does not work we might want to put that into
## the functions .setUp() and .tearDown(), which are executed before and after each
## test function.
data(minnbreast)
ped <- minnbreast[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(ped) <- c("family", "id", "father", "mother", "sex")
fam4 <- ped[ped$family == "4", ]
fad <- FAData(ped[ped$family %in% c(4:6), ])

## testing functions and methods related to matched controls and sampling.
test_matched_generation <- function(){
    ## data(minnbreast)
    ## ped <- minnbreast[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(ped) <- c("family", "id", "father", "mother", "sex")
    ## fam4 <- ped[ped$family == "4", ]
    ## fad <- FAData(ped[ped$family %in% c(4:20), ])

    ##
    genMatch <- FamAgg:::getGenerationMatched(pedigree(fad),
                                              id=c("4", "5", "56", "76", "19"))
    ## we get a list with the list names being the family id, and elements the
    ## ids of the individuals from the same generation.
    genMatch <- unlist(FamAgg:::getGenerationMatched(pedigree(fad),
                                                     id="4"), use.names=FALSE)
    checkEquals(sort(genMatch), sort(as.character(c(25, 4, 5, 26, 6, 27, 7,
                                                    28, 8, 29, 9, 10))))
    genMatch <- unlist(FamAgg:::getGenerationMatched(pedigree(fad), id=c("4", "2")),
                       use.names=FALSE)
    checkEquals(sort(genMatch), sort(as.character(c(25, 4, 5, 26, 6, 27, 7, 28,
                                                    8, 29, 9, 10, 1, 2))))
}

test_matched_sex_generation <- function(){
    genSexMatch <- unlist(FamAgg:::getGenerationSexMatched(pedigree(fad), id="4"),
                          use.names=FALSE)
    checkEquals(sort(genSexMatch), sort(as.character(c(4, 26, 27, 8, 9, 10))))
    ## when i use 2 it gets useless if the given ids have different sexes
    genSexMatch <- unlist(FamAgg:::getGenerationSexMatched(pedigree(fad),
                                                           id=c("4", "1")), use.names=FALSE)
    checkEquals(sort(genSexMatch), sort(as.character(c(25, 4, 5, 26, 6, 27, 7, 28, 8,
                                                       29, 9, 10, 1, 2))))
    ## it's different if both are from the same sex:
    genSexMatch <- unlist(FamAgg:::getGenerationSexMatched(pedigree(fad), id=c("25", "21")),
                          use.names=FALSE)
    checkEquals(sort(genSexMatch), sort(as.character(c(25, 5, 6, 7, 28, 29, 21))))
}

test_get_all <- function(){
    all <- unlist(FamAgg:::getAll(pedigree(fad), id="4"), use.names=FALSE)
    checkEquals(sort(all), sort(as.character(family(fad, family="4")[, "id"])))
}

test_sex_matched <- function(){
    females <- unlist(FamAgg:::getSexMatched(pedigree(fad), id="4"), use.names=FALSE)
    fam4 <- family(fad, "4")
    checkEquals(sort(females), sort(as.character(fam4[fam4$sex == "F", "id"])))
    males <- unlist(FamAgg:::getSexMatched(pedigree(fad), id="1"), use.names=FALSE)
    checkEquals(sort(males), sort(as.character(fam4[fam4$sex == "M", "id"])))
}

test_external_matched <- function(){
    ## use sex for matching; can compare that then to sex matching.
    extmatch <- unlist(FamAgg:::getExternalMatched(pedigree(fad), id="4",
                                                   match.using=fad$sex), use.names=FALSE)
    sexmatch <- unlist(FamAgg:::getSexMatched(pedigree(fad), id="4"),
                       use.names=FALSE)
    checkEquals(sort(extmatch), sort(sexmatch))
}


notrun_test_matched_generation_n_frequency <- function(){
    ## here we want to test getting generation matched controls for the ids and
    ## considering also the frequency of generation of the ids.
    ## data(minnbreast)
    ## ped <- minnbreast[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(ped) <- c("family", "id", "father", "mother", "sex")
    ## fam4 <- ped[ped$family == "4", ]
    ## fad <- FAData(ped[ped$family %in% c(4:20), ])

    genFreq <- FamAgg:::getGenerationMatchedFrequency(fam4, id=c("1", "5", "6", "22"))
    genFreq <- unlist(genFreq, use.names=FALSE)
    genFreq <- table(genFreq)
    ## the control ids:
    names(genFreq)
    ## the frequency
    as.numeric(genFreq)
}

notrun_test_sex_match <- function(){
    affIds <- c("4", "11")
    idL <- FamAgg:::getGenerationMatchedByIndividual(pedigree(fad), id=affIds)
    SexMatched <- FamAgg:::doSexMatch(pedigree(fad), idList=idL)
    haveS <- fad$sex[names(SexMatched)[1]]
    needS <- fad$sex[SexMatched[[1]]]
    checkEquals(all(as.character(needS)==as.character(haveS)), TRUE)
    haveS <- fad$sex[names(SexMatched)[2]]
    needS <- fad$sex[SexMatched[[2]]]
    checkEquals(all(as.character(needS)==as.character(haveS)), TRUE)

    ## sexmatched
    matches <- FamAgg:::getSexMatchedByIndividual(ped=family(fad, "4"), id=affIds)
    haveS <- fad$sex[names(matches)[1]]
    needS <- fad$sex[matches[[1]]]
    checkEquals(all(as.character(needS)==as.character(haveS)), TRUE)
    ## next
    haveS <- fad$sex[names(matches)[2]]
    needS <- fad$sex[matches[[2]]]
    checkEquals(all(as.character(needS)==as.character(haveS)), TRUE)

}

notrun_test_stratsample <- function(){
    ## hm, check stratsample...
    Test <- family(fad, "4")
    ## require(sampling)
    ## strata(data=Test, stratanames="sex", size=c(M=1, F=4))
    ## Test <- Test[Test$sex=="M", ]
    ## strata(data=Test, stratanames="sex", size=c(M=2, F=3), method="srswor")

    ## Test <- family(fad, "4")
    ## Test <- Test[Test$sex=="F", ]
    ## strata(data=Test, stratanames="sex", size=c(M=2, F=3), method="srswor")

    ## Test <- family(fad, "4")
    ## strata(data=Test, stratanames="sex", size=c(M=2, F=3), method="srswor")
    ## Test <- Test[order(as.character(Test$sex)), ]
    ## strata(data=Test, stratanames="sex", size=c(M=2, F=3), method="srswor")
    ## ## Oh man. that's not good...
    ## ## what if we had NAs???
    ## ## have to exclude those guys... dammit
    ## Test[1:3, "sex"] <- NA
    ## strata(data=Test, stratanames="sex", size=c(M=2, F=3), method="srswor")
    ## what with str
    library(survey)
    idx <- stratsample(Test$sex, c(M=4, F=3))
    Test$sex[idx]
    Test$sex[stratsample(Test$sex, c(M=4, F=3))]

    Test$sex[stratsample(Test$sex, c(F=3, M=6))]

    ## make 1 male 2 females:
    affIds <- c("5", "8", "20")
    counts <- table(Test[affIds, "sex"])
    idx <- stratsample(Test$sex, counts)
    Test[idx, "sex"]

    ## only females
    affIds <- c("4", "8", "20")
    counts <- table(Test[affIds, "sex"])
    counts <- counts[counts > 0]
    idx <- stratsample(Test$sex, counts)
    Test[idx, "sex"]

    ## do it on the full table
    ##ped <- ped[!is.na(ped$sex), ]
    ##system.time(ped$sex[stratsample(ped$sex, c(M=4, F=3))])
    ##system.time(strata(ped, stratanames="sex", size=c(M=4, F=3), method="srswor"))
    ## stratsample is way faster.
}



