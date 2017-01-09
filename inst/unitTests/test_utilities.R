pedFile <- system.file("txt/Large.ped", package="FamAgg")
largePed <- read.table(pedFile, sep=";")[, 1:5]
colnames(largePed) <- c("family", "id", "father", "mother", "sex")
rownames(largePed) <- as.character(largePed$id)
largeFad <- FAData(pedigree=largePed)

data(minnbreast)
## head(minnbreast)
minPed <- minnbreast[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(minPed) <- c("family", "id", "father", "mother", "sex")
minFad <- FAData(pedigree=minPed[minPed$family %in% 4:10, ])

test_ped2graph <- function(){
    allped <- largePed

    pedids <- c("3511", "3225", "3508", "3550", "3507", "3510", "3224",
                "3052", "3250", "2969", "3407")
    ped <- largePed[pedids, ]
    ped[!(ped$father %in% ped$id), "father"] <- 0
    ped[!(ped$mother %in% ped$id), "mother"] <- 0
    plot(pedigree(id=ped$id, dadid=ped$father, momid=ped$mother, sex=ped$sex))
    igr <- ped2graph(ped)
    igr <- add_vertices(igr, 3)
    plot(igr)


    ## distance matrix: this lists us all that we can get to following the direction
    dout <- distances(igr, mode="out")
    ## got the way back: from children -> parents
    din <- distances(igr, mode="in")
    ## use the ROWS in that.
    testids <- c("3511", "3507", "3224")
    subD <- din[testids, , drop=FALSE]
    subD <- subD[, !apply(subD, MARGIN=2, function(x){
        return(all(FamAgg:::infOrZero(x)))
    }) , drop=FALSE]

    ## hugegraph
    IGR <- ped2graph(allped)
    Dist <- distances(IGR)
    any(is.infinite(Dist[pedids, pedids])) ## means, they are connected!
    testids <- c("3224", "3511", "255")
    any(is.infinite(Dist[testids, testids]))
    any(is.infinite(Dist))
    ## apparently, that thing is totally connected...
    return(TRUE)
}

test_subPedigree <- function(){
    ped <- largePed
    testids <- c("3511", "3225", "3508", "3550", "3507", "3510", "3224")
    ## get the subPedigree for these
    subped <- subPedigree(ped, id=testids)

    checkEquals(sum(testids %in% as.character(subped$id)), length(testids))

    return(TRUE)
}


test_getCommonAncestor <- function(){
    do.plot <- FALSE
    ped <- largePed
    fad <- largeFad
    ##testids <- c("3511", "3225", "3508")  ## no common graph
    testids <- c("3511", "3225", "3508", "3550", "3507", "3510", "3224")
    ancs <- c("3052", "2969")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))

    ## plotting the pedigree for this...
    if(do.plot){
        switchPlotfun("ks2paint")
        plotPed(fad, id=testids, device="plot", prune=TRUE)
    }

    ## 1 generation up:
    testids <- c("3508", "3507")
    ancs <- c("3250", "3225")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    if(do.plot)
        plotPed(fad, id=testids, device="plot", prune=TRUE)
    ## OK

    ## 2 generations up:
    testids <- c("3550", "3511")
    ancs <- c("3052", "2969")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    if(do.plot)
        plotPed(fad, id=testids, device="plot", prune=TRUE)
    ## OK

    ## x generations?
    testids <- c("3224", "3250")
    getCommonAncestor(ped, id=testids)
    if(do.plot)
        plotPed(fad, id=c(testids, "2042"), device="plot", prune=TRUE)
    ## lovely pedigree ;)
    ## switchPlotfun()
    ## plotPed(fad, id=c(testids, "2042"), filename="~/Desktop/test.pdf", prune=TRUE)

    ## ids from different generations:
    testids <- c(3508, 3224, 3225)
    ancs <- c("3052", "2969")
    getCommonAncestor(ped, id=testids)
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    ## OK

    ## doesn't work if one is already the ancestor...
    testids <- c(3508, 3052)
    ancs <- c("3052", "2969")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    if(do.plot)
        plotPed(fad, id=c(testids, ancs), device="plot", prune=TRUE)

    ## test an example that can not work.
    testids <- c(3016, 3349)
    ancs <- getCommonAncestor(ped, testids)
    ## HECK?
    if(do.plot){
        plotPed(fad, id=c(testids), device="plot", prune=TRUE)
        plotPed(fad, id=c(testids, ancs), device="plot", prune=TRUE)
    }

    ## ????
    testids <- c(3016, 2962)
    ancs <- getCommonAncestor(ped, testids)
    if(do.plot){
        plotPed(fad, id=c(testids), device="plot", prune=TRUE)
        plotPed(fad, id=c(testids, ancs), device="plot", prune=TRUE)
    }

    ##************************************************
    ##*
    ##* test it with the minnesota breast cancer set:
    ##*
    ##table(fad$family)
    fad <- minFad
    ped <- pedigree(fad)

    if(do.plot)
        plotPed(fad, family="6", device="plot")
    testids <- c("103", "104", "95", "99")
    ancs <- c("80", "81")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    if(do.plot){
        plotPed(fad, id=c(testids, ancs), device="plot")
        plotPed(fad, id=c(testids, ancs), device="plot", prune=TRUE)
    }

    ## other example:
    testids <- c(101, 83)
    ancs <- c("107", "84")
    checkEquals(sort(getCommonAncestor(ped, id=testids)), sort(ancs))
    if(do.plot)
        plotPed(fad, id=c(testids), device="plot", prune=TRUE)

    ## can not find one.
    testids <- c("88", "93", "114")
    checkEquals(getCommonAncestor(ped, id=testids), NA)
    return(TRUE)
}

test_count_generations <- function(){
    ped <- minPed
    expect <- 3
    names(expect) <- "2"
    checkEquals(countFenerations(ped=ped, id=2), expect)

    expect <- 2
    names(expect) <- "11"
    checkEquals(countGenerations(ped=ped, id=11, direction="up"), expect)
}

test_find_founders <- function(){
    fad <- minFad
    ped <- pedigree(fad)
    family <- "4"
    checkEquals(findFounders(ped, family=family), c("1", "2"))
    ## for FAData
    fad <- FAData(ped[ped$family=="4", ])
    checkEquals(findFounders(fad, family=family), c("1", "2"))
    checkEquals(findFounders(pedigree(fad, return.type="pedigree"),
                             family=family), c("1", "2"))
}

test_find_foundersForId <- function() {
    ## Check if we get the same results for findFounders for family and for
    ## findFounders with an id within the family.
    fad <- minFad
    ped <- pedigree(fad)
    fnds <- findFounders(ped, id = 269)
    fnds_2 <- findFounders(ped, family = 10)
    checkEquals(fnds, fnds_2)
    fnds <- findFounders(ped, id = 20)
    fnds_2 <- findFounders(ped, family = 4)
    checkEquals(fnds, fnds_2)
    ## The same for FAData
    fnds <- findFounders(fad, id = 269)
    fnds_2 <- findFounders(fad, family = 10)
    checkEquals(fnds, fnds_2)
    fnds <- findFounders(fad, id = 20)
    fnds_2 <- findFounders(fad, family = 4)
    checkEquals(fnds, fnds_2)
}

test_find_siblings <- function(){
    ped <- minPed
    checkEquals(getSiblings(ped, id="11"), c("11", "12", "13"))
    ## for FAData
    fad <- FAData(ped[ped$family=="4", ])
    checkEquals(getSiblings(fad, id="11"), c("11", "12", "13"))
    checkEquals(getSiblings(pedigree(fad, return.type="pedigree"),
                            id="11"), c("11", "12", "13"))
}

test_do_get_ancestors <- function(){
    ped <- minPed
    checkEquals(getAncestors(ped, id="11"), c("5", "26", "1", "2"))
    ## for FAData
    fad <- FAData(ped[ped$family=="4", ])
    checkEquals(getAncestors(fad, id="11"), c("5", "26", "1", "2"))
    checkEquals(getAncestors(pedigree(fad, return.type="pedigree"),
                             id="11"), c("5", "26", "1", "2"))
    expect <- as.character(c(3, 4, 24, 25))
    checkEquals(sort(getAncestors(fad, c(22, 3), max.generations=1)), sort(expect))
}

test_do_get_children <- function(){
    ped <- minPed
    checkEquals(getChildren(ped, id="3"), c("21", "22", "23"))
    ## for FAData
    fad <- FAData(ped[ped$family=="4", ])
    checkEquals(getChildren(fad, id="3"), c("21", "22", "23"))
    checkEquals(getChildren(pedigree(fad, return.type="pedigree"), id="3"), c("21", "22", "23"))
    expect <- as.character(c(4, 5, 6, 7, 8, 9, 10, 17))
    checkEquals(sort(getChildren(fad, id=c(1, 28, 2, 8), max.generation=1)), sort(expect))
}

test_count_generations <- function(){
    ped <- minPed
    expect <- 1
    names(expect) <- "28"
    checkEquals(countGenerations(ped, id="28"), expect)
    ## for FAData
    fad <- FAData(ped[ped$family=="4", ])
    checkEquals(countGenerations(fad, id="28"), expect)
    checkEquals(countGenerations(pedigree(fad, return.type="pedigree"),
                                 id="28"), expect)
}

test_estimate_generations <- function(){
    ped <- minPed
    fam4 <- ped[ped$family == "4", ]
    fad <- FAData(ped[ped$family %in% c(4:20), ])

    ## testing stuff...
    from8 <- c(-1, -1, 1, 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,1,0,0,0,0,0)
    names(from8) <- as.character(1:length(from8))
    ## gens <- FamAgg:::doGetGenerationFrom(fam4, id=8)
    ## checkEquals(from8, gens[names(from8)])
    gens2 <- generationsFrom(fam4, id=8)
    checkEquals(from8, gens2[names(from8)])

    ## artificially create a second branch for fam4:
    addon <- data.frame(family=4,
                        id=50:59,
                        father=c(0, 0, 50, 50, 0, 54, 0, 54, 57, 57),
                        mother=c(0, 0, 51, 51, 0, 53, 0, 53, 56 ,56),
                        sex=c("M", "F", "M", "F", "M", "F", "F", "M", "F", "M")
                        )
    fam4mod <- rbind(fam4, addon)
    fam4mod[fam4mod$id == "25", c("father", "mother")] <- c(54, 53)
    checkEquals(findFounders(fam4mod, family="4"), c("50", "51"))

    addon2 <- data.frame(family=4,
                         id=70:76,
                         father=c(0, 70, 70, 70, 0, 71, 71),
                         mother=c(0, 20 ,20, 20, 0, 74, 74),
                         sex=c("M", "M", "F", "M", "F", "M", "F"))
    fam4mod2 <- rbind(fam4mod, addon2)
    checkEquals(findFounders(fam4mod2, 4), c("1", "2"))

    ## these modified families are interesting to test...
    fam4from6 <- c(-1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0)
    names(fam4from6) <- as.character(1:29)
    checkEquals(fam4from6, generationsFrom(fam4, id="6")[names(fam4from6)])

    addonfrom6 <- c(-2, -2, -1, -1, -1, 0, 0, 0, 1, 1)
    names(addonfrom6) <- as.character(50:59)
    fam4modfrom6 <- c(fam4from6, addonfrom6)
    checkEquals(fam4modfrom6, generationsFrom(fam4mod, id="6")[names(fam4modfrom6)])

    addon2from6 <- c(1, 2, 2, 2, 2, 3, 3)
    names(addon2from6) <- as.character(70:76)
    fam4mod2from6 <- c(fam4modfrom6, addon2from6)
    checkEquals(fam4mod2from6, generationsFrom(fam4mod2, id="6")[names(fam4mod2from6)])

    ## estimating generations...
    Gens <- estimateGenerations(fam4)[[1]]
    checkEquals(Gens[names(fam4from6)], fam4from6+1)
    ## what for fam4mod
    findFounders(fam4mod, 4)
    Gens <- estimateGenerations(fam4mod)[[1]]
    checkEquals(Gens[names(fam4modfrom6)], fam4modfrom6+2)
    ## what for fam4mod2
    findFounders(fam4mod2, 4)
    Gens <- estimateGenerations(fam4mod2)[[1]]
    checkEquals(Gens[names(fam4mod2from6)], fam4mod2from6+1)
}

## getFounders, getSingletons
test_getFounders_getSingletons <- function(){
    ped <- minPed
    fad <- FAData(ped)
    checkEquals(getFounders(fad), getFounders(ped))
    checkEquals(getSingletons(fad), getSingletons(ped))
    ## Now check if they are really childless...
    FounderPed <- ped[ped$id %in% getFounders(ped), ]
    checkEquals(sum(FounderPed[, c("father", "mother")]), 0)
    ## Check if childless founders really have no children.
    clf <- getSingletons(fad)
    checkEquals(any(ped$father %in% clf), FALSE)
    checkEquals(any(ped$mother %in% clf), FALSE)
}

## Test the difference between removeSingletons and doPrunePed
test_singletons_vs_prune <- function(){
    ped <- minPed
    ## using prune
    system.time(
        pedPruned <- FamAgg:::doPrunePed(ped)
    )
    system.time(
        pedRemSin <- removeSingletons(ped)
    )
    pedPruned <- FamAgg:::sanitizePed(pedPruned)
    rownames(pedRemSin) <- pedRemSin$id
    rownames(pedPruned) <- pedPruned$id
    checkEquals(pedPruned, pedRemSin)
    ## OK, it's the same; now let's check what happens if we want to add missing
    ## mates on them.
    system.time(
        pedPruned <- FamAgg:::doPrunePed(ped, addMissingMates=TRUE)
    )
    checkEquals(nrow(pedPruned), nrow(pedRemSin))
    ## Eventually I might want to use the removeSingletons instead of the prune.
}

test_removeSingletons <- function() {
    ped <- minFad
    sings <- getSingletons(ped)
    pedSub <- removeSingletons(ped)
    sings_2 <- getSingletons(pedSub)
    checkTrue(length(sings_2) == 0)
    checkTrue((length(sings) + nrow(pedigree(pedSub))) == nrow(pedigree(ped)))
    trait <- sample(c(0, 1), nrow(pedigree(ped)), replace = TRUE)
    names(trait) <- ped$id
    trait(ped) <- trait
    ## Do the subsetting with the trait:
    pedSub <- removeSingletons(ped)
    checkEquals(trait(pedSub), trait[pedSub$id])
}




