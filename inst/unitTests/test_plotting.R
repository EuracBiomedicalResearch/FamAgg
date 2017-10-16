data(minnbreast)
mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
rownames(mbsub) <- mbsub$id
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- FamAgg:::.PEDCN
PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
fad <- FAData(PedDf)
## trait...
far <- fad
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id
trait(far) <- tcancer

## OK, 0.0.6
## test the argument checking function of the plot functions.
test_plot_arg_numeric_logical <- function(){
    ## test numeric input allowing only 0 and 1
    testvec <- c("0", "1", "0", "0", "0", "1")
    checkEquals(as.numeric(testvec), FamAgg:::checkArgNumericLogical(testvec, 6))
    testvec <- c(testvec, "z")
    checkException(FamAgg:::checkArgNumericLogical(testvec, 7))
    testvec <- rep(TRUE, 6)
    checkEquals(as.numeric(testvec), FamAgg:::checkArgNumericLogical(testvec, 6))
    ## check if we get a vector of 6 0s
    checkEquals(rep(0, 6), FamAgg:::checkArgNumericLogical(ninds=6))
}

test_plot_arg_are_twins <- function(){
    testvec <- c(NA, NA, "m_1", NA, "m_1", NA)
    testveccomp <- testvec
    testveccomp[is.na(testveccomp)] <- ""
    checkEquals(FamAgg:::checkArgAreTwins(testvec, 6), testveccomp)
    ## empty
    checkEquals(FamAgg:::checkArgAreTwins(ninds=6), rep("", 6))
}

test_plot_arg_consanguineous <- function(){
    testvec <- c(NA, NA, NA, "a", NA, "a")
    testveccomp <- testvec
    testveccomp[is.na(testveccomp)] <- ""
    checkEquals(FamAgg:::checkArgConsanguineous(testvec, 6), testveccomp)
    ## check error:
    testvec <- c(testvec, "b")
    checkException(FamAgg:::checkArgConsanguineous(testvec, 7))
    ## empty.
    checkEquals(FamAgg:::checkArgConsanguineous(ninds=6), rep("", 6))
}

test_kinship2_plot <- function(){
    ## data(minnbreast)
    ## mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
    ## rownames(mbsub) <- mbsub$id
    ## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(PedDf) <- FamAgg:::.PEDCN
    ## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    ## far <- as(FAData(pedigree=PedDf), "FAResult")
    ## ## trait...
    ## tcancer <- mbsub$cancer
    ## names(tcancer) <- mbsub$id
    ## trait(far) <- tcancer

    ped <- family(far, family="4", return.type="pedigree")
    ped <- ped["4"]
    plot(ped)

    ## directly plotting...
    fam <- family(far, family="4", return.type="data.frame")
    FamAgg:::ks2paint(family=fam$family, individual=fam$id, father=fam$father,
                      mother=fam$mother, gender=fam$sex, affected=fam$affected,
                      dev="plot")
    ## random deads
    is.deceased <- sample(c(TRUE, FALSE), nrow(fam), replace=TRUE)
    FamAgg:::ks2paint(family=fam$family, individual=fam$id, father=fam$father,
                      mother=fam$mother, gender=fam$sex, affected=fam$affected,
                      is.deceased=is.deceased, dev="plot")
}

test_haplopaint_txt_plot <- function() {
    fam <- family(far, family = "4", return.type = "data.frame")
    is.deceased <- sample(c(TRUE, FALSE), nrow(fam), replace=TRUE)
    suppressWarnings(
        tmp <- FamAgg:::haplopaint(family = fam$family, individual = fam$id,
                                   father = fam$father, mother = fam$mother,
                                   gender = fam$sex, affected = fam$affected,
                                   is.deceased = is.deceased, dev = "txt")
    )
    checkTrue(is(tmp, "character"))
    res <- read.table(tmp, sep = "\t", as.is = TRUE, header = TRUE,
                      comment.char = "")
    checkEquals(as.integer(is.deceased)[fam$id %in% res$INDIVID], res$DEAD)
    pf <- options()$FamAgg$plotfun
    FA <- options()$FamAgg
    FA$plotfun <- "haplopaint"
    options(FamAgg = FA)
    suppressWarnings(
        tmp <- plotPed(far, family = "5", device = "txt")
    )
    checkTrue(is(tmp, "character"))
    res <- read.table(tmp, sep = "\t", as.is = TRUE, header = TRUE,
                      comment.char = "")
    checkTrue(ncol(res) == 18)
    FA$plotfun <- pf
    options(FamAgg = FA)
}

notrun_test_plot <- function(){
    ## ## again, make a simple pedigree with trait.
    ## data(minnbreast)
    ## mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
    ## rownames(mbsub) <- mbsub$id
    ## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(PedDf) <- FamAgg:::.PEDCN
    ## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    ## far <- as(FAData(pedigree=PedDf), "FAResult")
    ## ## trait...
    ## tcancer <- mbsub$cancer
    ## names(tcancer) <- mbsub$id
    ## trait(far) <- tcancer

    checkException(FamAgg:::buildHaplopaintDataframe(family=fam$family,
                                                     individual=fam$id,
                                                     father=fam$father,
                                                     mother=fam$mother,
                                                     gender=fam$sex,
                                                     affected=fam$affected,
                                                     is.deceased=c(FALSE, TRUE)))

    ## plain function to generate...
    fam <- family(far, family="4")
    rownames(fam) <- as.character(fam$id)
    ## removing childless founders?
    ## df <- buildHaplopaintDataframe(family=fam$family, individual=fam$id, father=fam$father,
    ##                                mother=fam$mother, gender=fam$sex, affected=fam$affected)
    ## kins <- kinship(id=fam$id, momid=fam$mother, dadid=fam$father, sex=as.numeric(fam$sex))
    ## cols <- colSums(kins)
    ## if(any(cols <= 0.5)){
    ##     kins <- kins[cols > 0.5, cols > 0.5]
    ##     fam <- fam[colnames(kins), ]
    ##     warning(paste0("Removed ", sum(cols <= 0.5), " childless founders!"))
    ## }
    ## plain one...
    res <- FamAgg:::haplopaint(family=fam$family, individual=fam$id, father=fam$father,
                               mother=fam$mother, gender=fam$sex, affected=fam$affected)
    ## add more stuff... random deads...
    deads <- sample(c(TRUE, FALSE), size=nrow(fam), replace=TRUE)
    res <- FamAgg:::haplopaint(family=fam$family, individual=fam$id, father=fam$father,
                               mother=fam$mother, gender=fam$sex, affected=fam$affected,
                               is.deceased=deads)
}


## testing some larger pedigrees... whether that works at all...
notrun_test_larger_ped <- function(){
    ## data(minnbreast)
    ## mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
    ## PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    ## colnames(PedDf) <- FamAgg:::.PEDCN
    ## PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
    ## far <- as(FAData(pedigree=PedDf), "FAResult")
    ## ## trait...
    ## tcancer <- mbsub$cancer
    ## names(tcancer) <- mbsub$id

    tcancershuffle <- sample(tcancer, length(tcancer))
    trait(far) <- tcancershuffle

}



