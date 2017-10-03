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

test_dot_binomialTest <- function() {
    ## Working:
    res <- .binomialTest(fad, trait = tcancer)
    checkEquals(res$total_phenotyped[1], sum(!is.na(tcancer)))
    checkEquals(res$total_affected[1], sum(tcancer, na.rm = TRUE))
    checkEquals(as.character(unique(fad$family)), res$family)
    ## Check errors
    checkException(.binomialTest(fad, trait = tcancer, alternative = "other"))
    checkException(.binomialTest(fad))
    checkException(.binomialTest(fad, trait = 1:10))
    ## no affected.
    tr <- rep(0, nrow(fad@pedigree))
    res <- .binomialTest(fad, trait = tr)
    checkEquals(unique(res$pvalue), 1)
    checkEquals(unique(res$affected), 0)
    checkEquals(as.character(unique(fad$family)), res$family)
    ## One affected in one family.
    tr[13] <- 1
    res <- .binomialTest(fad, trait = tr)
    ## NAs
    tr[] <- NA
    res <- .binomialTest(fad, trait = tr)
    checkEquals(as.character(unique(fad$family)), res$family)
    checkTrue(all(is.na(res$pvalue)))
    checkTrue(all(res$phenotyped == 0))
    ## only one family.
    tmp <- fad
    tmp@pedigree$family <- 1
    checkException(.binomialTest(tmp, trait = tcancer))
    suppressWarnings(
        res <- .binomialTest(tmp, trait = tcancer, prob = 1/16)
    )
    checkEquals(nrow(res), 1)
    ## Different alternative.
    res_2 <- .binomialTest(tmp, trait = tcancer, prob = 1/16,
                           alternative = "less", global = TRUE)
    checkTrue(res_2$pvalue != res$pvalue)
    ## Set one family to NA and check that the pvalue is NA for them.
    tr <- tcancer
    tr[fad$family == 5] <- NA
    res <- .binomialTest(fad, trait = tr)
    checkTrue(is.na(res$pvalue[res$family == "5"]))
    checkTrue(all(!is.na(res$pvalue[res$family != "5"])))
    
    ## Set all families except one to NA.
    tr <- tcancer
    tr[] <- NA
    tr[fad$family == 5] <- tcancer[fad$family == 5]
    res <- .binomialTest(fad, trait = tr, global = TRUE, prob = 1/16)
    checkTrue(nrow(res) == 1)

    tr <- tcancer
    tr[] <- NA
    tr[fad$family %in% c(5, 7)] <- tcancer[fad$family %in% c(5, 7)]
    res <- .binomialTest(fad, trait = tr)
    checkTrue(!is.na(res$pvalue[res$family == "5"]))
    checkTrue(!is.na(res$pvalue[res$family == "7"]))
    checkTrue(all(is.na(res$pvalue[!(res$family %in% c("5", "7"))])))
}

test_validateFABinTestResults <- function() {
    tmp <- new("FABinTestResults")
    checkTrue(validObject(tmp))
    checkTrue(FamAgg:::.validateFABinTestResults(tmp))
    ## Modify the object
    tmp@result <- tmp@result[, 1:3, drop = FALSE]
    checkException(validObject(tmp))
    checkTrue(is(FamAgg:::.validateFABinTestResults(tmp), "character"))
}

test_binomialTest <- function() {
    res <- binomialTest(fad, trait = tcancer, traitName = "Breast Cancer")
    checkTrue(validObject(res))
    checkEquals(as.character(unique(fad$family)), res@result$family)
    checkTrue(all(res@result$pvalue < 1))
    res_tab <- result(res)
    checkTrue(all(res_tab$trait_name == "Breast Cancer"))
    checkEquals(res_tab$family, res_tab[order(res_tab$pvalue), "family"])
    checkEquals(res_tab$pvalue, sort(res@result$pvalue))
    trait(res) <- tcancer
    checkEquals(nrow(res@result), 0)
    ## Can run the analysis on this?
    res <- binomialTest(res, trait = tcancer)
    res_2 <- binomialTest(fad, trait = tcancer)
    checkEquals(res, res_2)

    res_tab <- result(res, method = "BH")
    res_tab_2 <- result(res, method = "bonferroni")
    checkTrue(all(res_tab$padj < res_tab_2$padj))
    
    ## all NA trait.
    tr <- rep(NA, length(fad$id))
    names(tr) <- fad$id
    checkException(binomialTest(fad, trait = tr))

    ## All except one family NA
    tr[fad$family == 8] <- tcancer[fad$family == 8]
    checkException(binomialTest(fad, trait = tr))
    res <- binomialTest(fad, trait = tr, global = TRUE, prob = 1/16)
    checkEquals(nrow(res@result), 1)
    checkTrue(res@result$pvalue[1] < 0.05)
    res_tab <- result(res)
    checkTrue(nrow(res_tab) == 1)
    
    res_2 <- binomialTest(fad, trait = tcancer, global = TRUE, prob = 1/16)
    checkTrue(res_2@result$pvalue > res@result$pvalue)
    
    ## Errors
    checkException(binomialTest(fad))
    checkException(binomialTest(fad, trait = 1:4))
    checkException(binomialTest(fad, trait = rep(4, length(fad$id))))
}
