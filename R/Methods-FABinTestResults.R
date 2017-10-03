## Results for the binomial test.

.validateFABinTestResults <- function(object) {
    msg <- character()
    ## Check that we've got all required columns
    cn <- colnames(object@result)
    exp_cn <- c("total_phenotyped", "total_affected", "family", "phenotyped",
                "affected", "pvalue", "prob")
    if (!all(exp_cn %in% cn))
        msg <- c(msg,
                 paste0("Expected columns: ", paste(exp_cn, collapse = ", "),
                 " but got: ", paste(cn, collapse = ", ")))
    if (length(msg))
        msg
    else
        TRUE
}
setValidity("FABinTestResults", .validateFABinTestResults)

setMethod("show", "FABinTestResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    cat(paste0(" * Number of rows of the result data.frame: ",
               nrow(object@result),".\n"))
})

#' @description Performs a binomial test on the provided data.
#'
#' @details `NA`s are removed prior to testing. If the pedigree does not provide
#'     a family id or if the whole pedigree consists of a single family, by
#'     default a *global* test is conducted. This requires parameter `prob` to
#'     be defined (i.e. the probability of being affected from the whole
#'     population).
#'
#' @param x `FAData` or `data.frame` with a column `"family"`.
#'
#' @param trait `logical` or `integer` defining the trait. We assume that these
#'     are in the correct format, i.e. 0,1,NA or TRUE,FALSE,NA.
#'
#' @param global `logical` defining whether a global or a per-family test should
#'     be conducted. If no families are available in 
#'
#' @param prob `numeric(1)` providing the probability of being affected in the
#'     whole population.
#'
#' @param alternative `character(1)` indicates the alternative hypothesis and
#'     must be one of `"greater"` (default), `"less"` or `"two.sided"`.
#' 
#' @return `data.frame` with the results of the binomial test. One row for each
#' family (or a single row if `global = TRUE`). Columns are:
#' - `"total_phenotyped"`
#' - `"total_affected"`
#' - `"family"`: the ID (number) of the family.
#' - `"phenotyped"`
#' - `"affected"`
#' - `"pvalue"`: the p-value.
#' - `"prob"`: the probability of being affected. Either a *local* probability
#'    calculated based on all affected and phenotyped in the whole pedigree or
#'    a *global* (population) probability that was provided with argument
#'    `prob`.
#' 
#' @md
#'
#' @author Johannes Rainer
#'
#' @noRd
.binomialTest <- function(x, trait = NULL, global = is.null(x$family),
                          prob = NULL, alternative = c("greater", "less",
                                                       "two.sided")) {
    alternative <- match.arg(alternative)
    if (is.null(x$family) | global)
        fams <- rep_len("full pedigree", nrow(x@pedigree))
    else
        fams <- x$family
    ## Check that length of trait matches length of fams.
    if (length(fams) != length(trait))
        stop("length of 'trait' has to match the number of individuals in ",
             "the pedigree")
    fam_ids <- as.character(unique(fams))
    is_na <- is.na(trait)
    ## Kick out NAs.
    trait <- trait[!is_na]
    if (length(trait)) {
        ## Only perform the tests if we have something to test!
        fams <- fams[!is_na]
        trait_by_fam <- split(trait, f = fams)
        ## Now, if we have a single family assume we do a global test.
        if (length(trait_by_fam) == 1 & !global) {
            warning("Have a single family with phenotyped individuals. Will ",
                    "switch to a global test.")
            global <- TRUE
        }
        ## global = TRUE requires definition of prob.
        if (global & is.null(prob))
            stop("Argument 'prob' has to be specified if a global test is ",
                 "performed")
        total_phenotyped <- length(trait)
        total_affected <- sum(trait)
        fam_phenotyped <- lengths(trait_by_fam)
        fam_affected <- unlist(lapply(trait_by_fam, sum))
        if (is.null(prob))
            prob <- total_affected / total_phenotyped
        pvals <- .binom.test(fam_affected, fam_phenotyped, prob = prob,
                             alternative = alternative)
        ## Compile results data.frame
        res <- data.frame(total_phenotyped = total_phenotyped,
                          total_affected = total_affected,
                          family = names(trait_by_fam),
                          phenotyped = fam_phenotyped,
                          affected = fam_affected,
                          pvalue = pvals,
                          prob = prob,
                          check.names = FALSE, stringsAsFactors = FALSE)
    } else {
        if (is.null(prob))
            prob <- NA
        trait_by_fam <- c()
        res <- data.frame(total_phenotyped = integer(),
                          total_affected = integer(),
                          family = character(),
                          phenotyped = integer(),
                          affected = integer(),
                          pvalue = numeric(),
                          prob = numeric(),
                          check.names = FALSE, stringsAsFactors = FALSE)
    }
    ## Add also results for families without affected.
    if (any(!(fam_ids %in% names(trait_by_fam)))) {
        na_fams <- fam_ids[!(fam_ids %in% names(trait_by_fam))]
        na_data <- data.frame(total_phenotyped = total_phenotyped,
                              total_affected = total_affected,
                              family = na_fams,
                              phenotyped = 0,
                              affected = 0,
                              pvalue = NA,
                              prob = prob,
                              check.names = FALSE, stringsAsFactors = FALSE)
        res <- rbind(res, na_data)
        
    }
    ## Re-order them.
    rownames(res) <- as.character(res$family)
    res[fam_ids, ]
}

## binomial test. Since that's a function it can be run on any objects
## inheriting FAData.
binomialTest <- function(object, trait, traitName, global = FALSE,
                         prob = NULL, alternative = c("greater", "less",
                                                      "two.sided")) {
    alternative <- match.arg(alternative)
    if (!is(object, "FAData"))
        stop("'object' is expected to be a 'FAData' object but I got '",
             class(object), "'")
    if (!is(object, "FABinTestResults"))
        object <- as(object, "FABinTestResults")
    if (missing(trait))
        trait <- trait(object)
    if (!length(trait))
        stop("'trait' is missing")
    ## Setting the trait. This ensures also that we don't provide wrong format.
    trait(object) <- trait
    if (!missing(traitName))
        object@traitname <- traitName
    ## Run the test.
    object@result <- .binomialTest(object, trait = trait, global = global,
                                   prob = prob, alternative = alternative)
    if (validObject(object))
        object
}

setReplaceMethod("trait", "FABinTestResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@result <- data.frame(total_phenotyped = integer(),
                                total_affected = integer(),
                                family = character(),
                                phenotyped = integer(),
                                affected = integer(),
                                pvalue = numeric(),
                                prob = numeric(),
                                check.names = FALSE, stringsAsFactors = FALSE)
    object@traitname <- character()
    if (validObject(object))
        object
})

setMethod("result", "FABinTestResults", function(object, method = "BH") {
    method <- match.arg(method, p.adjust.methods)
    traitName <- object@traitname
    if (!length(traitName))
        traitName <- NA
    res <- object@result
    padj <- numeric()
    if (nrow(res))
        padj <- p.adjust(res$pvalue, method = method)
    res <- cbind(trait_name = rep(traitName, nrow(res)),
                 res,
                 padj = padj, stringsAsFactors = FALSE)
    if (nrow(res))
        res <- res[order(res$pvalue), ]
    res
})

#' @description `.binom.test` used the `binom.test` function to calculate the
#'     p-value and supports different alternatives.
#' 
#' @md
#' @noRd
.binom.test <- function(x, size, prob, alternative = c("greater", "less",
                                                           "two.sided")) {
    alternative <- match.arg(alternative)
    mapply(x, size, FUN = function(a, b)
        binom.test(a, b, p = prob, alternative = alternative)$p.value)
}
#' @description `.pbinom.test` uses `pbinom` and `dbinom` to calculate the
#'     p-value and does only support `alternative = "greater"`. It is however
#'     considerably faster than `.binom.test`.
#'
#' @md
#' @noRd
.pbinom.test <- function(x, size, prob, alternative = "greater") {
    if (alternative != "greater")
        stop("Currently only 'alternative' 'greater' supported")
    pbinom(x, size, prob, lower.tail = FALSE) + dbinom(x, size, prob)
}

