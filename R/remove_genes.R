
#' Remove genes based on different criteria
#'
#' Remove genes with no variability (method = "constant"), missing values
#' (method = "missing") or low expression (method = "detection.pvalue",
#' method = "zero" or method = "edgeR").
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay to be used for
#' filtering.
#' @param method Method to determine genes to be removed: "constant", "missing",
#' "detection.pvalue", "zero", "edgeR" (using
#' \code{\link[edgeR]{filterByExpr}}).
#' @param freq Numeric. If more than freq*100 \% of the samples fulfill
#' criterion, the gene is removed (not used if method = "constant").
#' @param verbose Logical. Should number of removed genes be reported?
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with genes removed
#'
#' @importFrom edgeR filterByExpr
#' @export
#'
#' @examples
#' data("se.probeset")
#' print(se.probeset)
#'
#' ## remove probesets with missing values
#' se.probeset.red = remove_genes(se = se.probeset,
#'                                assay = "exprs.log",
#'                                method = "missing",
#'                                freq = 0)
#' print(se.probeset.red)
#'
#' ## remove probesets with low expression
#' se.probeset.red = remove_genes(se = se.probeset,
#'                                assay = "detection.pvalue",
#'                                method = "detection.pvalue",
#'                                freq = 0.75)
#' print(se.probeset.red)

remove_genes <- function(se,
                         assay = 1,
                         method,
                         freq = 0.25,
                         verbose = FALSE) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    if (freq < 0 | freq > 1) {
        stop("freq needs to be between 0 and 1 (inclusive)!")
    }

    expr = assays(se)[[assay]]

    if (method == "constant") {
        sd = apply(expr, 1, sd)
        ind.rm = which(sd == 0)
    } else if (method == "edgeR") {
        n = floor(ncol(se) * freq)
        group = c(rep(1, n), rep(2, ncol(se) - n))
        info.keep = filterByExpr(y = expr,
                                 group = group,
                                 min.count = 10,
                                 min.total.count = 15)
        ind.rm = which(!info.keep)
    } else {
        if (method == "missing") {
            crit = apply(expr, 1, function(x) {sum(is.na(x))}) / ncol(expr)
        } else if (method == "zero") {
            crit = apply(expr, 1, function(x) {sum(x == 0)}) / ncol(expr)
        } else if (method == "detection.pvalue") {
            range = range(as.numeric(expr), na.rm = TRUE)
            if (min(range) < 0 | max(range) > 1) {
                stop(paste("assay", assay, "does not contain P values!"))
            }
            crit = apply(expr, 1, function(x) {sum(x >= 0.05)}) / ncol(expr)
        } else {
            stop(paste("method", method, "not known!"))
        }

        if (freq == 1) {
            ind.rm = which(crit == freq)
        } else {
            ind.rm = which(crit > freq)
        }
    }
    if (length(ind.rm) > 0) {
        if (verbose) {
            print(paste(length(ind.rm), "genes removed"))
        }
        se = se[-ind.rm, ]
    }
    return(se)

}
