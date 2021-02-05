
#' Remove samples based on different criteria
#'
#' Remove samples with large frequency of missing values
#' (method = "missing") or low expression (method = "detection.pvalue" or
#' method = "zero").
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay to be used for
#' filtering.
#' @param method Method to determine samples to be removed: "missing",
#' "detection.pvalue", "zero".
#' @param freq Numeric. If more than freq*100 \% of the genes fulfill criterion,
#' the sample is removed.
#' @param verbose Logical. Should number of removed samples be reported?
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with samples removed
#'
#' @export
#'
#' @examples
#' library(recount)
#' data("rse_gene_SRP009615")
#'
#' rse_gene_SRP009615 = remove_samples(se = rse_gene_SRP009615,
#'                                     method = "zero",
#'                                     freq = 0.4)

remove_samples <- function(se,
                           assay = 1,
                           method,
                           freq = 0.5,
                           verbose = FALSE) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    if (freq < 0 | freq > 1) {
        stop("freq needs to be between 0 and 1 (inclusive)!")
    }

    expr = assays(se)[[assay]]

    if (method == "missing") {
        crit = apply(expr, 2, function(x) {
            sum(is.na(x))}) / nrow(expr)
    } else if (method == "zero") {
        crit = apply(expr, 2, function(x) {
            sum(x == 0, na.rm = TRUE)}) / nrow(expr)
    } else if (method == "detection.pvalue") {
        range = range(as.numeric(expr), na.rm = TRUE)
        if (min(range) < 0 | max(range) > 1) {
            stop(paste("assay", assay, "does not contain P values!"))
        }
        crit = apply(expr, 2, function(x) {
            sum(x >= 0.05, na.rm = TRUE)}) / nrow(expr)
    } else {
        stop(paste("method", method, "not known!"))
    }

    if (freq == 1) {
        ind.rm = which(crit == freq)
    } else {
        ind.rm = which(crit > freq)
    }

    if (length(ind.rm) > 0) {
        if (verbose) {
            print(paste(length(ind.rm), "samples removed"))
        }
        se = se[, -ind.rm]
    }
    return(se)

}
