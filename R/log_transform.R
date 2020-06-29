
#' Log transformation
#'
#' Performs log transformation of expression values (log2). Values below zero
#' will be set to zero and if zero values are available a pseudocount will be
#' added.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' expression data to be log transformed.
#' @param pseudocount Numeric. Value added to expression values before log
#' transformation (relevant for count data containing zero values).
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with normalized log transformed expression values in additional assay
#' called <assay>.log.
#'
#' @export
#'
#' @examples
#' data("se.probeset")
#'
#' se.probeset = log_transform(se = se.probeset)

log_transform <- function(se,
                          assay = 1,
                          pseudocount = 1) {

    expr = assays(se)[[assay]]

    if (any(expr < 0, na.rm = TRUE)) {
        expr[expr < 0] = 0
    }
    if (any(expr == 0, na.rm = TRUE)) {
        expr = expr + pseudocount
    }
    expr.log = log2(expr)
    if (any(is.infinite(expr.log))) {
        stop("infinite expression values after log transformation!")
    }
    if (any(is.na(expr.log)) && sum(is.na(expr.log)) > sum(is.na(expr))) {
        stop("missing expression values after log transformation!")
    }

    name.new = paste(ifelse(is.numeric(assay),
                            names(assays(se))[assay],
                            assay),
                     "log", sep = ".")
    if (name.new %in% names(assays(se))) {
        warning(paste("assay", name.new, "already existed!"))
    }

    assays(se)[[name.new]] = expr.log

    return(se)
}
