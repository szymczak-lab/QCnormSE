
#' Log transformation
#'
#' Performs log transformation of expression values (log2). Values below zero 
#' will be set to zero and if zero values are available 
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' expression data to be log transformed.
#' @param pseudocount Numeric. Value added to expression values before log
#' transformation (relevant for count data containing zero values).
#' 
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object.
#' 
#' @export

log_transform <- function(se, 
                          assay = 1,
                          pseudocount = 1) {
    
    expr = assays(se)[[assay]]
    
    if (any(expr < 0)) {
        expr[expr < 0] = 0
    }
    if (any(expr == 0)) {
        expr = expr + pseudocount
    }
    expr = log2(expr)
    if (any(is.na(expr) | is.infinite(expr))) {
        stop("missing expression values after log transformation!")
    }
    assays(se)[[assay]] = expr
    
    return(se)
}
