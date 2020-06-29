
#' Aggregate expression values based on new identifiers
#'
#' For statistical analysis original gene identifiers (e.g. vendor specific
#' probe set identifiers) often need to be mapped to new gene identifiers (e.g.
#' Ensembl gene identifiers or HGNC gene symbols). This function aggregates
#' expression values of original idenfiers that map to the same new gene
#' identifier by e.g. selecting the one with the largest average expression
#' across all samples.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay used for
#' aggregating.
#' @param col.new Character or integer. Name or number of column in rowData
#' to be used as new gene identifier.
#' @param sep Character. Separator for multiple gene identifiers or names
#' (default: "///" used by GEO).
#' @param method Method to use for aggregating: "max_median" (default),
#' "max_mean", "sum"
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with aggregated expression values
#'
#' @import SummarizedExperiment
#' @importFrom stats IQR median quantile
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' data("se.probeset")
#'
#' ## restrict to subset of probesets (for illustration only)
#' genes = c("DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y",
#'           "RPS4Y1", "TXLNG2P", "UTY", "XIST")
#' ind = unlist(sapply(genes, function(g) {
#'     grep(g, rowData(se.probeset)$Gene.symbol)}))
#' se.probeset = se.probeset[ind, ]
#' print(se.probeset)
#'
#' ## aggregate by gene symbol
#' se.gene = aggregate_by_new_id(se = se.probeset,
#'                               col.new = "Gene.symbol",
#'                               sep = "///")
#' print(se.gene)

# @seealso \code{\link{get_annotation}}

aggregate_by_new_id <- function(se,
                                assay = 1,
                                col.new = "symbol",
                                sep = "///",
                                method = "max_median") {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    if (!(col.new %in% colnames(rowData(se)))) {
        stop(paste("column", col.new, "not found in rowData!"))
    }
    anno = rowData(se)

    ## remove rows with missing or multiple information in col.new
    ind.na = which(is.na(anno[, col.new]) | anno[, col.new] == "")
    ind.multiple = grep(sep, anno[, col.new])
    ind.rm = union(ind.na, ind.multiple)
    if (length(ind.rm) > 0) {
        anno = anno[-ind.rm, ]
        expr = expr[-ind.rm, ]
    }

    ## extract info for unique genes
    tab = table(anno[, col.new])
    genes.unique = names(tab)[tab == 1]
    genes.unique.original = rownames(anno)[which(anno[, col.new] %in%
                                                     genes.unique)]
    expr.new.unique = expr[genes.unique.original, ]
    id.new.unique = genes.unique.original

    ## aggregate info for non-unique genes
    genes.new = names(tab)[tab > 1]
    no.genes = length(genes.new)

    expr.new = matrix(nrow = no.genes,
                      ncol = ncol(expr))
    id.new = rep(NA, no.genes)
    for (i in seq_len(no.genes)) {
        genes.original = rownames(expr)[which(anno[, col.new] == genes.new[i])]
        expr.temp = expr[genes.original, , drop = FALSE]

        if (grepl("max", method)) {
            if (method == "max_median") {
                avg = apply(expr.temp, 1, median, na.rm = TRUE)
            } else if (method == "max_mean") {
                avg = apply(expr.temp, 1, mean, na.rm = TRUE)
            } else {
                stop(paste("method", method, "not defined!"))
            }
            ind.max = which.max(avg)
            expr.new[i, ] = expr.temp[ind.max, ]
            id.new[i] = genes.original[ind.max]
        } else if (method == "sum") {
            expr.new[i, ] = apply(expr.temp, 2, sum)

            ## possible improvement: keep all original identifiers and
            ## corresponding annotation
            id.new[i] = genes.original[1]
        } else {
            stop(paste("method", method, "not defined!"))
        }
    }
    expr.new = rbind(expr.new.unique,
                     expr.new)
    id.new = c(id.new.unique, id.new)
    dimnames(expr.new) = list(id.new,
                              colnames(expr))

    assays.list = list(expr.new)
    name = ifelse(is.numeric(assay),
                  names(assays(se))[assay], assay)
    names(assays.list) = name
    se.new = SummarizedExperiment(assays = assays.list,
                                  colData = colData(se),
                                  rowData =
                                      rowData(se)[id.new, ])

    rownames(se.new) = rowData(se.new)[, col.new]
    return(se.new)
}
