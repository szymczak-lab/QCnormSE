
#' Merge SummarizedExperiments objects
#'
#' Combine several SummarizedExperiments objects into a single
#' SummarizedExperiments object based on common genes and common colnames in
#' colData.
#'
#' @param se.l
#' list of \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' objects
#' @param info Vector of same length as se.l. Contains additional information
#' about each SE object (e.g. study name) which will be replicated for each
#' sample in the corresponding SE object and stored in a column named info.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with all samples but no rowData() information
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' data("se.gene")
#'
#' ## split SE object into two SE objects (for illustration only)
#' ## and store in list
#' se.l = list(part.1 = se.gene[, 1:10],
#'             part.2 = se.gene[, 11:26])
#'
#' se.combined = combine_se_objects(se.l = se.l)
#' all.equal(dim(se.gene), dim(se.combined))
#'
#' ## add information about parts
#' se.combined = combine_se_objects(se.l = se.l,
#'                                  info = c("part.1", "part.2"))
#' table(se.combined$info)

combine_se_objects <- function(se.l,
                               info = NULL) {

    ## identify common genes
    genes = unlist(lapply(se.l, rownames))
    tab = table(genes)
    genes.all = names(tab)[tab == length(se.l)]
    if (length(genes.all) == 0) {
        stop("no overlapping genes found!")
    }

    ## identify common columns in colData
    colnames = unlist(lapply(se.l, function(x) {
        colnames(colData(x))}))
    tab = table(colnames)
    colnames.all = names(tab)[tab == length(se.l)]
    if (length(colnames.all) == 0) {
        stop("no overlapping colnames in colData found!")
    }

    ## merge SEs
    se.all = NULL
    for (i in seq_len(length(se.l))) {
        se = se.l[[i]][genes.all, ]
        colData(se) = colData(se)[, colnames.all]
        rowData(se) = NULL
        metadata(se) = list()
        if (i == 1) {
            se.all = se
        } else {
            se.all = cbind(se.all, se)
        }
    }

    if (!is.null(info)) {
        if (length(info) != length(se.l)) {
            stop("info and se.l need to have the same length!")
        }
        info = rep(info,
                   times = vapply(se.l,
                                  FUN = ncol,
                                  FUN.VALUE = numeric(1)))
        if ("info" %in% colnames.all) {
            warning("column info will be overwritten!")
        }
        se.all$info = info
    }

    metadata(se.all) = list()
    return(se.all)
}
