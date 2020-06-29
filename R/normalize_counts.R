
#' Normalize raw read counts of RNA-seq experiments
#'
#' Generates normalized read counts of RNA-seq data sets using several
#' functions provided in the Bioconductor package \pkg{edgeR}.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing raw
#' read counts.
#' @param method Method to use for normalizing: "tmm" (default), "cpm",
#' "rpkm".
#' @param col.gene.length Character or integer. Name or number of column in
#' rowData containing gene length information.
#' @param log Logical. Should log2 transformed values be returned?
#' @param prior.count Numeric. Average count to be added to each observation
#' to avoid taking log of zero.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with normalized read counts in additional assay called <assay>.norm.
#'
#' @importFrom edgeR cpm rpkm calcNormFactors
#' @export
#'
#' @examples
#' library(recount)
#' data("rse_gene_SRP009615")
#'
#' rse_gene_SRP009615 = normalize_counts(se = rse_gene_SRP009615,
#'                                       method = "tmm")

normalize_counts <- function(se,
                             assay = 1,
                             method = "tmm",
                             col.gene.length = "bp_length",
                             log = TRUE,
                             prior.count = 2) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    counts = assays(se)[[assay]]

    if (method == "tmm") {
        norm.factors = calcNormFactors(counts,
                                       method = "TMM")
        counts.norm = cpm(y = counts,
                          log = log,
                          lib.size = norm.factors * colSums(counts),
                          prior.count = prior.count)

    } else if (method == "cpm") {
        counts.norm = cpm(y = counts,
                          lib.size = colSums(counts),
                          log = log,
                          prior.count = prior.count)

    } else if (method == "rpkm") {
        anno = rowData(se)
        if (!(col.gene.length %in% colnames(anno))) {
            stop(paste("column", col.gene.length,
                       "not available in row data!"))
        }
        gene.length = anno[, col.gene.length]
        if (!is.numeric(gene.length)) {
            stop(paste("column", col.gene.length,
                       "not numeric!"))
        }

        counts.norm = rpkm(y = counts,
                           gene.length = gene.length,
                           lib.size = NULL,
                           log = log,
                           prior.count = prior.count)
    } else {
        stop(paste("method", method, "not defined!"))
    }

    name.new = paste(ifelse(is.numeric(assay),
                            names(assays(se))[assay],
                            assay),
                     "norm", sep = ".")
    if (name.new %in% names(assays(se))) {
        warning(paste("assay", name.new, "already exists!"))
    }

    assays(se)[[name.new]] = counts.norm

    # assays.l = c(assays(se),
    #              list(counts.norm))
    # names(assays.l) = c(names(assays(se)),
    #                     name.new)
    # assays(se) = assays.l
    return(se)
}
