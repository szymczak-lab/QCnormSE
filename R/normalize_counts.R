
#' Normalize raw read counts of RNA-seq experiments
#'
#' Generates normalized read counts of RNA-seq data sets using several functions
#' provided in the Bioconductor package \pkg{DESeq2}.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay.raw Character or integer. Name or number of assay containing raw
#' read counts.
#' @param assay.norm Character. Name of new assay containing normalized
#' read counts.
#' @param method Method to use for normalizing: "vst" (default),
#' "rlog", "normTransform"
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with normalized read counts
#'
#' @import DESeq2
#' @export

normalize_counts <- function(se,
                             assay.raw = 1,
                             assay.norm = "counts.norm",
                             method = "vst") {

    ## improvements:
    ## add edgeR normalization method

    if (is.character(assay.raw) && !(assay.raw %in% names(assays(se)))) {
        stop(paste("assay", assay.raw, "not found!"))
    }
    counts = assays(se)[[assay.raw]]
    if (method == "vst") {
        counts.norm = vst(object = counts,
                          blind = TRUE)
    } else if (method == "rlog") {
        counts.norm = rlog(object = counts,
                           blind = TRUE)
        rownames(counts.norm) = rownames(counts)
    } else if (method == "normTransform") {
        dds = DESeqDataSet(se,
                           design = ~ 1)
        dds.norm = normTransform(object = dds,
                                 f = log2,
                                 pc = 1)
        counts.norm = assays(dds.norm)[[1]]
    } else {
        stop(paste("method", method, "not defined!"))
    }

    assays.l = c(assays(se),
                 list(counts.norm))
    names(assays.l) = c(names(assays(se)),
                        assay.norm)
    assays(se) = assays.l
    return(se)
}
