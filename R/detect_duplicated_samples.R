
#' Identification of duplicated samples
#'
#' Detects duplicated samples with extremely high correlation of expression
#' values using the function \code{\link[doppelgangR]{outlierFinder}}.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' log transformed expression values.
#' @param cor.method Character. Correlation coefficient to be used ("pearson",
#' "spearman" (default) or "kendall").
#' @param ... Additional arguments passed to the
#' \code{\link[doppelgangR]{outlierFinder}} function.
#'
#' @return List with the following components:
#' \itemize{
#' \item info: data.frame with identifiers of duplicated samples and their
#' correlation or NULL
#' \item plot: Histogram as returned by \code{\link[ggpubr]{gghistogram}}
#' }
#'
#' @importFrom doppelgangR outlierFinder
#' @importFrom methods as
#' @importFrom stats cor
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' detect_duplicated_samples(se = se.gene)

detect_duplicated_samples <- function(se,
                                      assay = 1,
                                      cor.method = "spearman",
                                      ...) {

    ## convert to ExpressionSet using only specified assay
    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    ## estimate pairwise correlation
    cor.m = cor(expr,
                method = cor.method,
                use = "pairwise.complete.obs")

    ## outlier detection using function in doppelgangR package
    cor.m[lower.tri(cor.m)] = NA
    diag(cor.m) = NA

    res.doppel = outlierFinder(similarity.mat = cor.m,
                               ...)
    ind.out = which(res.doppel$outlierFinder.res$doppel == TRUE)

    if (length(ind.out) > 0) {
        info.out = res.doppel$outlierFinder.res[ind.out, , drop = FALSE]
    } else {
        info.out = NULL
    }

    ## use tidyr to convert matrix to data.frame with info for each pair
    #https://stackoverflow.com/questions/45825685/correlations-for-pairs-of-combinations
    info.cor = data.frame(row = rownames(cor.m)[row(cor.m)[upper.tri(cor.m)]],
                          col = colnames(cor.m)[col(cor.m)[upper.tri(cor.m)]],
                          cor = cor.m[upper.tri(cor.m)],
                          stringsAsFactors = FALSE)

    g = gghistogram(data = info.cor$cor,
                    y = "..density..",
                    xlab = "pairwise correlations",
                    fill = "lightgray",
                    bins = 20,
                    alpha = 0,
                    add_density = TRUE,
                    rug = TRUE)

    if (length(ind.out) > 0) {
        g = g + geom_vline(aes(xintercept = info.out$similarity),
                           color = "red",
                           size = 0.75)
    }
#    print(g)

    return(list(info = info.out,
                plot = g))
}
