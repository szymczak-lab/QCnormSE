
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
#' @param title Character. Title of the plot.
#' @param use.fast Logical. Should fast implementation of covariance matrix
#' estimation provided in the R package coop (default: TRUE).
#' @param remove.genes Logical. Should genes with missing values be removed?
#' (default: TRUE).
#' @param cor.method Character. Correlation coefficient to be used ("pearson",
#' "spearman" or "kendall"). Ignored if use.fast = TRUE.
#' @param use Character. Method to deal with missing values as in
#' \code{\link[stats]{cor}} and \code{\link[coop]{covar}} (default:
#' "everything").
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
#' @importFrom coop covar
#' @importFrom doppelgangR outlierFinder
#' @importFrom methods as
#' @importFrom stats cor sd
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' detect_duplicated_samples(se = se.gene)

detect_duplicated_samples <- function(se,
                                      assay = 1,
                                      title = NULL,
                                      use.fast = TRUE,
                                      remove.genes = TRUE,
                                      cor.method = "pearson",
                                      use = "everything",
                                      ...) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    ## remove genes with missing values
    if (remove.genes && any(is.na(expr))) {
        se.reduced = remove_genes(se = se,
                                  assay = assay,
                                  method = "missing",
                                  freq = 0)
        expr = assays(se.reduced)[[assay]]
        print(paste("removed", nrow(se) - nrow(se.reduced),
                    "genes with missing values!"))
    }

    ## estimate pairwise correlation
    if (use.fast) {
        expr.std = apply(expr, 2, function(x) {x / sd(x, na.rm = TRUE)})
        cor.m = covar(expr.std,
                      use = use)
        dimnames(cor.m) = list(colnames(expr.std),
                               colnames(expr.std))
    } else {
        cor.m = cor(expr,
                    method = cor.method,
                    use = use)
    }

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
                    title = title,
                    fill = "lightgray",
                    bins = 20,
                    alpha = 0,
                    add_density = TRUE,
                    rug = TRUE)

    if (length(ind.out) > 0) {
        g = g + geom_vline(xintercept = info.out$similarity,
                           color = "red",
                           size = 0.75)
    }
#    print(g)

    return(list(info = info.out,
                plot = g))
}
