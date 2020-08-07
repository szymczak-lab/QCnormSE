
#' Check for batch effects
#'
#' Pairwise associations between each of the first three components of a MDS or
#' PCA analysis and defined phenotype variables is tested using F tests in a
#' linear model. P-values are visualized in a heatmap (called prince plot in
#' the R package swamp).
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param col.test Character or integer vector. Column(s) in colData() with
#' phenotype information to be tested.
#' @param res.pca List. Output of \code{\link{calculate_mds_pca}}.
#' @param title Character. Title of the plot.
#'
#' @return List with the following components:
#' \itemize{
#' \item pval: Matrix with P-values between phenotype variables in rows and
#' components in columns
#' \item r2: Matrix with absolute adjusted r^2 values between phenotype
#' variables in rows and components in columns
#' \item plot: Plot with heatmaps as returned from the
#' \code{\link[ggpubr]{ggarrange}} function
#' }
#'
#' @importFrom ggplotify as.grob
#' @importFrom ggpubr ggarrange
#' @importFrom stats lm pf
#'
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' res.pca = calculate_mds_pca(se = se.gene,
#'                             method = "pca")
#'
#' col.test = c("Age.of.patient",
#'              "Body.surface.area",
#'              "Duration.of.psoriasis",
#'              "Induration",
#'              "Overall.erythema",
#'              "Scaling",
#'              "Sex",
#'              "scan.date")
#'
#' check_batch_effects(se = se.gene,
#'                     res.pca = res.pca,
#'                     col.test = col.test)

check_batch_effects <- function(se,
                                res.pca,
                                col.test = NULL,
                                title = NULL) {

    pheno = colData(se)
    if (!is.null(col.test)) {
        if (!all(col.test %in% colnames(colData(se)))) {
            stop("not all col.test variables available in colData!")
        }
        pheno = pheno[, col.test, drop = FALSE]
    }
    pheno = pheno[, sort(colnames(pheno)), drop = FALSE]
    scores = res.pca$scores

    ## linear regression for each PC and phenotype variable
    ## P-value of F statistic and adjusted R-squared
    pval = abs.adj.r.squared = matrix(ncol = ncol(scores),
                                  nrow = ncol(pheno),
                                  dimnames = list(colnames(pheno),
                                                  colnames(scores)))
    for (i in seq_len(nrow(pval))) {
        if (length(unique(na.omit(pheno[, i]))) < 2) {
            warning(paste0("only one level found in ",
                           colnames(pheno)[i], "!\n"))
        } else {
            for (j in seq_len(ncol(pval))) {
                fit = lm(scores[, j] ~ pheno[, i])
                s = summary(fit)
                pval[i, j] = pf(q = s$fstatistic[1],
                                df1 = s$fstatistic[2],
                                df2 = s$fstatistic[3],
                                lower.tail = FALSE)
                abs.adj.r.squared[i, j] = abs(s$adj.r.squared)
            }
        }
    }

    ## remove rows with NA
    ind.na = which(apply(pval, 1, function(x) {any(is.na(x))}))
    if (length(ind.na) > 0) {
        pval = pval[-ind.na, , drop = FALSE]
        abs.adj.r.squared = abs.adj.r.squared[-ind.na, , drop = FALSE]
    }
    if (nrow(pval) == 0) {
        warning("no variables with P-values!")
        return(NULL)
    }

    ## heatmaps
    hm.pval = plot_heatmap(matrix = pval,
                           type = "pval",
                           title = title)
    hm.r2 = plot_heatmap(matrix = abs.adj.r.squared,
                         type = "r2",
                         title = title)

    ## combine
    hm = ggarrange(as_ggplot(as.grob(hm.pval)),
                   as_ggplot(as.grob(hm.r2)),
                   nrow = 1, ncol = 2)

    return(list(pval = pval,
                r2 = abs.adj.r.squared,
                plot = hm))
}



#' Plot heatmap
#'
#' Plot heatmap with predefined color scheme based on type.
#'
#' @param matrix Matrix. Values to be plotted.
#' @param type Character. Type of heatmap ('pval' for P-values and
#' 'r2' for absolute adjusted r^2).
#' @param title Character. Title of the plot.
#'
#' @return Heatmap as \code{\link[ComplexHeatmap]{Heatmap-class}} object
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#'
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' res.pca = calculate_mds_pca(se = se.gene,
#'                             method = "pca")
#'
#' col.test = c("Age.of.patient",
#'              "Body.surface.area",
#'              "Duration.of.psoriasis",
#'              "Induration",
#'              "Overall.erythema",
#'              "Scaling",
#'              "Sex",
#'              "scan.date")
#'
#' res = check_batch_effects(se = se.gene,
#'                           res.pca = res.pca,
#'                           col.test = col.test)
#'
#' # plot P-values
#' plot_heatmap(matrix = res$pval,
#'              type = "pval")
#'
#' # plot absolute adjusted r^2
#' plot_heatmap(matrix = res$r2,
#'              type = "r2")

plot_heatmap <- function(matrix,
                         type,
                         title = NULL) {

    if (type == "pval") {
        ## categorize P-values
        breaks = c(0, 10^-4, 10^-3, 10^-2, 0.05, 1)
        matrix = apply(matrix, c(1, 2), function(x) {
            cut(x, breaks = breaks)
        })

        col = c("white", "yellow", "orange", "red", "darkred")
        names(col) = c("(0.05,1]",
                       "(0.01,0.05]",
                       "(0.001,0.01]",
                       "(0.0001,0.001]",
                       "(0,0.0001]")

        name = "P-value"

    } else if (type == "r2") {
        col = colorRamp2(c(0, 1),
                         c("white", "blue"))
        name = "|r^2|"

    } else {
        stop(paste("type", type, "not known! Use 'pval' or 'r2."))
    }

    hm = Heatmap(matrix,
                 rect_gp = gpar(col = "black"),
                 name = name,
                 col = col,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 heatmap_legend_param = list(border = "black"),
                 row_names_side = "left",
                 column_names_side = "top",
                 column_names_rot = 0,
                 column_title = title)

    return(hm)

}


#' Define batches based on scan date
#'
#' Groups samples into batches based on scan date. Samples run within a short
#' time interval can be defined to belong to the same batch.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param col.scan.date Character or Integer. Column in colData() with scan
#' dates (default: scan.date as generated by the function
#' \code{\link{extract_scan_date}}.
#' @param diff.ignore Numeric. Time difference (in days) defined as negligible
#' (default: 1).
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with batch information in additional column in colData()
#'
#' @export
#'
#' @examples
#' data("se.gene")
#' se.gene = define_batches(se = se.gene,
#'                          col.scan.date = "scan.date")

define_batches <- function(se,
                           col.scan.date = "scan.date",
                           diff.ignore = 1) {

    if (!(col.scan.date %in% colnames(colData(se)))) {
        stop(paste("column", col.scan.date, "not available in column data!"))
    }

    if ("batch" %in% colnames(colData(se))) {
        stop(paste("column batch already available in column data!"))
    }

    scan.date = colData(se)[, col.scan.date]
    if (!is(scan.date, "Date")) {
        scan.date = as.Date(scan.date)
    }
    names(scan.date) = colnames(se)

    scan.date = sort(scan.date)
    diff = c(0, as.numeric(diff(scan.date)))
    diff[diff == diff.ignore] = 0
    batch = as.numeric(as.factor(cumsum(diff)))
    names(batch) = names(scan.date)
    se$batch = batch[colnames(se)]
    return(se)
}


