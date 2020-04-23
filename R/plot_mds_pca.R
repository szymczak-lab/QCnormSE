
#' MDS or PCA
#'
#' Performs dimensionality reduction using multi-dimensional scaling (MDS) or
#' principal component analysis (PCA) on the sample level.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' expression data to be used for dimensionality reduction.
#' @param method Method to use: "pca" (default), "mds".
#' @param dist Distance method to be used for MDS: see method argument in
#' \code{\link[stats]{dist}} function.
#' @param center Logical. Should the data be centered by mean? (default: TRUE).
#' @param scale Logical. Should the data be scaled by standard deviation?
#' (default: FALSE).
#'
#' @return data.frame with components of PCA or MDS analysis.
#'
#' @importFrom stats cmdscale dist na.omit prcomp
#' @export

calculate_mds_pca <- function(se,
                              assay = 1,
                              method = "pca",
                              dist = "euclidean",
                              center = TRUE,
                              scale = FALSE) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    if (center) {
        mean = apply(expr, 1, mean)
        expr = expr - mean
    }
    if (scale) {
        sd = apply(expr, 1, sd)
        expr = expr / sd
    }

    if (method == "pca") {
        res = prcomp(na.omit(t(expr)),
                     center = FALSE,
                     scale = FALSE)
        scores = res$x[, 1:3]

    } else if (method == "mds") {

        # if (dist == "gower") {
        #     d = daisy(t(expr), metric = "gower", stand = FALSE)
        # } else {
        #     d = Dist(t(expr), method = dist)
        # }
        d = dist(t(expr))
        fit = cmdscale(d, eig = TRUE, k = 3)
        scores = fit$points

    } else {
        stop(paste("method", method, "not known!"))
    }
    scores = data.frame(scores)
#    colnames(scores) = paste("component",
#                             seq_len(ncol(scores)),
#                             sep = ".")
    return(scores)
}

#' MDS or PCA plot
#'
#' Plots components estimated with the function \code{\link{calculate_mds_pca}}.
#' Color and shape of each sample can be set based on different variables.
#'
#' @param scores data.frame. Components of PCA or MDS analysis as estimated by
#' \code{\link{calculate_mds_pca}}.
#' @param var.color Character or integer vector. Variable used to determine
#' color. If NULL black color will be used for all samples.
#' @param palette Color palette to be used (default palette from
#' \code{\link[ggpubr]{get_palette}}).
#' @param var.shape Character or integer vector. Variable used to determine
#' shape. If NULL filled circles will be used for all samples.
#' @param title Character. Title of the plot. If NULL title will be set based on
#' method.
#' @param return.outliers Logical. Return info about outlier samples.
#' @param factor Numeric. Parameter of the function
#' \code{\link[aplpack]{compute.bagplot}}. (default: 5)
#'
#' @return If return.outliers = TRUE, data.frame with identifiers of outlier
#' samples for each pairwise combination of components.
#'
#' @import ggpubr ggplot2
#' @importFrom aplpack compute.bagplot
#' @export

plot_mds_pca <- function(scores,
                         var.color = NULL,
                         palette = NULL,
                         var.shape = NULL,
                         title = NULL,
                         return.outliers = TRUE,
                         factor = 5) {

    ## set title to method (set based on colnames of scores)
    method = ifelse(grepl("^PC", colnames(scores)[1]),
                    "PCA", "MDS")
    if (is.null(title)) title = method

    ## plots
    plots.l = list(plot_mds_pca_2d(scores = scores,
                                   dim = 1:2,
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape),
                   plot_mds_pca_2d(scores = scores,
                                   dim = c(1, 3),
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape),
                   plot_mds_pca_2d(scores = scores,
                                   dim = 2:3,
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape))

    legend = ifelse(is.null(var.color) & is.null(var.shape),
                    "none", "bottom")
    g = ggarrange(plotlist = plots.l,
                  ncol = 3,
                  nrow = 1,
                  legend = legend,
                  common.legend = TRUE)

    g = annotate_figure(g,
                        top = text_grob(title,
                                        face = "bold",
                                        size = 14))
    print(g)


    if (return.outliers) {
        info.out = rbind(get_outliers_mds_pca_2d(scores = scores,
                                                 dim = 1:2,
                                                 factor = factor),
                         get_outliers_mds_pca_2d(scores = scores,
                                                 dim = c(1, 3),
                                                 factor = factor),
                         get_outliers_mds_pca_2d(scores = scores,
                                                 dim = 2:3,
                                                 factor = factor))
        return(info.out)
    }

}

#' MDS or PCA plot (2D)
#'
#' Plots two sepcified components estimated with the function
#' \code{\link{calculate_mds_pca}}. Color and shape of each sample can be set
#' based on different variables.
#'
#' @param scores data.frame. Components of PCA or MDS analysis as estimated by
#' \code{\link{calculate_mds_pca}}.
#' @param dim Numeric vector (2). Numbers of components to be plotted.
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param var.color Character. Name of variable used to determine
#' color. If NULL black color will be used for all samples.
#' @param palette Color palette to be used (default palette used in
#' \code{\link[ggpubr]{ggscatter}}).
#' @param var.shape Character. Name of variable used to determine
#' shape. If NULL filled circles will be used for all samples.
#'
#' @import ggpubr ggplot2
#' @export

plot_mds_pca_2d <- function(scores,
                            dim,
                            se,
                            var.color = NULL,
                            palette = NULL,
                            var.shape = NULL) {

    xlab = paste("Component", dim[1])
    ylab = paste("Component", dim[2])

    ## set color and shape
    info = data.frame(prepare_var_for_plot(se = se,
                                           var = var.color),
                      prepare_var_for_plot(se = se,
                                           var = var.shape))
    colnames(info) = c(ifelse(is.null(var.color), "color", var.color),
                       ifelse(is.null(var.shape), "shape", var.shape))
    scores = cbind(info, scores)

    g = ggscatter(scores,
                  x = colnames(scores)[dim[1] + 2],
                  y = colnames(scores)[dim[2] + 2],
                  color = colnames(scores)[1],
                  shape = colnames(scores)[2],
                  size = 3,
                  xlab = xlab,
                  ylab = ylab)

    if (is.null(var.color)) {
        g = g + guides(color = FALSE)
    }
    if (is.null(var.shape)) {
        g = g + guides(shape = FALSE)
    }
    ## set colors
    if (!is.null(palette)) {
        if (is.factor(scores[, 1])) {
            g = change_palette(g, palette)
        } else {
            g = g + gradient_color(palette)
        }
    } else {
        if (length(unique(scores[, 1])) == 1) {
            g = change_palette(g, "black")
        }
    }

    return(g)

}



#' internal function used by plot_mds_pca
#'
#' if var = NULL returns vector of 1s
#'
#' @keywords internal

prepare_var_for_plot <- function(se,
                                 var = NULL) {

    if (is.null(var) || !(var %in% colnames(colData(se)))) {
        values = rep(1, ncol(se))

    } else {
        values = colData(se)[, var]
    }

    if (length(unique(values)) <= 10) {
        values = as.factor(values)
    }

    return(values)
}


#' internal function used by plot_mds_pca
#'
#' determines outliers in two dimensional plot
#'
#' @importFrom aplpack compute.bagplot
#' @keywords internal

get_outliers_mds_pca_2d <- function(scores,
                                    dim,
                                    factor = 5) {

    method = ifelse(grepl("^PC", colnames(scores)[1]),
                    "PCA", "MDS")
    criterion = paste(method, "dim", paste(dim, collapse = "."), sep = ".")

    info.bagplot = compute.bagplot(x = scores[, dim[1]],
                                            y = scores[, dim[2]],
                                            factor = factor)
    out = info.bagplot$pxy.outlier
    if (!is.null(out)) {

        ind.out = apply(out, 1, function(x) {
            which(scores[, dim[1]] == x[1] & scores[, dim[2]] == x[2])
        })
        info.out = data.frame(id = rownames(scores)[ind.out],
                              criterion = rep(criterion, length(ind.out)),
                              stringsAsFactors = FALSE)
    } else {
        info.out = NULL
    }
    return(info.out)
}
