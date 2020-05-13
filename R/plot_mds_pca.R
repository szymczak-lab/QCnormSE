
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
#' @return list with components scores (data.frame with components of PCA or
#' MDS analysis) and var.explained (vector with explained variance; only for
#' PCA).
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
        var.explained = (res$sdev^2)[1:3] / sum(res$sdev^2)

    } else if (method == "mds") {

        # if (dist == "gower") {
        #     d = daisy(t(expr), metric = "gower", stand = FALSE)
        # } else {
        #     d = Dist(t(expr), method = dist)
        # }
        d = dist(t(expr))
        fit = cmdscale(d, eig = TRUE, k = 3)
        scores = fit$points
        var.explained = NULL

    } else {
        stop(paste("method", method, "not known!"))
    }
    scores = data.frame(scores)
#    colnames(scores) = paste("component",
#                             seq_len(ncol(scores)),
#                             sep = ".")
    return(list(scores = scores,
                var.explained = var.explained))
}

#' MDS or PCA plot
#'
#' Plots components estimated with the function \code{\link{calculate_mds_pca}}.
#' Color and shape of each sample can be set based on different variables.
#'
#' @param res data.frame. Output of \code{\link{calculate_mds_pca}}.
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param var.color Character or integer vector. Variable used to determine
#' color. If NULL black color will be used for all samples.
#' @param palette Color palette to be used (default palette from
#' \code{\link[ggpubr]{get_palette}}).
#' @param var.shape Character or integer vector. Variable used to determine
#' shape. If NULL filled circles will be used for all samples.
#' @param shape.values Vector with symbols. Needs to provide a symbol for each
#' unique value of var.shape.
#' @param title Character. Title of the plot. If NULL title will be set based on
#' method.
#' @param return.outliers Logical. Return info about outlier samples.
#' @param factor Numeric. Parameter of the function
#' \code{\link[aplpack]{compute.bagplot}}. (default: 5)
#' @param ellipse Logical. Should ellipses around points be drawn? (default:
#' FALSE).
#' @param ellipse.type Character. Type of ellipse as given in
#'
#' @return If return.outliers = TRUE, data.frame with identifiers of outlier
#' samples for each pairwise combination of components.
#'
#' @import ggpubr ggplot2
#' @importFrom aplpack compute.bagplot
#' @export

plot_mds_pca <- function(res,
                         se,
                         var.color = NULL,
                         palette = NULL,
                         var.shape = NULL,
                         shape.values = NULL,
                         title = NULL,
                         return.outliers = TRUE,
                         factor = 5,
                         ellipse = FALSE,
                         ellipse.type = "convex") {

    ## set title to method (set based on colnames of scores)
#    method = ifelse(grepl("^PC", colnames(scores)[1]),
    method = ifelse(!is.null(res$var.explained), "PCA", "MDS")
    if (is.null(title)) title = method

    ## plots
    plots.l = list(plot_mds_pca_2d(res = res,
                                   se = se,
                                   dim = 1:2,
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape,
                                   shape.values = shape.values,
                                   title = "",
                                   ellipse = ellipse,
                                   ellipse.type = ellipse.type),
                   plot_mds_pca_2d(res = res,
                                   se = se,
                                   dim = c(1, 3),
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape,
                                   shape.values = shape.values,
                                   title = "",
                                   ellipse = ellipse,
                                   ellipse.type = ellipse.type),
                   plot_mds_pca_2d(res = res,
                                   se = se,
                                   dim = 2:3,
                                   var.color = var.color,
                                   palette = palette,
                                   var.shape = var.shape,
                                   shape.values = shape.values,
                                   title = "",
                                   ellipse = ellipse,
                                   ellipse.type = ellipse.type))

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
        info.out = rbind(get_outliers_mds_pca_2d(scores = res$scores,
                                                 dim = 1:2,
                                                 factor = factor),
                         get_outliers_mds_pca_2d(scores = res$scores,
                                                 dim = c(1, 3),
                                                 factor = factor),
                         get_outliers_mds_pca_2d(scores = res$scores,
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
#' @param res data.frame. Output of \code{\link{calculate_mds_pca}}.
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
#' @param shape.values Vector with symbols. Needs to provide a symbol for each
#' unique value of var.shape.
#' @param title Character. Main title of plot.
#' @param ellipse Logical. Should ellipses around points be drawn? (default:
#' FALSE).
#' @param ellipse.type Character. Type of ellipse as given in
#' \code{\link[ggpubr]{ggscatter}} (default: "convex").
#'
#' @import ggpubr ggplot2
#' @export

plot_mds_pca_2d <- function(res,
                            dim,
                            se,
                            var.color = NULL,
                            palette = NULL,
                            var.shape = NULL,
                            shape.values = NULL,
                            title = NULL,
                            ellipse = FALSE,
                            ellipse.type = "convex") {
    scores = res$scores

    xlab = paste("Component", dim[1])
    ylab = paste("Component", dim[2])
    if (!is.null(res$var.explained)) {
        var = round(res$var.explained * 100, digits = 2)
        xlab = paste0(xlab, " (", var[dim[1]], "%)")
        ylab = paste0(ylab, " (", var[dim[2]], "%)")

        if (is.null(title)) {
            title = "PCA"
        }
    } else {
        if (is.null(title)) {
            title = "MDS"
        }
    }

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
                  ylab = ylab,
                  title = title,
                  ellipse = ellipse,
                  ellipse.type = ellipse.type)

    if (!is.null(shape.values)) {
        g = g + scale_shape_manual(values = shape.values)
    }

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

    if (is.null(var)) {
        values = rep(1, ncol(se))

    } else {
        if (!(var %in% colnames(colData(se)))) {
            stop(paste("var", var, "not found in colnames of colData!"))
        }
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
