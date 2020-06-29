
#' Plot distribution of expression values for each sample
#'
#' Use boxplots, violin plots or quantile plots (in form of parallel coordinate
#' plots of the quantiles) to compare the distribution of expression values for
#' each sample.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay to be used for
#' plotting.
#' @param method Method to use for plotting: "quantileplot" (default),
#' "boxplot", "violinplot".
#' @param coef Numeric. Used in outlier definition (median +/- coef * IQR)
#' @param title Character. Title of the plot.
#'
#' @return List with the following components:
#' \itemize{
#' \item info: data.frame with information about outlier samples or NULL
#' \item plot: plot as returned by \code{\link[GGally]{ggparcoord}} (quantile
#' plot), \code{\link[ggplot2]{ggplot}} (boxplot) or
#' \code{\link[ggpubr]{ggviolin}} (violin plot)
#' }
#'
#' @import ggplot2 ggpubr tidyr
#' @importFrom GGally ggparcoord
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' ## quantile plot
#' plot_distribution(se = se.gene,
#'                   method = "quantileplot")
#'
#' ## boxplot
#' plot_distribution(se = se.gene,
#'                   method = "boxplot")
#'
#' ## violinplot
#' library("ggpubr")
#' plot_distribution(se = se.gene,
#'                   method = "violinplot")

plot_distribution <- function(se,
                              assay = 1,
                              method = "quantileplot",
                              coef = 5,
                              title = paste("Distribution of expression",
                                            "values in each sample")) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    quant = t(apply(expr, 2, quantile, na.rm = TRUE))
    colnames(quant) = c("min", "lower", "middle", "upper", "max")

    if (method == "quantileplot") {
        prob = seq(0, 1, 0.1)
        quant.for.plot = data.frame(quantile =
                                        factor(prob,
                                               levels =
                                                   sort(prob,
                                                        decreasing = TRUE)),
                                    apply(expr, 2, quantile,
                                          probs = prob, na.rm = TRUE))
        showPoints = ifelse(ncol(se) < 50, TRUE, FALSE)

        g = ggparcoord(data = quant.for.plot,
                       columns = 2:ncol(quant.for.plot),
                       groupColumn = 1,
                       scale = "globalminmax",
                       showPoints = showPoints) +
            #            scale_color_viridis(discrete = TRUE) +
            scale_colour_viridis_d() +
            xlab("samples") +
            ylab("expression value") +
            ggtitle(title) +
            theme_pubr() +
            theme(legend.position = "right",
                  #              axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())

    } else if (method == "boxplot") {

        df = data.frame(group = rownames(quant),
                        quant)
        g = ggplot(data = df, aes_string("group")) +
            geom_boxplot(aes_string(ymin = "min",
                              lower = "lower",
                              middle = "middle",
                              upper = "upper",
                              ymax = "max"),
                         # geom_boxplot(aes_(ymin = ~min,
                         #     lower = ~lower,
                         #     middle = ~middle,
                         #     upper = ~upper,
                         #     ymax = ~max),
                         # data = df,
                         stat = "identity",
                         fill = "lightgray") +
            coord_flip() +
            labs(y = "expression value") +
            theme_pubr() +
            theme(axis.title.y = element_blank(),
                  axis.ticks.y = element_blank())


    } else if (method == "violinplot") {
        expr.long = na.omit(gather(data.frame(expr, check.names = FALSE),
                                   "sample",
                                   "value",
                                   colnames(expr)))

        g = ggviolin(data = expr.long,
                     x = "sample",
                     y = "value",
                     ylab = "expression value",
                     xlab = "",
                     add = "median_iqr",
                     fill = "lightgray",
                     orientation = "horiz")

    } else {
        stop(paste("method", method, "not known!"))
    }

#    print(g)

    ## define outliers
    out.med = get_outliers(values = quant[, c("middle")],
                           coef = coef)
    out.iqr = get_outliers(values = quant[, c("upper")] - quant[, c("lower")],
                           coef = coef)
    out.all = union(out.med, out.iqr)
    if (length(out.all) > 0) {
        if (is.null(colnames(se))) {
            stop("se needs to have colnames")
        }
        info.out = data.frame(id = colnames(se)[out.all],
                              criterion = rep("distribution", length(out.all)),
                              stringsAsFactors = FALSE)
    } else {
        info.out = NULL
    }

    return(list(info = info.out,
                plot = g))
}

