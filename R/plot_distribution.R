
#' Plot distribution of expression values for each sample
#'
#' Use boxplots, violin plots or quantile plots (in form of parallel coordinate
#' plots of the quantiles) to compare the distribution of expression values for
#' each sample.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay to be used for
#' plotting.
#' @param method Method to use for plotting: "quantileplot" (default),
#' "boxplot", "violinplot".
#' @param return.outliers Logical. Return info about outlier samples.
#' @param coef Numeric. Used in outlier definition (median +/- coef * IQR)
#'
#' @return If return.outliers = TRUE, data.frame with identifiers of outlier
#' samples
#'
#' @import ggplot2 ggpubr tidyr
#' @importFrom GGally ggparcoord
#' @export

plot_distribution <- function(se,
                              assay = 1,
                              method = "quantileplot",
                              return.outliers = TRUE,
                              coef = 5) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se)[[assay]]

    quant = t(apply(expr, 2, quantile))
    colnames(quant) = c("min", "lower", "middle", "upper", "max")

    if (method == "quantileplot") {
        prob = seq(0, 1, 0.1)
        quant.for.plot = data.frame(quantile = factor(prob),
                                    apply(expr, 2, quantile,
                                          probs = prob))
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
        expr.long = gather(data.frame(expr),
                           "sample",
                           "value",
                           colnames(expr))

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

    print(g)

    ## define outliers
    out.med = get_outliers(values = quant[, c("middle")],
                           coef = coef)
    out.iqr = get_outliers(values = quant[, c("upper")] - quant[, c("lower")],
                           coef = coef)
    out.all = union(out.med, out.iqr)
    if (length(out.all) > 0) {
        info.out = data.frame(id = colnames(se)[out.all],
                              criterion = rep("distribution", length(out.all)),
                              stringsAsFactors = FALSE)
    } else {
        info.out = NULL
    }

    if (return.outliers) return(info.out)
}

#' internal function used by plot_distribution
#'
#' determines indices of outliers in values
#'
#' @keywords internal

get_outliers <- function(values, coef = 5) {
    med = median(values, na.rm = TRUE)
    iqr = IQR(values, na.rm = TRUE)

    return(which(values < med - coef * iqr |
                     values > med + coef * iqr))
}
