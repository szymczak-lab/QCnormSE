
# internal function used by plot_distribution
#
# determines indices of outliers in values
#
#' @keywords internal

get_outliers <- function(values, coef = 5, tail = "both") {
    med = median(values, na.rm = TRUE)
    iqr = IQR(values, na.rm = TRUE)

    lower = which(values < med - coef * iqr)
    upper = which(values > med + coef * iqr)

    if (!(tail %in% c("both", "upper", "lower"))) {
        stop("tail must be set to either 'both', 'upper' or 'lower'!")
    }
    res = switch(tail,
                 both = union(upper, lower),
                 upper = upper,
                 lower = lower)
    return(res)
}
