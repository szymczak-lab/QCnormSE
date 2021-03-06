
#' Extract GEO information
#'
#' Extracts additional sample information from GEO into separate columns of
#' colData(). Numeric variables are converted from character
#' to numeric and syntactically valid variable names are used for the new
#' columns.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param cols Character or integer. Name or number of column(s) in colData
#' with GEO information (optional, will otherwise be extracted assuming column
#' names contain term 'characteristics').
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with GEO information in additional columns in colData()
#'
#' @importFrom methods is
#' @export
#'
#' @examples
#' library(recount)
#' data("rse_gene_SRP009615")
#'
#' rse_gene_SRP009615 = extract_geo_characteristics(se = rse_gene_SRP009615)

extract_geo_characteristics <- function(se, cols = NULL) {

    pheno = colData(se)

    ## columns with information about characteristics
    if (is.null(cols)) cols = grep("characteristics", colnames(pheno))

    ## GEO information is stored as list in single column in SE objects from
    ## recount (will be converted to matrix)
    if (length(cols) == 1 && is(pheno[, cols], "CharacterList")) {
        info.geo.original = t(sapply(pheno[, cols], unlist))
    } else {
        info.geo.original = pheno[cols, , drop = FALSE]
    }

    if (all(is.na(info.geo.original))) {
        warning("no GEO information available!")
    } else {

        ## extract information
        info.geo = do.call(cbind, apply(info.geo.original, 2, convert))

        ## convert colnames to syntactically valid variable names
        colnames(info.geo) = base::make.names(colnames(info.geo))

        colData(se) = cbind(pheno,
                            info.geo)
    }
    return(se)

}


# internal function used by extract_geo_characteristics
#
# extracts name and information from single column
# converts to numeric if possible
#
#' @importFrom Hmisc all.is.numeric
#' @keywords internal

convert <- function(col.char) {
    res = vapply(col.char,
                  FUN = function(x) {
                      unlist(strsplit(x, ": "))[2]},
                  FUN.VALUE = character(1))

    res[which(res == "NA")] = NA
    res = all.is.numeric(res,
                         what = "vector",
                         extras = NA)
    name = unlist(strsplit(col.char[1], ": "))[1]
    df = data.frame(res, stringsAsFactors = FALSE)
    colnames(df) = name
    return(df)
}
