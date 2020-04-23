
#' Extract GEO information
#'
#' Extracts additional sample information from GEO into separate columns of
#' colData(). Numeric variables are converted from character
#' to numeric and syntactically valid variable names are used for the new
#' columns.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param cols Character or integer. Name or number of column(s) in colData with
#' GEO information (optional, will otherwise be extracted assuming column names
#' contain term 'characteristics').
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with GEO information in additional columns in colData()
#'
#' @importFrom methods is
#' @export

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

    ## extract information
    info.geo = do.call(cbind, apply(info.geo.original, 2, convert))

    ## convert colnames to syntactically valid variable names
    colnames(info.geo) = base::make.names(colnames(info.geo))

    colData(se) = cbind(pheno,
                        info.geo)
    return(se)

}


#' internal function used by extract_geo_characteristics
#'
#' extracts name and information from single column
#' converts to numeric if possible
#'
#' @importFrom Hmisc all.is.numeric
#' @keywords internal

convert <- function(col.char) {
    res = sapply(col.char, function(x) {
        unlist(strsplit(x, ": "))[2]})
    res[which(res == "NA")] = NA
    res = all.is.numeric(res,
                         what = "vector",
                         extras = NA)
    name = unlist(strsplit(col.char[1], ": "))[1]
    df = data.frame(res, stringsAsFactors = FALSE)
    colnames(df) = name
    return(df)
}
