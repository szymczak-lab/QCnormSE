
## todo:
## y chromosomal genes
## or analysis based on X and Y chromosomal genes
## references to selected genes

#' Sex check
#'
#' Predicts sex using k nearest neighbour classification based on expression
#' of sex specific genes.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' expression data to be used for sex check.
#' @param symbol.column Character. Column in rowData containing gene symbols
#' (default: "symbol").
#' @param sex.column Character. Column in colData containing sex information.
#' @param col Character vector. Colours to be used for coloring sex in MDS plot
#' (default: red and blue).
#' @param genes Character vector. Symbols of sex specific genes to be used in
#' classification (default: genes collected from the literature, see details).
#' @param plot Logical. Should MDS plot be generated (default: TRUE).
#'
#' @return data.frame with information about samples with wrong sex or NULL.
#'
#' @importFrom caret knn3
#' @importFrom S4Vectors unstrsplit
#' @importFrom stats predict
#' @export
#'
#' @examples
#' data("se.gene")
#'
#' check_sex(se = se.gene,
#'           symbol.column = "Gene.symbol",
#'           sex.column = "Sex")

check_sex <- function(se,
                      assay = 1,
                      symbol.column = "symbol",
                      sex.column = "sex",
                      col = c("red", "blue"),
                      genes = c("DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "RPS4Y1",
                                "TXLNG2P", "UTY", "XIST"),
                      plot = TRUE) {

    ## extract sex information
    pheno = colData(se)

    if (!(sex.column %in% colnames(pheno))) {
        stop(paste("column", sex.column, "not found in colData!"))
    }
    sex = pheno[, sex.column]
    se$sex = sex

    ## remove samples with missing sex information
    ind.rm = which(is.na(se$sex))
    if (length(ind.rm) > 0) {
        warning(paste("removed", length(ind.rm), "samples with missing",
                      "sex information!"))
        se = se[, -ind.rm]
    }

    ## extract genes
    anno = rowData(se)
    if (!(symbol.column %in% colnames(anno))) {
        stop(paste("column", symbol.column, "not found in rowData!"))
    }
    anno$symbol.char = unstrsplit(as.character(anno[, symbol.column]),
                                  sep = "|")
    genes.use = rownames(anno)[which(anno$symbol.char %in% genes)]

    if (length(genes.use) < 2) {
        stop("at least 2 genes are needed to predict sex!")
    }

    if (plot) {
        res.mds = calculate_mds_pca(se = se[genes.use, ],
                                    assay = assay,
                                    method = "mds",
                                    dist = "euclidean")
        g = plot_mds_pca_2d(res = res.mds,
                            se = se,
                            dim = c(1, 2),
                            var.color = "sex",
                            title = "MDS (sex specific genes)",
                            palette = col)
        print(g)
    }

    ## classification using knn
    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se[genes.use, ])[[assay]]

    res.knn = knn3(x = t(expr),
                   y = factor(se$sex),
                   k = 5)
    sex.pred = predict(res.knn,
                       t(expr),
                       type = "class")

    ind.wrong = which(sex.pred != se$sex)
    if (length(ind.wrong) > 0) {
        info.wrong = data.frame(id = colnames(se)[ind.wrong],
                                sex.pheno = se$sex[ind.wrong],
                                sex.predicted = sex.pred[ind.wrong],
                                stringsAsFactors = FALSE)
    } else {
        info.wrong = NULL
    }
    return(info.wrong)

}
