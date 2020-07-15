
#' Sex check
#'
#' Predicts sex using k nearest neighbour classification based on expression
#' of sex specific genes.
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay containing
#' expression data to be used for sex check.
#' @param info.sex.genes data.frame. Information about sex chromosomal genes
#' (see \code{info.sex.genes}).
#' @param gene.column Character. Column in rowData containing gene symbols
#' (default: "symbol") or "rownames" if rownames should be used.
#' @param sex.column Character. Column in colData containing sex information.
#' @param col Character vector. Colours to be used for coloring sex in MDS plot
#' (default: red and blue).
#' @param genes Character vector. Symbols of sex specific genes to be used in
#' classification (default: genes collected from the literature, see details).
#' @param title Character. Title of the plot. If NULL title will be set to
#' "MDS (sex specific genes)".
#' @param ellipse Logical. Should ellipses around points be drawn? (default:
#' TRUE).
#' @param ellipse.type Character. Type of ellipse as given in
#' \code{\link[ggpubr]{ggscatter}} (default: "norm").
#'
#' @return List with the following components:
#' \itemize{
#' \item info: data.frame with information about samples with wrong sex or NULL
#' \item plot: MDS plot as returned by \code{plot_mds_pca_2d}
#' }
#'
#' @importFrom caret knn3
#' @importFrom S4Vectors unstrsplit
#' @importFrom stats predict
#' @export
#'
#' @examples
#' data("se.gene")
#' data("info.sex.genes")
#'
#' check_sex(se = se.gene,
#'           info.sex.genes = info.sex.genes,
#'           gene.column = "Gene.symbol",
#'           sex.column = "Sex")

check_sex <- function(se,
                      assay = 1,
                      gene.column = "symbol",
                      info.sex.genes,
                      sex.column = "sex",
                      col = c("red", "blue"),
                      title = NULL,
                      ellipse = TRUE,
                      ellipse.type = "norm") {

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
    if (gene.column == "rownames") {
        anno$gene.to.use = rownames(anno)
    } else {
        if (!(gene.column %in% colnames(anno))) {
            stop(paste("column", gene.column, "not found in rowData!"))
        }
        if (is.factor(anno[, gene.column])) {
            anno$gene.to.use = as.character(anno[, gene.column])
        } else {
            anno$gene.to.use = anno[, gene.column]
        }
    }

    # does not work with SimpleCharacterList (e.g. used in ERP009768)
    anno$gene.char = unstrsplit(#as.character(anno[, gene.column]),
        anno$gene.to.use,
        sep = "|")

    ## identify column of info.sex.genes to use
    ## (maximum number of matches with gene.column)
    cols = c("ensembl_gene_id",
             "hgnc_symbol",
             "entrezgene_id")
    match.l = apply(info.sex.genes[, cols], 2,
                      FUN = function(x) {
                        intersect(anno$gene.char, x)
                      })
    if (length(match.l) == 0) {
        stop("none of the genes in se overlap with genes in info.sex.genes!")
    }
    no.match = vapply(match.l,
                      FUN = length,
                      FUN.VALUE = numeric(1))
    ind.use = which.max(no.match)

    genes.use = rownames(anno)[which(anno$gene.char %in% match.l[[ind.use]])]

    if (length(genes.use) < 2) {
        stop("at least 2 genes are needed to predict sex!")
    }

    se.use = se[genes.use, ]
    se.use = remove_genes(se = se.use,
                          assay = assay,
                          method = "missing",
                          freq = 0)

    res.mds = calculate_mds_pca(se = se.use,
                                assay = assay,
                                method = "mds",
                                dist = "euclidean")
    if (is.null(title)) {
        title = "MDS (sex specific genes)"
    }
    g = plot_mds_pca_2d(res = res.mds,
                        se = se.use,
                        dim = c(1, 2),
                        var.color = "sex",
                        title = title,
                        palette = col,
                        ellipse = ellipse,
                        ellipse.type = ellipse.type)$plot
    #        print(g)

    ## classification using knn
    if (is.character(assay) && !(assay %in% names(assays(se.use)))) {
        stop(paste("assay", assay, "not found!"))
    }
    expr = assays(se.use)[[assay]]

    res.knn = knn3(x = t(expr),
                   y = factor(se.use$sex),
                   k = 5)
    sex.pred = predict(res.knn,
                       t(expr),
                       type = "class")

    ind.wrong = which(sex.pred != se.use$sex)
    if (length(ind.wrong) > 0) {
        info.wrong = data.frame(id = colnames(se)[ind.wrong],
                                sex.pheno = se$sex[ind.wrong],
                                sex.predicted = sex.pred[ind.wrong],
                                stringsAsFactors = FALSE)
    } else {
        info.wrong = NULL
    }
    return(list(info = info.wrong,
                plot = g))

}
