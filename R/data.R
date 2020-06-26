#' Example data set (microarray, probeset level)
#'
#' Microarray data from the GEO data set GSE6710 on the probeset
#' level restricted to the first five samples.
#'
#' @format A
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} with
#' 22283 probesets and 5 samples.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6710}
"se.probeset"

#' Example data set (microarray, gene level)
#'
#' Microarray data from the GEO data set GSE6710 on the gene
#' level.
#'
#' @format A
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} with
#' 12502 genes and 26 samples.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6710}
"se.gene"

#' Genes on the X and Y chromosome
#'
#' Information about human genes on the sex chromosomes (Ensembl release 100
#' GRCh37.p13)
#'
#' @format A data.frame with the following columns:
#' \itemize{
#' \item ensembl_gene_id: Ensembl gene identifier
#' \item chromosome_name: chromosome ('X' or 'Y')
#' \item gene_biotype: gene type (e.g. 'protein_coding')
#' \item hgnc_symbol: gene symbol (HGNC)
#' \item entrezgene_id: NCBI gene identifier
#' }
"info.sex.genes"
