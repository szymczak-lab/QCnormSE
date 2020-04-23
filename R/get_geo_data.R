
#' Export gene expression data sets from GEO
#'
#' Downloads normalized gene expression data from Gene Expression Omnibus (GEO)
#' and stores it as a SummarizedExperiment object. Information about detection
#' P values and dates of scanning can also be extracted.
#'
#' @param accession Character. Identifier of GEO data set (GSEXXX).
#' @param temp.dir Character. Destination directory to downloaded array files.
#' @param platform Character. GEO identifier of platform that should be
#' selected. Only relevant, if study profiled samples on different platforms.
#' @param detection.pvalue Logical. Should information about detection P values
#' be extracted (if available)?
#' @param scan.date Logical. Should information about scan dates be extracted
#' (if raw data files are available)?
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object.
#'
#' @import GEOquery
#' @importFrom stats IQR
#' @importFrom BiocGenerics annotation
#' @export

get_geo_data <- function(accession,
                         temp.dir = tempdir(),
                         platform = NULL,
                         detection.pvalue = FALSE,
                         scan.date = FALSE) {

  #require("GEOquery")
  #require("org.Hs.eg.db")

    ## check that accession is a Series record
    if (!grepl("^GSE", accession)) {
        stop("acsession needs to be a Series record (GSExxx)!")
    }

  ## download from GEO
  temp.l = getGEO(GEO = accession,
                  GSEMatrix = TRUE,
                  AnnotGPL = TRUE,
                  destdir = temp.dir)

  ## some studies contain samples from different platforms which are stored as
  ## separate ExpressionSet objects
  if (length(temp.l) > 1) {
        if (is.null(platform)) {
            stop(paste("more than one ExpressionSet returned for accession",
                       accession, "and no platform specified!"))
        } else {
            info.platform = sapply(temp.l, annotation)
            ind = which(info.platform == platform)
            if (length(ind) == 0) {
                stop(paste("platform", platform, "does not match information",
                           "about downloaded platforms:",
                           paste(info.platform, collapse = ", ")))
            } else {
                eset = temp.l[[ind]]
            }
        }

  } else {
      eset = temp.l[[1]]
  }

  ## convert to SummarizedExperiment object
  se = as(eset, "SummarizedExperiment")

  ## perform log transformation if not already performed
  expr = assays(se)[[1]]
  iqr = IQR(as.numeric(expr))
  if (iqr > 100) {
      print("performing log transformation ...")
      se = log_transform(se = se,
                         assay = 1,
                         pseudocount = 1)
  }

  ## extract detection p-values
  if (detection.pvalue) {
    se = extract_detection_pvalue(accession = accession,
                                  se = se)
  }

  if (scan.date) {
    se = extract_scan_date(se = se,
                           temp.dir = temp.dir)

  }
  return(se)

  # prepare.sum.exp(se.gene = se.gene, assay = "expr", col.symbol = "symbol",
  #                 se.file = file.path(se.dir, paste0("se_", accession, ".rds")),
  #                 remove.zeros = FALSE, type.detection = type.detection)
}


#' Extract scan date
#'
#' Extracts date of scanning from original array files (e.g. CEL for Affymetrix
#' and idat for Illumina arrays).
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param temp.dir Character. Destination directory to downloaded array files.
#' @param tryFormats Character vector. Format strings for date to try.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with scan date added as column scan.date in colData.
#'
#' @import GEOquery
#' @importFrom affyio get.celfile.dates
#' @importFrom illuminaio readIDAT
#'
#' @export

extract_scan_date <- function(se,
                              temp.dir = tempdir(),
                              tryFormats = c("%Y-%m-%d",
                                             "%Y/%m/%d",
                                             "%m/%d/%Y",
                                             "%d-%b-%Y")) {

    accession.samples = se$geo_accession

    scan.date.all = NULL
    for (i in 1:length(accession.samples)) {

        array.file = dir(temp.dir,
                         pattern = accession.samples[i],
                         full.names = TRUE)

        ## avoid downloading of existing file
        if (length(array.file) == 0) {
            getGEOSuppFiles(GEO = accession.samples[i],
                            makeDirectory = FALSE,
                            baseDir = temp.dir,
                            fetch_files = TRUE)
            array.file = dir(temp.dir,
                             pattern = accession.samples[i],
                             full.names = TRUE)
        }

        #    ## example Illumina file
        #    array.file = system.file("extdata", "idat", "4343238080_A_Grn.idat",
        #                             package = "IlluminaDataTestFiles")

        if (length(array.file) == 0) {
            stop(paste("no array file for sample", accession.samples[i],
                       "found!"))
        }
        if (length(array.file) > 1) {
            stop(paste("more than one array file for sample",
                       accession.samples[i], "found!"))
        }

        if (grepl("cel", array.file, ignore.case = TRUE)) {
            scan.date = as.character(get.celfile.dates(filenames = array.file))
        } else if (grepl("idat", array.file, ignore.case = TRUE)) {
            idat = readIDAT(file = array.file)
            ind = which(idat$RunInfo[, "Name"] == "Scan")[1]
            if (length(ind) == 0) {
                stop(paste("no scan date for file", array.file, "found!"))
            }
            scan.date = idat$RunInfo[ind, "Date"]
            scan.date = unlist(strsplit(scan.date, " "))[1]
#            scan.date = as.Date(scan.date, format = "%m/%d/%Y")
        } else {
            ## Agilent file (based on GSE32062 with platform GPL6480)
            temp = unlist(strsplit(readLines(array.file, n = 3)[[3]], "\t"))
            scan.date = unlist(strsplit(temp[3], " "))[1]
#            stop(paste("no method available to extract date from filetype",
#                       array.file))
        }
        scan.date.all[i] = scan.date
    }

    se$scan.date = as.Date(scan.date.all,
                          tryFormats = tryFormats)
    return(se)

}


#' Extract detection P value
#'
#' Extracts detection P value information from each sample file.
#'
#' @param accession Character. Identifier of GEO data set (GSEXXX).
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with detection P value added as assay with the name detection.pvalue.
#'
#' @export

extract_detection_pvalue <- function(accession, se) {

    gse = getGEO(GEO = accession,
                 GSEMatrix = FALSE)
    gsmlist = GSMList(gse)
    info.col = Columns(gsmlist[[1]])

    ## identify column with detection P value
    ind.det.pval = unique(as.numeric(apply(info.col, 2, function(x) {
        grep("detect|p-value|pvalue", x, ignore.case = TRUE)})))

    if (length(ind.det.pval) == 0) {
        stop("no information about detection P value available!")
    } else {
        if (length(ind.det.pval) > 1) {
            stop(paste("more than one column with detection P value available:",
                       paste(info.col[ind.det.pval, "Column"],
                             collapse = ", ")))
        } else {
            print(paste0("extract additional column ",
                         info.col[ind.det.pval, "Column"]))

            pval = do.call("cbind", lapply(gsmlist, function(x) {
                tab = Table(x)
                rownames(tab) = tab$ID_REF
                if (!all(rownames(se) %in% rownames(tab))) {
                    stop("Some probes not found in detection P value information!")
                }
                return(tab[rownames(se), ind.det.pval])
            }))
            rownames(pval) = rownames(se)

            assays.l = c(assays(se),
                            list(detection.pvalue = pval))
            assays(se) = assays.l
        }
    }
    return(se)
}

