
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
#'
#' @examples
#' se = get_geo_data(accession = "GSE6710")

get_geo_data <- function(accession,
                         temp.dir = tempdir(),
                         platform = NULL,
                         detection.pvalue = FALSE,
                         scan.date = FALSE) {

    #require("GEOquery")
    #require("org.Hs.eg.db")

    ## check that accession is a Series record
    if (!grepl("^GSE", accession)) {
        stop("accession needs to be a Series record (GSExxx)!")
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
            info.platform = vapply(X = temp.l,
                                   FUN = annotation,
                                   FUN.VALUE = character(1))
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
    if (nrow(se) == 0) {
        stop("no expression data found!")
    }

    ## perform log transformation if not already performed
    expr = assays(se)[[1]]
    iqr = IQR(as.numeric(expr),
              na.rm = TRUE)
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

    ## remove .ch1* from colnames
    colnames(colData(se)) = gsub("\\.ch1[0-9_.]*$", "",
                                 colnames(colData(se)))

    return(se)

}


#' Extract scan date
#'
#' Extracts date of scanning from original array files (e.g. CEL for Affymetrix
#' and idat for Illumina arrays).
#'
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param col.sample.id Character or integer vector. Column in colData() with
#' sample information.
#' @param temp.dir Character. Destination directory to downloaded array files.
#' @param tryFormats Character vector. Format strings for date to try.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with scan date added as column scan.date in colData.
#'
#' @import GEOquery
#' @importFrom illuminaio readIDAT
#'
#' @export
#'
#' @examples
#' data("se.probeset")
#'
#' se.probeset = extract_scan_date(se = se.probeset)

extract_scan_date <- function(se,
                              col.sample.id = "geo_accession",
                              temp.dir = tempdir(),
                              tryFormats = c("%Y-%m-%d",
                                             "%Y/%m/%d",
                                             "%m/%d/%Y",
                                             "%m/%d/%y",
                                             "%m-%d-%Y",
                                             "%d-%b-%Y")) {

    if (!(col.sample.id %in% colnames(colData(se)))) {
        stop(paste("column", col.sample.id, "not available in colData!"))
    }
    accession.samples = as.character(colData(se)[, col.sample.id])
    if (!all(grepl("^GSM", accession.samples))) {
        stop(paste("column", col.sample.id,
                   "must contain valid GEO sample ids (GSMxxx)!"))
    }

    files.available = vapply(
        X = accession.samples,
        FUN = function(x) {
            temp = getGEOSuppFiles(GEO = x,
                                   makeDirectory = FALSE,
                                   fetch_files = FALSE)
            ifelse(is.null(temp), FALSE, TRUE)
        },
        FUN.VALUE = logical(1))

    if (all(!files.available)) {
        warning("no raw files with scan date available!")
    } else {

        ## create temp.dir if it does not exist
        if (!file.exists(temp.dir)) {
            dir.create(temp.dir,
                       recursive = TRUE)
        }

        scan.date.all = NULL
        for (i in seq_len(length(accession.samples))) {

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

            ### example Illumina file
            #array.file = system.file("extdata", "idat",
            #                         "4343238080_A_Grn.idat",
            #                         package = "IlluminaDataTestFiles")

            if (length(array.file) == 0) {
                stop(paste("no array file for sample", accession.samples[i],
                           "found!"))
            }
            if (length(array.file) > 1) {
                warning(paste("more than one array file for sample",
                              accession.samples[i],
                              "found, thus using first file\n"))
                array.file = array.file[1]
            }

            if (grepl("cel", array.file, ignore.case = TRUE)) {
                scan.date = get.celfile.date(file = array.file)
            } else if (grepl("idat", array.file, ignore.case = TRUE)) {
                idat = readIDAT(file = array.file)
                ind = which(idat$RunInfo[, "Name"] == "Scan")[1]
                if (length(ind) == 0) {
                    warning(paste("no scan date for file", array.file,
                                  "found!"))
                    scan.date = NA
                } else {
                    scan.date = idat$RunInfo[ind, "Date"]
                    scan.date = unlist(strsplit(scan.date, " "))[1]
                    #scan.date = as.Date(scan.date, format = "%m/%d/%Y")
                }
            } else {
                ## Agilent file (based on GSE32062 with platform GPL6480)
                lines = readLines(array.file, n = 3)

                ## identify index of date
                lines.2.v = unlist(strsplit(lines[2], "\t"))
                ind.scan = grep("Scan_Date", lines.2.v)
                if (length(ind.scan) == 0) {
                    warning(paste("no scan date for file", array.file,
                                  "found!"))
                    scan.date = NA
                } else {
                    if (length(ind.scan) > 1) {
                        warning(paste("more than one scan date for file",
                                      array.file,
                                      "found, use first one!"))
                        ind.scan = ind.scan[1]
                    }
                    lines.3.v = unlist(strsplit(lines[3], "\t",
                                                useBytes = TRUE))
                    scan.date = unlist(strsplit(lines.3.v[ind.scan], " "))[1]
                }
            }
            scan.date.all[i] = scan.date
        }

        scan.date.formats = lapply(tryFormats, function(x) {
            as.Date(scan.date.all,
                    tryFormats = x,
                    optional = TRUE)
        })
        no.na = vapply(X = scan.date.formats,
                       FUN = function(x) {sum(is.na(x))},
                       FUN.VALUE = integer(1))
        ind.min = which.min(no.na)
        if (no.na[ind.min] > sum(is.na(scan.date.all))) {
            warning("some dates could not be correctly converted!")
        }

        se$scan.date = scan.date.formats[[ind.min]]
    }
    return(se)

}


# internal function used by extract_scan_date()
#
# extracts date from Affymetrix CEL files
# based on code from the affio function get.celfile.dates()
# but modified since sometimes tm$ScanDate is empty
#
#' @importFrom affyio read.celfile.header
#'
#' @keywords internal

get.celfile.date <- function(file) {

    tmp = read.celfile.header(file, info = "full")
    if (length(tmp$ScanDate) == 0) {
        date = NA
    } else {
        date = strsplit(tmp$ScanDate, "T| ")[[1]][1]
        date = as.character(as.Date(date,
                                    tryFormats = c("%Y-%m-%d",
                                                   "%m/%d/%y",
                                                   "%Y/%m/%d"),
                                    optional = TRUE))
    }
    return(date)
}


#' Extract detection P value
#'
#' Extracts detection P value information from each sample file.
#'
#' @param accession Character. Identifier of GEO data set (GSEXXX).
#' @param se
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with detection P value added as assay with the name detection.pvalue.
#'
#' @export
#'
#' @examples
#' data("se.probeset")
#'
#' se.probeset = extract_detection_pvalue(se = se.probeset,
#'                                        accession = "GSE6710")

extract_detection_pvalue <- function(accession, se) {

    gse = getGEO(GEO = accession,
                 GSEMatrix = FALSE,
                 getGPL = FALSE,
                 parseCharacteristics = FALSE)
    gsmlist = GSMList(gse)[se$geo_accession]
    info.col = Columns(gsmlist[[1]])

    ## identify column with detection P value
    ind.det.pval = unique(as.numeric(apply(info.col[1], 2, function(x) {
        grep("detect|p-value|pvalue|pval", x, ignore.case = TRUE)})))

    if (length(ind.det.pval) == 0) {
        warning("no information about detection P value available!")
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
                probes.both = intersect(rownames(se), rownames(tab))
                temp = tab[probes.both, ind.det.pval]
                probes.miss = setdiff(rownames(se), rownames(tab))
                if (length(probes.miss) > 0) {
                    temp = c(temp, rep(NA, length(probes.miss)))
                }
                names(temp) = c(probes.both, probes.miss)
                return(temp[rownames(se)])
            }))
            rownames(pval) = rownames(se)

            assays.l = c(assays(se),
                         list(detection.pvalue = pval))
            assays(se) = assays.l
        }
    }
    return(se)
}

