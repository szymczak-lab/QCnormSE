## code to prepare datasets

#detach("package:QCnormSE", unload = TRUE)
library(QCnormSE)
library(SummarizedExperiment)

## SE microarray probeset level
se.probeset.all = get_geo_data(accession = "GSE6710",
                               scan.date = TRUE,
                               detection.pvalue = TRUE)
se.probeset.all = extract_scan_date(se = se.probeset.all)

## reduce number of probesets
genes = c("DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y",
          "RPS4Y1", "TXLNG2P", "UTY", "XIST")
ind = unlist(sapply(genes, function(g) {
    grep(g, rowData(se.probeset.all)$Gene.symbol)}))
ind.use = union(ind, 1:1988)

se.probeset = se.probeset.all[ind.use, 1:5]

usethis::use_data(se.probeset, overwrite = TRUE)


## SE microarray gene level
se.probeset = se.probeset.all[ind.use, ]

se.gene = aggregate_by_new_id(se = se.probeset,
                              assay = "exprs.log",
                              col.new = "Gene.symbol",
                              sep = "///")

se.gene$group = ifelse(grepl("lesional", se.gene$title),
                       "lesional", "uninvolved")

usethis::use_data(se.gene, overwrite = TRUE)
