## code to prepare datasets

#detach("package:QCnormSE", unload = TRUE)
library(QCnormSE)
library(SummarizedExperiment)
library(biomaRt)

## SE microarray probeset level
se.probeset.all = get_geo_data(accession = "GSE6710",
                               scan.date = TRUE,
                               detection.pvalue = TRUE)
se.probeset.all = extract_scan_date(se = se.probeset.all)

## modify some of the phenotype information
pheno = colData(se.probeset.all)
pheno$Age.of.patient = as.numeric(gsub(" years", "", pheno$Age.of.patient))
pheno$Duration.of.psoriasis = as.numeric(gsub(" years", "",
                                              pheno$Duration.of.psoriasis))
pheno$Body.surface.area = as.numeric(gsub(" percent", "",
                                          pheno$Body.surface.area))
pheno$Induration = gsub(" \\(classified by the investigators according to the following six grades: clear  -  minimal  -  mild  -  moderate  -  severe  -  very severe\\)",
                        "", pheno$Induration)
pheno$Scaling = gsub(" \\(classified by the investigators according to the following six grades: clear  -  minimal  -  mild  -  moderate  -  severe  -  very severe\\)",
                        "", pheno$Scaling)
colData(se.probeset.all) = pheno

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


## gene annotation

# load mart object from ensembl
# 2020-06-26: release 100
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

#listAttributes(mart)
chrom = c("X", "Y")

# extract gene symbols, Ensembl ids and ENTREZ ids on chromosome X and Y
info.sex.genes = getBM(attributes = c("ensembl_gene_id",
                                      "chromosome_name",
                                      "gene_biotype",
                                      "hgnc_symbol",
                                      "entrezgene_id"),
                       filters = "chromosome_name",
                       values = chrom,
                       mart = mart)
info.sex.genes$hgnc_symbol[which(info.sex.genes$hgnc_symbol == "")] = NA

usethis::use_data(info.sex.genes,
                  overwrite = TRUE)

