
# QCnormSE 0.99.1

## bug fixes

check_sex()
- sex.column can now also be a factor
- error if none of the genes in se overlap with genes in info.sex.genes

define_batches()
- ensure that scan.date is of class Date

extract_scan_date()
- col.sample.id can now also be a factor
- additional possible date format for CEL files
- different positions of scan date information possible in Agilent raw files

detect_duplicate_samples()
- more than one duplicate shown as vertical red line

plot_mds_pca_2d()
- remove unused levels of var.color and var.shape

## new features

new function combine_se_objects()

## modifications

check_batch_effects()
- additional heatmap showing absolute adjusted r^2 values

check_sex()
- rownames can be used as gene identifier
- add ellipses to MDS plot

detect_duplicate_samples()
- improve performance by using default values of function stats::cor() and 
provide access to fast implementation in R package coop
- add option to set title of plot

extract_scan_date()
add option to specify column with sample information


# QCnormSE 0.99.0

- initial commit
