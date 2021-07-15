
# QCnormSE 0.99.4.9000

## bug fix

check_sex()
- set k to odd number <= minimal group size

extract_scan_date()
- use correct sample file

# QCnormSE 0.99.4

## modifications

remove_genes()  
- additional arguments in edgeR method

## bug fix

remove_genes()  
- correct passing of freq argument in edgeR method

# QCnormSE 0.99.3

## modifications

remove_samples()
- if freq = 1 remove all samples for which criterion is fullfilled for all 
genes

detect_duplicated_samples()  
- add parameters to optionally remove genes with missing values and to specify 
method to deal with missing values in cor and covar functions

# QCnormSE 0.99.2

## bug fixes

combine_se_objects()  
- set metadata of combined SE object to empty list 

extract_scan_date()  
- create temp.dir if it does not exist

get_outliers_mds_pca_2d() (internal function)  
- convert from list to vector if duplicate values occur in outlier detection

plot_heatmap()  
- include lower limits of intervals when defining color categories of P-values 

remove_genes()  
- if freq = 1 remove all genes for which criterion is fullfilled for all 
samples

## modifications

aggregate_by_new_id() 
- keep all assays  
- remove method "sum" since not all assays can be kept if expression values of
several probes are summed up

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
