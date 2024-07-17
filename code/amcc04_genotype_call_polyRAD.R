library(VariantAnnotation)
library(polyRAD)
library(qqman)
library(dplyr)


# Check is "plyr" is loaded or not
if(any(grepl("package:dplyr", search()))) detach("package:dplyr") else message("dplyr not loaded")


packageVersion("dplyr")

getwd()
setwd("/Users/dz359/PycharmProjects/BI/alfalfa_Brian_Heathcliffe_Longxi_16plates/data/")

data <- readVcf(file = "DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted.vcf")

# Compress and index VCF
Rsamtools::bgzip("DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted.vcf")
Rsamtools::indexTabix("DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted.vcf.bgz", format = "vcf")

# Index reference genome
Rsamtools::indexFa("/Users/dz359/PycharmProjects/BI/demo_DArTag_validation_QC/data/ref_genome/XinJiangDaYe_set1_monoploid.fa")

# Importing the data, with very loose filtering

?VCF2RADdata
mydata <- VCF2RADdata("DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted.vcf.bgz",
                      refgenome = "/Users/dz359/PycharmProjects/BI/demo_DArTag_validation_QC/data/ref_genome/XinJiangDaYe_set1_monoploid.fa",
                      phaseSNPs = FALSE,
                      min.ind.with.reads = 0,
                      min.ind.with.minor.allele = 0,
                      possiblePloidies = list(4))
# Merging rare haplotypes...
# 2869 markers retained out of 3000 originally.

# heterozygosity-based quality statistic
hh <- HindHe(mydata)

rownames(hh)

mydata <- SubsetByTaxon(mydata, rownames(hh))

mydataPopStruct <- IteratePopStruct(mydata, overdispersion = 6)

genomat <- GetWeightedMeanGenotypes(mydataPopStruct)

str(genomat)

RADdata2VCF(mydataPopStruct, file = "DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted_called.vcf")
