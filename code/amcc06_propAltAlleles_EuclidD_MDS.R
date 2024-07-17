# 2021.10.7
# cmdscale: Classical (Metric) Multidimensional Scaling
# Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis (Gower, 1966).
# Multidimensional scaling takes a set of dissimilarities and returns a set of points such that the distances between the points are approximately equal to the dissimilarities. (It is a major part of what ecologists call ‘ordination’.)
# A set of Euclidean distances on n points can be represented exactly in at most n−1 dimensions. 

setwd("/Users/dz359/PycharmProjects/BI/alfalfa_Brian_Heathcliffe_Longxi_16plates/data/")

#args <- commandArgs(trailingOnly = TRUE)
propAltGenosFile = "DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted_called_propAlt.txt"
propAltGenos <- read.delim(propAltGenosFile, row.names=1,header=TRUE)
propAltGenosMat=t(data.matrix(propAltGenos))
euclidD <- dist(propAltGenosMat)
fit <- cmdscale(euclidD,eig=TRUE, k=4)

outp = "DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted_called_propAlt_EuclidD.txt"
write.table(fit$points, outp, sep="\t")



#=========================================================
# Modify the euclidD output generated above and Generate MDS plots
#=========================================================
library(ggplot2)
fit <- read.csv("DAl22-7011_Allele_match_counts_collapsed_combined_modi_Heathcliffe_sorted_called_propAlt_EuclidD_modi_passport.csv", header=TRUE)

head(fit)

pdf(file = "20220812_Alfalfa_Riday_mds_plots.pdf", title = "Alfalfa")
ggplot(fit, aes(x=V1, y=V2, color=Population, shape=DNA_or_leaf)) +
   geom_point()


dev.off()

dim(fit)
outliers <- subset(fit, fit$V1>7)
View(outliers)
write.table(outliers, "20220525-DAl21-6679_EuclidD-outliers.csv", sep=",")