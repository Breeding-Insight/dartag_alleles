library(ggplot2)

# Read and preprocess
data <- read.csv("/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v050_allele_stats.csv")
data$Position <- as.numeric(data$Position)
data$Mb_Bin <- as.numeric(data$Mb_Bin)

# 1. Bar plot: Distribution of number of loci per chromosome
# With the number of rows (loci)grouped by chromosome
ggplot(data, aes(x = Chromosome)) +
  geom_bar(fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Number of Loci per Chromosome",
       x = "Chromosome",
       y = "Loci Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


# 2. Boxplot: Distribution of allele counts across chromosomes
# Provides insight into variability in the number of alleles per chromesome
ggplot(data, aes(x = Chromosome, y = Number_Alleles, fill = Chromosome)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Alleles per Chromosome",
       x = "Chromosome",
       y = "Number of Alleles") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 3. Scatter plot: Locus position vs. Alleles
# visualize how the number of alleles varies along the chromosome (by genomic position)
ggplot(data, aes(x = Position, y = Number_Alleles, color = Chromosome)) +
  geom_point(alpha = 0.7) +  # Scatter plot with transparency
  theme_minimal() +
  labs(title = "Number of Alleles Across Genomic Positions",
       x = "Position (bp)",
       y = "Number of Alleles") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Chromosome, scales = "free_x")  # Create one plot per chromosome


# 4. Density plot: Distribution of allele counts
ggplot(data, aes(x = Number_Alleles, fill = Chromosome)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density Plot of Allele Counts by Chromosome",
       x = "Number of Alleles",
       y = "Density")


ggplot(data, aes(x = Mb_Bin, y = Chromosome, fill = Number_Alleles)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Heatmap of Allele Count Across Chromosomes",
       x = "Mb Bin",
       y = "Chromosome",
       fill = "Number of Alleles")
