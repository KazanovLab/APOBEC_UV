#TOTAL NUMBER OF CLUSTERS(ALL SAMPLES) PER CHROMOSOME 
library(ggplot2)

sample1 <- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_0c7aca3f-e006-4de3-afc2-20b4f727d4fd.txt")
sample2<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_2b142863-b963-4cc9-8f8f-c72503c93390.txt")
sample3<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_7d2a22eb-7344-4cba-ad7d-94c3f9ef3d7c.txt")
sample4<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_8c619cbc-9e91-4716-9711-5236e55d8f46.txt")
sample5<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_45a7949d-e63f-4956-866c-df51257032de.txt")
sample6<-read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_448fe471-3f4e-4dc8-a4e0-6f147dc93abe.txt")
sample7<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_94108975-b7a0-40ba-ad39-e44cc62e8cc0.txt")
sample8<- read.table("C:/Users/bkurt/Desktop/trialapo1/merge_results/SigProfilerClusters_cda1a403-16b6-487c-a82a-c377d1d0f89d.txt")
sample9 <- read.table("C:/Users/bkurt/Desktop/trialapo2/merge_results/SigProfilerClusters_abd2d959-d5ed-4eb3-9759-67eb1aa23325.txt")
sample10<- read.table("C:/Users/bkurt/Desktop/trialapo2/merge_results/SigProfilerClusters_acc629cb-ad03-4cec-9b21-922e4932ef3e.txt")
sample11<- read.table("C:/Users/bkurt/Desktop/trialapo2/merge_results/SigProfilerClusters_d4615ca0-b5c7-4a5c-8593-bd50034a78ae.txt")
sample12<- read.table("C:/Users/bkurt/Desktop/trialapo2/merge_results/SigProfilerClusters_df8a913c-5160-4fc5-950d-7c890e24e820.txt")

samplesc <- list(sample1, sample2, sample3, sample4, sample5, sample6,
                 sample7, sample8, sample9, sample10, sample11, sample12)

totnummut<- c("17400", "27489", "10935", "10312", "25277", "36959",
              "65247","41148", "15372", "15783", "15674", "37962")
totnummut <- as.numeric(totnummut)
clustered_s <- list()          # To hold raw cluster data for each sample
results <- list()              # To hold normalized cluster data for each sample
total_clusters <- numeric()    # To store total clusters for all samples


for (i in seq_along(samplesc)) {
  sample_data <- samplesc[[i]]
  sample_data <- sample_data[-1,]
  clustered_data <- sample_data[sample_data$V9 != "NONE", ]
  chromosome_cluster_counts <- aggregate(
    V9 ~ V1,  # V9 contains cluster information, V1 is the chromosome
    data = clustered_data,
    FUN = length
  )
  colnames(chromosome_cluster_counts) <- c("Chromosome", "Cluster_Number")
  chromosome_cluster_counts$Cluster_Density <- (chromosome_cluster_counts$Cluster_Number / totnummut[i]) * 10000
  results[[i]] <- chromosome_cluster_counts
  # Sum normalized counts to get the total for the sample
  total_clusters[i] <- sum(chromosome_cluster_counts$Normalized_Counts)
}

# Set up the PNG device with 600 DPI
png("barplots_600dpi.png", width = 12, height = 8, units = "in", res = 600)  # Adjust width and height as needed

# Arrange plots in a grid
par(mfrow = c(3, 4))  # 3x4 grid

# Iterate over results and generate barplots
chromosome_order <- c(as.character(1:22), "X", "Y")
adjusted_results <- list()

for (i in seq_along(results)) {
  sample_data <- results[[i]]
  sample_data$Chromosome <- factor(sample_data$Chromosome, levels = chromosome_order)
  
  # Add missing chromosomes with Cluster_Number and Cluster_Density set to 0
  missing_chromosomes <- setdiff(chromosome_order, sample_data$Chromosome)
  if (length(missing_chromosomes) > 0) {
    missing_data <- data.frame(
      Chromosome = missing_chromosomes,
      Cluster_Number = 0,
      Cluster_Density = 0
    )
    sample_data <- rbind(sample_data, missing_data)
  }
  
  # Sort data by chromosome order
  sample_data <- sample_data[order(factor(sample_data$Chromosome, levels = chromosome_order)), ]
  adjusted_results[[i]] <- sample_data
  
  # Generate the barplot
  barplot(
    sample_data$Cluster_Density,
    names.arg = sample_data$Chromosome,
    main = paste("Sample", i),
    xlab = "Chromosome",
    ylab = "Cluster Density",
    col = adjustcolor("#4169E1", alpha.f = 0.3),  # Adjust alpha transparency
    las = 2,  # Rotate x-axis labels for readability
    cex.names = 0.8  # Adjust label size
  )
}

# Close the PNG device
dev.off()




# Adjusted results list now includes missing chromosomes set to 0


sample_names <- paste0("Sample", seq_along(results))
combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  cbind(Samples = sample_names[i], results[[i]])
}))
library(writexl)
write_xlsx(combined_results, path = "C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/plots and data-apobec/samspecchrom.xlsx")

# Define chromosome order and sample names
chromosome_order <- c(as.character(1:22), "X", "Y")
sample_names <- paste0("Sample", seq_along(results))

# Reorder and fill data for all samples
for (i in seq_along(results)) {
  sample_data <- results[[i]]
  
  # Convert Chromosome column to factor with the desired order
  sample_data$Chromosome <- factor(
    sample_data$Chromosome,
    levels = chromosome_order
  )
  
  # Sort data by the custom chromosome order
  sample_data <- sample_data[order(sample_data$Chromosome), ]
  
  # Update the results with the sorted data
  results[[i]] <- sample_data
}

# Combine all samples into one data frame
df <- do.call(rbind, lapply(seq_along(results), function(i) {
  sample_data <- as.data.frame(results[[i]])  
  sample_data$Sample <- sample_names[i]       
  return(sample_data)
}))

# Fill NAs with 0 for Cluster_Density
library(dplyr)
library(tidyr)
df <- df %>%
  complete(Sample, Chromosome, fill = list(Cluster_Density = 0))  # Ensures all combinations are filled

# Optionally reorder samples by total cluster density
df$Sample <- reorder(df$Sample, df$Cluster_Density, FUN = sum)

# Create the heatmap
library(ggplot2)
heatmap <- ggplot(df, aes(x = Chromosome, y = Sample, fill = Cluster_Density)) +
  geom_tile(color = "white") +  # Add white grid lines
  scale_fill_gradientn(
    colors = c("lightyellow", "orange", "red", "darkred"),  # Color palette
    limits = c(0, 200),  # Adjust scale limits based on your data
    name = "Cluster Density"
  ) +
  labs(
    title = "Heatmap of Cluster Density by Chromosomes Across Samples",
    x = "Chromosome",
    y = "Sample"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),  # Center and size the title
    axis.text.x = element_text(angle = 45, hjust = 1, size =20),
    axis.text.y = element_text(size =20)
)

# Print the heatmap
print(heatmap)

# Save the heatmap with 600 DPI
ggsave(
  filename = "apoh2.png",
  plot = heatmap,
  width = 10, height = 8, units = "in", dpi = 600
)


