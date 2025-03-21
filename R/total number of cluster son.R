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
  clustered_s[[i]] <- sample_data[["V9"]][sample_data[["V9"]] != "NONE"]
  cluster_counts <- table(clustered_s[[i]])
  cluster_counts <- cluster_counts[cluster_counts > 1]
  sorted_cluster_counts <- sort(cluster_counts, decreasing = TRUE)
  cluster_summary <- as.data.frame(sorted_cluster_counts)
  colnames(cluster_summary) <- c("Cluster", "Freq")
  cluster_size_summary <- as.data.frame(table(cluster_summary$Freq))
  colnames(cluster_size_summary) <- c("Cluster_Size", "Number_of_Clusters")
  cluster_size_summary$Normalized_Counts <- (cluster_size_summary$Number_of_Clusters / totnummut[i]) * 10000
  results[[i]] <- cluster_size_summary
  total_clusters[i] <- sum(cluster_size_summary$Normalized_Counts)
}

plot_data <- data.frame(
  Sample = paste0("Sample ", seq_along(samplesc)),
  Cluster_Density = total_clusters
)

a <- ggplot(data = plot_data) +
  aes(x = Cluster_Density, y = reorder(Sample, Cluster_Density)) +
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6, fill = "royalblue") +  # Fill the bar with a single color
  geom_text(aes(label = totnummut), hjust = -0.1, size = 4, nudge_y = 0) +  # Add text near bars
  labs(title = "Cluster Density Across Samples",
       x = "Cluster Density (per 10K nt)",
       y = "Sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))  # Center align title
ggsave("apobec.png", plot = a, width = 12, height = 8, dpi = 600)


sample_names <- paste0("Sample", seq_along(results))
combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  cbind(Samples = sample_names[i], results[[i]])
}))

exdata <- data.frame(
  Sample = paste0("Sample ", seq_along(samplesc)),
  Total_Number_of_Clusters = total_number_of_clusters,
  Total_Number_of_Mutations = totnummut,
  Cluster_Density = total_clusters
)

library(writexl)
write_xlsx(exdata, path = "C:/Users/bkurt/Desktop/exdata.xlsx")
