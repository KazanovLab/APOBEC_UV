library(ggplot2)

sample1 <- read.table("C:/Users/bkurt/Desktop/trial/merge_results/SigProfilerClusters_1ac15380-04a2-42dd-8ade-28556a570e80.txt")
sample2<- read.table("C:/Users/bkurt/Desktop/trial/merge_results/SigProfilerClusters_1cd0acf2-3116-4dfa-a063-0a435b9f6da3.txt")
sample3<- read.table("C:/Users/bkurt/Desktop/trial/merge_results/SigProfilerClusters_1d4a091d-fe65-49c0-8810-5a95243b108a.txt")
sample4<- read.table("C:/Users/bkurt/Desktop/trial2/merge_results/SigProfilerClusters_12f1ae2f-2666-45be-9742-f502d212373d.txt")
sample5<- read.table("C:/Users/bkurt/Desktop/trial2/merge_results/SigProfilerClusters_20e02396-e676-412d-9724-44a428919cdb.txt")
sample6<-read.table("C:/Users/bkurt/Desktop/trial2/merge_results/SigProfilerClusters_22d67778-61fc-4f15-95b8-7e7c6cc7112b.txt")
sample7<- read.table("C:/Users/bkurt/Desktop/trial2/merge_results/SigProfilerClusters_25e20393-752b-4796-9001-0e22ee04c586.txt")
sample8<- read.table("C:/Users/bkurt/Desktop/trial2/merge_results/SigProfilerClusters_9e0009d1-c993-4247-9706-88ee84591dec.txt")
sample9<- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_08b5d0e4-4661-460e-a9f7-f2e687414711.txt")
sample10 <- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_0ab4d782-9a50-48b9-96e4-6ce42b2ea034.txt")
sample11<- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_2e76891c-b620-4cc0-9315-6f1217b09b1e.txt")
sample12<- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_4e8396f7-9506-4401-96b6-bb2e89557d59.txt")
sample13<-  read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_7e22401d-f4cd-44c5-8a01-b08a439e5a31.txt")
sample14<- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_7edc42d3-d08e-4360-a3e1-aeb57cfc6640.txt")
sample15<- read.table("C:/Users/bkurt/Desktop/trial3/merge_results/SigProfilerClusters_7f031d71-3cb7-4744-86bd-a3beecfe166e.txt")

samplesc <- list(sample1, sample2, sample3, sample4, sample5, sample6,
                 sample7, sample8, sample9, sample10, sample11, sample12,
                 sample13, sample14, sample15)
totnummut<- c("8272", "9079", "161472", "97305", "4974", "96022",
              "155414","15352", "80034", "236839", "54048", "28409","52151",
              "8996","2995")
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
  Total_Normalized_Clusters = total_clusters
)


b <- ggplot(data = plot_data) +
  aes(x = Total_Normalized_Clusters, y = reorder(Sample, Total_Normalized_Clusters)) +
  geom_bar(stat = "identity", fill = "darkblue", alpha = 0.6, width = 0.6) +
  geom_text(aes(label = totnummut), hjust = -0.1, size = 4, nudge_y = 0) + # Add text near bars
  labs(title = "Cluster Density Across Samples",
       x = "Cluster Density (per 10K nt)",
       y = "Sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("melanoma.png", plot = b, width = 12, height = 8, dpi = 600)



sample_names <- paste0("Sample", seq_along(results))
combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  cbind(Samples = sample_names[i], results[[i]])
}))
library(writexl)
write_xlsx(combined_results, path = "C:/Users/bkurt/Desktop/cresults.xlsx")
