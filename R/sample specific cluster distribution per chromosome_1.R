#TOTAL NUMBER OF CLUSTERS(ALL SAMPLES) PER CHROMOSOME 
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
chrom_length<- c("248,956,422", "133,797,422", "135,086,622", "133,275,309", "114,364,328",
                 "107,043,718", "101,991,189", "90,338,345", "83,257,441", "80,373,285", "58,617,616",
                 "242,193,529", "64,444,167", "46,709,983", "50,818,468", "50,818,468", "190,214,555",
                 "181,538,259", "170,805,979", "159,345,973", "145,138,636", "138,394,717", "156,040,895","57,227,415")
chrom_length <- as.numeric(gsub(",", "", chrom_length))
results <- list()
total_clusters <- numeric()   
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
  chromosome_cluster_counts$Normalized_Counts <- (chromosome_cluster_counts$Cluster_Number / chrom_length[i]) * 10000000
  # Store results for the current sample
  results[[i]] <- chromosome_cluster_counts
  # Sum normalized counts to get the total for the sample
  total_clusters[i] <- sum(chromosome_cluster_counts$Normalized_Counts)
}

par(mfrow = c(3, 5))

for (i in seq_along(results)) {
  sample_data <- results[[i]]
  sample_data$Chromosome <- reorder(
    sample_data$Chromosome,
    -sample_data$Normalized_Counts
  )
  barplot(
    sample_data$Normalized_Counts[order(-sample_data$Normalized_Counts)],
    names.arg = sample_data$Chromosome[order(-sample_data$Normalized_Counts)],
    main = paste("Sample", i),
    xlab = "Chromosome",
    ylab = "The Total Number of Clusters",
    col = "skyblue",
    las = 2,  # Rotate x-axis labels for readability
    cex.names = 0.8  # Adjust label size
  )
}

sample_names <- paste0("Sample", seq_along(results))
combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  cbind(Samples = sample_names[i], results[[i]])
}))
library(writexl)
write_xlsx(combined_results, path = "C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/plots and data/samspecchrom.xlsx")

