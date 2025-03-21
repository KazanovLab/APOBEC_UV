library(dplyr)
setwd("C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/newfiles")
files <- list.files(pattern = "sample\\d+_with_uv_and_fraction\\.txt$")
all_data <- list()
for (file in files) {
  sample_number <- sub("sample(\\d+)_with_uv_and_fraction\\.txt", "\\1", file)
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  clustered_df <- df %>%
    filter(CLUSTER != "NonClustered", CLUSTER_TYPE != "NONE") %>%
    group_by(CLUSTER) %>%
    summarize(
      Sample = paste0("sample", sample_number),
      Chromosome = first(CHROM),
      Start_Pos = min(POS),
      Num_Mutations = n(),
      UV_Fraction = mean(UV_Fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Store result if not empty
  if (nrow(clustered_df) > 0) {
    all_data[[sample_number]] <- clustered_df
  }
}

# Combine all results if there is data
if (length(all_data) > 0) {
  final_data <- bind_rows(all_data)
  write.table(final_data, "merged_clustered_mutations.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  print("Merged file created successfully!")
} else {
  print("No clustered mutations found in any file.")
}
print(all_data)
View(all_data) #list


library(ggplot2)                                                                                                                                               
library(dplyr) #if you want to continue from this line 
setwd("C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/newfiles")
data <- read.table("merged_clustered_mutations.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert Sample to a factor for better plotting
data$Sample <- as.factor(data$Sample)
data$Sample <- factor(data$Sample, levels = sort(unique(data$Sample), decreasing = FALSE))

# Plot
p1 <- ggplot(data, aes(x = as.factor(Num_Mutations), y = UV_Fraction)) +
  geom_jitter(size=1.5, alpha=0.7,  color = "purple") +  # Dot plot with jitter
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 1, size = 1.5, color = "black") +  # Mean bars
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.6 , color = "black") +  # Error bars (SD)
  labs(title = "UV Fraction vs. Cluster Size (All Samples)", 
       x = "Cluster Size", y = "UV Fraction") +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14))

ggsave("UV_fraction_vs_cluster_size_all_samples.png", p1, width = 7, height = 5)

