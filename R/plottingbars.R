library(ggplot2)
library(dplyr)
library(gtools)  # Load gtools for correct sorting

setwd("C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/newfiles")

files <- list.files(pattern = "sample\\d+_with_uv_and_fraction\\.txt$")

# Use mixedsort() to sort numerically
sorted_files <- mixedsort(files)

results <- list()

for (file in sorted_files) {
  sample_number <- as.numeric(sub("sample(\\d+)_with_uv_and_fraction\\.txt", "\\1", file))  
  print(paste("Processing:", file, "Sample:", sample_number))  # Debugging output
  
  results[[as.character(sample_number)]] <- read.table(file, header = TRUE, sep = "\t", 
                                                       stringsAsFactors = FALSE, check.names = FALSE)
}

sorted_sample_numbers <- sort(as.numeric(names(results)))

combined_df <- bind_rows(lapply(sorted_sample_numbers, function(sample_num) {
  results[[as.character(sample_num)]] %>% 
    mutate(Sample = sample_num)  # Add Sample column
}))
table(combined_df$Sample)

sfractions <- combined_df %>%
  group_by(Sample) %>%
  summarize(UV_Fraction = sum(`UV-Induced` == 1) / n(), .groups = "drop")
#checking
sample5_uv_mutations <- combined_df %>%
  filter(Sample == 5, `UV-Induced` == 1) %>%
  nrow()
sample5_uv_mutations#correct

sfractions <- as.data.frame(sfractions)

puv <- ggplot(data = sfractions) +
  aes(x = reorder(Sample, UV_Fraction), y = UV_Fraction)+
  geom_bar(stat = "identity", fill="deepskyblue3", width = 0.7)+
  geom_text(aes(label = table(combined_df$Sample)), vjust = -0.5) +  # Add text labels on top of bars
  theme_bw()+
  labs(title = "Sample-Specific Fraction of UV-induced Mutations ", x = "Sample", y = "The Fraction of UV-induced Mutations")+
  theme(plot.title = element_text(hjust = 0.5))
puv
ggsave("Sample-specific UV_fraction.png", puv, width = 9, height = 8)

