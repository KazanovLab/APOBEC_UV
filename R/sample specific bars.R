library(ggplot2)
library(dplyr)
library(gtools)  # Load gtools for correct sorting

setwd("C:/Users/bkurt/Desktop/apobec_data/newfiles")

files <- list.files(pattern = "sample\\d+_with_apobec_and_fraction\\.txt$")

# Use mixedsort() to sort numerically
sorted_files <- mixedsort(files)

results <- list()

for (file in sorted_files) {
  sample_number <- as.numeric(sub("sample(\\d+)_with_apobec_and_fraction\\.txt", "\\1", file))  
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
  summarize(APOBEC_Fraction = sum(`APOBEC` == 1) / n(), .groups = "drop")

sfractions <- as.data.frame(sfractions)

sample1_apobec_mutations <- combined_df %>%
  filter(Sample == 1, `APOBEC` == 1) %>%
  nrow()
sample1_apobec_mutations#correct

pAPO <- ggplot(data = sfractions) +
  aes(x = reorder(Sample, APOBEC_Fraction), y = APOBEC_Fraction)+
  geom_bar(stat = "identity", fill="palegreen3", width = 0.7)+
  geom_text(aes(label = table(combined_df$Sample)), vjust = -0.5) +  # Add text labels on top of bars
  theme_bw()+
  labs(title = "Sample-Specific Fraction of APOBEC-induced Mutations", x = "Sample", y = "The Fraction of APOBEC-induced mutations")+
  theme(plot.title = element_text(hjust = 0.5))
pAPO
ggsave("Sample-specific APOBEC-induced Mutation Fraction.png", pAPO, width = 9, height = 8)

