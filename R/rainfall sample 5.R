somatic.mutations <- read.table(file="C:/Users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/rainfall plot data/SigProfilerClusters_20e02396-e676-412d-9724-44a428919cdb.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
somatic.mutations <- setNames(somatic.mutations, c("chr", "pos", "ıd", "ref", "alt", "qual", "filter", "info", "cluster", "cluster type"))
somatic.mutations$chr <- as.character(somatic.mutations$chr)
head(somatic.mutations)
# Filter mutations for chromosome 11
chr11.mutations <- somatic.mutations[somatic.mutations$chr == "11", ]

# Sort by start position
chr11.mutations <- chr11.mutations[order(chr11.mutations$pos), ]

# Calculate intermutation distances
chr11.mutations$pos <- as.numeric(chr11.mutations$pos)
chr11.mutations$intermutation_distance <- c(NA, diff(chr11.mutations$pos))

# Remove the first row with NA distance (optional for clean plotting)
chr11.mutations <- chr11.mutations[!is.na(chr11.mutations$intermutation_distance), ]

# Take the log10 of the intermutation distances
chr11.mutations$log10_distance <- log10(chr11.mutations$intermutation_distance)

# Plot the rainfall plot
library(ggplot2)
ggplot(chr11.mutations, aes(x = pos, y = log10_distance, color = `cluster type`)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Rainfall Plot for Chromosome 11",
       x = "Genomic Position (bp)",
       y = "Log10(Intermutation Distance)",
       color = "Mutation Type") +
  theme(plot.title = element_text(hjust = 0.5))


library(scales)  # For comma() function

ggplot(chr11.mutations, aes(x = pos, y = log10_distance, color = `cluster type`)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Rainfall Plot for Chromosome 11",
    x = "Genomic Position (bp)",
    y = "Log10(Intermutation Distance)",
    color = "Mutation Type"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma)  # Convert scientific notation to standard numbers

#substitution type a göre deneyelim
# Create a substitution type column
chr11.mutations$substitution <- paste(chr11.mutations$ref, chr11.mutations$alt, sep = ">")

ggplot(chr11.mutations, aes(x = pos, y = log10_distance, color = substitution)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_minimal() +
  theme_bw() +
  labs(
    title = "Rainfall Plot for Chromosome 11 in Sample 5",
    x = "Genomic Position (bp)",
    y = "Log10(Intermutation Distance)",
    color = "Mutation Type"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma)

#finding the most observed substitution type: 
# Count the occurrences of each substitution type
substitution_counts <- table(chr11.mutations$substitution)

# Sort the counts in descending order
sorted_counts <- sort(substitution_counts, decreasing = TRUE)

# Display the sorted counts
print(sorted_counts)

# Find the most observed substitution type
most_observed <- names(sorted_counts)[1]
most_observed_count <- sorted_counts[1]

cat("The most observed substitution type is", most_observed, "with", most_observed_count, "occurrences.\n")



















# Load necessary libraries
library(ggplot2)
library(scales)
library(RColorBrewer)

# Define a more visually appealing color palette
custom_colors <- c(
  "C>A" = "#1b9e77",    # Green
  "C>T" = "#d95f02",    # Orange
  "C>G" = "#7570b3",    # Purple
  "T>A" = "#e7298a",    # Pink
  "T>C" = "#66a61e",    # Lime Green
  "T>G" = "#e6ab02"     # Yellow
)

# Generate the plot
ggplot(chr11.mutations, aes(x = pos, y = log10_distance, color = substitution)) +
  geom_point(size = 2.5, alpha = 0.9) +  # Larger points with better transparency
  theme_minimal() +
  labs(
    title = "Rainfall Plot for Chromosome 11",
    x = "Genomic Position (bp)",
    y = "Log10(Intermutation Distance)",
    color = "Mutation Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Larger title
    axis.title = element_text(size = 14),  # Larger axis titles
    axis.text = element_text(size = 12),  # Larger axis text
    legend.title = element_text(size = 12),  # Larger legend title
    legend.text = element_text(size = 10)  # Larger legend text
  ) +
  scale_x_continuous(labels = comma) +
  scale_color_manual(values = custom_colors)  # Apply the custom color palette
