# Load processed gene expression data
expression_data <- readRDS(
  "data/processed/expression_data.rds"
)

# Load sample information
sample_info <- read.csv(
  "data/processed/sample_info.csv"
)

# Check expression data dimensions
dim(expression_data)

# Check sample information dimensions
dim(sample_info)

# Check expression data structure
str(expression_data)

# View first few expression values
expression_data[1:6, 1:5]

# View first few sample details
head(sample_info)


hist(
  expression_data,
  breaks = 50,
  main = "Distribution of Gene Expression Values",
  xlab = "Gene Expression Values",
  col = "lightblue",
  border = "white"
)

boxplot(
  as.vector(expression_data),
  main = "Boxplot of Gene Expression Values",
  ylab = "Expression Values",
  col = "lightgreen",
  horizontal = TRUE
)

hist(
  log2(expression_data + 1),
  breaks = 50,
  main = "Distribution after Log2 Transformation",
  xlab = "Log2(Expression + 1)",
  col = "lightpink"
)


boxplot(
  log_expression,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Values",
  outline = FALSE,
  las = 2,
  col = "lightblue"
)


sample_medians <- apply(
  log_expression,
  2,
  median
)

summary(sample_medians)


sample_means <- apply(
  log_expression,
  2,
  mean
)

summary(sample_means)



# Save sample-wise statistics

sample_statistics <- data.frame(
  Sample = colnames(log_expression),
  Mean = sample_means,
  Median = sample_medians
)

write.csv(
  sample_statistics,
  "results/sample_statistics.csv",
  row.names = FALSE
)