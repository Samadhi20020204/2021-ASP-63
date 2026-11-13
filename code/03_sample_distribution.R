
# Day 3: Gene Expression Data Transformation and Sample-wise Analysis

# 1. Log2 transformation
log_expression <- log2(expression_data + 1)

# 2. Summary statistics of transformed data
summary(as.vector(log_expression))

# 3. Visualize original and transformed expression data
par(mfrow = c(1, 2))

hist(
  expression_data,
  breaks = 50,
  main = "Original Expression",
  xlab = "Expression Values"
)

hist(
  log_expression,
  breaks = 50,
  main = "Log2 Transformed Expression",
  xlab = "Log2 Expression"
)

par(mfrow = c(1, 1))

# 4. Sample-wise boxplot
boxplot(
  log_expression,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Values",
  outline = FALSE,
  las = 2,
  col = "lightblue"
)

# 5. Calculate sample-wise medians
sample_medians <- apply(
  log_expression,
  2,
  median
)

summary(sample_medians)

# 6. Calculate sample-wise means
sample_means <- apply(
  log_expression,
  2,
  mean
)

summary(sample_means)

# 7. Create sample statistics table
sample_statistics <- data.frame(
  Sample = colnames(log_expression),
  Mean = sample_means,
  Median = sample_medians
)

# 8. Export sample-wise statistics
write.csv(
  sample_statistics,
  "results/sample_statistics.csv",
  row.names = FALSE
)
