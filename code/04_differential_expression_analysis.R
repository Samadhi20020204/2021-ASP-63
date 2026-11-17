
hist(
  as.vector(expression_data),
  breaks = 50,
  main = "Distribution of Gene Expression Values",
  xlab = "Expression Values",
  col = "skyblue",
  border = "white"
)


boxplot(
  as.vector(expression_data),
  main = "Boxplot of Gene Expression Values",
  ylab = "Expression Values",
  col = "lightgreen"
)


# 4.3 Histogram Frequency Information

hist_data <- hist(
  as.vector(expression_data),
  breaks = 50,
  plot = FALSE
)

hist_data$breaks
hist_data$counts


# 4.4 Boxplot Statistics

boxplot.stats(
  as.vector(expression_data)
)


# 4.5 Log2 Expression Distribution

hist(
  log2(as.vector(expression_data) + 1),
  breaks = 50,
  main = "Log2 Gene Expression Distribution",
  xlab = "Log2 Expression Values",
  col = "lightblue"
)


class(expression_data)

range(expression_data)

summary(as.vector(expression_data))

class(log_expression)

dim(log_expression)

summary(as.vector(log_expression))

range(log_expression)


# 5. Verify Log2 Transformation

all.equal(
  log_expression,
  log2(expression_data + 1)
)


# 6. Sample-wise Expression Distribution Analysis

boxplot(
  log_expression,
  main = "Sample-wise Gene Expression Distribution",
  xlab = "Samples",
  ylab = "Log2 Expression Values",
  las = 2,
  outline = FALSE,
  col = "lightblue",
  cex.axis = 0.5
)


# 6. Sample-wise Expression Distribution Analysis

# Calculate median, mean and standard deviation

sample_distribution <- data.frame(
  Sample_ID = colnames(log_expression),
  Median = apply(log_expression, 2, median),
  Mean = apply(log_expression, 2, mean),
  SD = apply(log_expression, 2, sd)
)

# Display first 10 samples
cat("\n--- First 10 Samples ---\n")
print(head(sample_distribution, 10))

# Overall summary
cat("\n--- Overall Sample Statistics ---\n")
print(summary(sample_distribution[, 2:4]))

# Sample with minimum median
cat("\n--- Sample with Minimum Median ---\n")
print(sample_distribution[which.min(sample_distribution$Median), ])

# Sample with maximum median
cat("\n--- Sample with Maximum Median ---\n")
print(sample_distribution[which.max(sample_distribution$Median), ])


# 7. Check Dataset Information

# Check ExpressionSet class
cat("\n--- ExpressionSet Class ---\n")
print(class(eset))

# Check platform annotation
cat("\n--- Platform Annotation ---\n")
print(Biobase::annotation(eset))

# Check experiment title
cat("\n--- Experiment Title ---\n")
print(Biobase::experimentData(eset)@title)

# Check experiment description
cat("\n--- Experiment Description ---\n")
print(Biobase::experimentData(eset)@abstract)

# Check protocol information
cat("\n--- Protocol Information ---\n")
print(Biobase::protocolData(eset))



# 8. Inspect Original Dataset Processing Information

# Check GEO object
cat("\n--- GEO Object Class ---\n")
print(class(gse))

# Check GEO object structure
cat("\n--- GEO Object Names ---\n")
if (is.list(gse)) {
  print(names(gse))
} else {
  print(slotNames(gse))
}

# Check ExpressionSet additional information
cat("\n--- Additional Experiment Information ---\n")
print(Biobase::experimentData(eset)@other)

# Check sample metadata columns
cat("\n--- Sample Metadata Columns ---\n")
print(colnames(Biobase::pData(eset)))

# Check first 5 samples metadata
cat("\n--- First 5 Samples Metadata ---\n")
print(head(Biobase::pData(eset)[, 1:8], 5))


# 9. Inspect GEO Data Processing Information

sample_metadata <- Biobase::pData(eset)

cat("\n--- Data Processing Information ---\n")
print(unique(sample_metadata$data_processing))

cat("\n--- Sample Description ---\n")
print(head(sample_metadata$description, 5))

cat("\n--- Expression Data Range ---\n")
print(range(expression_data))

cat("\n--- Log Expression Data Range ---\n")
print(range(log_expression))


# 10. Prepare Paired Experimental Design

# Set sample IDs as row names
rownames(sample_info) <- sample_info$geo_accession


# Set group order
sample_info$group <- factor(
  sample_info$group,
  levels = c("Normal", "Cancer")
)

# Convert patient ID into factor
sample_info$patient_id <- factor(
  sample_info$patient_id
)

# Create paired design matrix
design <- model.matrix(
  ~ patient_id + group,
  data = sample_info
)

# Display group information
cat("\n--- Group Information ---\n")
print(table(sample_info$group))

# Check paired samples
cat("\n--- Paired Sample Verification ---\n")
print(table(sample_info$patient_id, sample_info$group))

# Display design matrix dimensions
cat("\n--- Design Matrix Dimensions ---\n")
print(dim(design))

# Display design matrix column names
cat("\n--- Design Matrix Columns ---\n")
print(colnames(design))

# Display first 6 rows
cat("\n--- First 6 Rows of Design Matrix ---\n")
print(head(design))


# 11. Differential Gene Expression Analysis using limma

# Install limma if it is not already installed
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install("limma")
}

library(limma)

# Check sample order
stopifnot(
  identical(colnames(log_expression), rownames(design))
)

# Fit linear model
fit <- lmFit(log_expression, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract differential expression results
deg_results <- topTable(
  fit,
  coef = "groupCancer",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Display top 10 differentially expressed probes
cat("\n--- Top 10 Differentially Expressed Probes ---\n")
print(head(deg_results, 10))

# Identify significant genes
significant_genes <- deg_results[
  deg_results$adj.P.Val < 0.05,
]

# Count significant genes
cat("\n--- Number of Significant Probes ---\n")
print(nrow(significant_genes))

# Display result dimensions
cat("\n--- DEG Results Dimensions ---\n")
print(dim(deg_results))


exists("log_expression")
exists("sample_info")
exists("design")
exists("fit")
exists("deg_results")


# Step 2: Fit linear model
fit <- lmFit(log_expression, design)

# Step 3: Empirical Bayes analysis
fit <- eBayes(fit)

# Step 4: Extract Differentially Expressed Genes
deg_results <- topTable(
  fit,
  coef = "groupCancer",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Step 5: Check results
head(deg_results)
dim(deg_results)
summary(deg_results)

# Step 6: Identify significant genes

significant_genes <- deg_results[
  deg_results$adj.P.Val < 0.05,
]

# Number of significant genes
nrow(significant_genes)

# Display first 10 significant genes
head(significant_genes, 10)


# Step 7: Separate upregulated and downregulated genes

upregulated_genes <- significant_genes[
  significant_genes$logFC > 0,
]

downregulated_genes <- significant_genes[
  significant_genes$logFC < 0,
]

# Count each group
nrow(upregulated_genes)
nrow(downregulated_genes)


# Step 8: Create Volcano Plot

# Calculate -log10 adjusted p-values
volcano_data <- deg_results

volcano_data$neg_log10_padj <- -log10(
  pmax(volcano_data$adj.P.Val, .Machine$double.xmin)
)

# Correct colour assignment

volcano_data$colour <- "grey"

volcano_data$colour[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC > 0
] <- "red"

volcano_data$colour[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC < 0
] <- "blue"

# Save corrected Volcano Plot
png("results/volcano_plot.png",
    width = 1200, height = 900, res = 120)

plot(
  volcano_data$logFC,
  volcano_data$neg_log10_padj,
  col = volcano_data$colour,
  pch = 20,
  cex = 0.6,
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  main = "Volcano Plot: Cancer vs Normal"
)

abline(h = -log10(0.05), lty = 2)
abline(v = 0, lty = 2)

legend(
  "topright",
  legend = c("Upregulated", "Downregulated", "Not Significant"),
  col = c("red", "blue", "grey"),
  pch = 20
)

dev.off()

write.csv(
  deg_results,
  "results/all_differential_expression_results.csv"
)

write.csv(
  significant_genes,
  "results/significant_genes.csv"
)

write.csv(
  upregulated_genes,
  "results/upregulated_genes.csv"
)

write.csv(
  downregulated_genes,
  "results/downregulated_genes.csv"
)


top50_genes <- rownames(
  significant_genes[1:50, ]
)

top50_expression <- log_expression[
  top50_genes,
]


# Create sample annotation

annotation_col <- data.frame(
  Group = sample_info$group
)

rownames(annotation_col) <- rownames(sample_info)

# Arrange annotation according to expression columns

annotation_col <- annotation_col[
  colnames(top50_expression),
  , drop = FALSE
]

# Save heatmap

png("results/heatmap_top50_genes.png",
    width = 1600, height = 1200, res = 120)

pheatmap(
  top50_expression,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()



# Step 14: Install annotation packages

if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
  BiocManager::install("hgu133a.db")
}

library(hgu133a.db)
library(AnnotationDbi)

# Step 15: Extract probe IDs

probe_ids <- rownames(deg_results)

# Map probe IDs to gene symbols

gene_symbols <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Map probe IDs to gene names

gene_names <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)


# Create annotated results table

annotated_results <- deg_results

annotated_results$ProbeID <- rownames(deg_results)

annotated_results$GeneSymbol <- unname(
  gene_symbols[rownames(deg_results)]
)

annotated_results$GeneName <- unname(
  gene_names[rownames(deg_results)]
)

# Arrange columns

annotated_results <- annotated_results[
  ,
  c("ProbeID", "GeneSymbol", "GeneName",
    "logFC", "AveExpr", "t", "P.Value",
    "adj.P.Val", "B")
]



# Display first 10 annotated genes

head(annotated_results, 10)

# Count successfully mapped probes

sum(!is.na(annotated_results$GeneSymbol))

# Count unmapped probes

sum(is.na(annotated_results$GeneSymbol))

# Check dimensions

dim(annotated_results)


write.csv(
  annotated_results,
  "results/annotated_gene_expression_results.csv",
  row.names = FALSE
)



# Get probe IDs from differential expression results
probe_ids <- rownames(deg_results)

# Map probe IDs to gene symbols
gene_symbols <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Map probe IDs to full gene names
gene_names <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)

# Create annotated results table
annotated_results <- deg_results

annotated_results$ProbeID <- rownames(deg_results)

annotated_results$GeneSymbol <- unname(
  gene_symbols[rownames(deg_results)]
)

annotated_results$GeneName <- unname(
  gene_names[rownames(deg_results)]
)

# Arrange columns
annotated_results <- annotated_results[
  ,
  c(
    "ProbeID",
    "GeneSymbol",
    "GeneName",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


head(annotated_results, 10)

# Count successfully annotated probes
sum(!is.na(annotated_results$GeneSymbol))

# Count probes without gene symbols
sum(is.na(annotated_results$GeneSymbol))

# Check dimensions
dim(annotated_results)



write.csv(
  annotated_results,
  "results/annotated_gene_expression_results.csv",
  row.names = FALSE
)

file.exists("results/annotated_gene_expression_results.csv")


significant_annotated <- annotated_results[
  annotated_results$adj.P.Val < 0.05,
]

write.csv(
  significant_annotated,
  "results/annotated_significant_genes.csv",
  row.names = FALSE
)

dim(significant_annotated)
head(significant_annotated, 10)


# Upregulated genes in Cancer
upregulated_annotated <- significant_annotated[
  significant_annotated$logFC > 0,
]

# Downregulated genes in Cancer
downregulated_annotated <- significant_annotated[
  significant_annotated$logFC < 0,
]

# Save both datasets
write.csv(
  upregulated_annotated,
  "results/annotated_upregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  downregulated_annotated,
  "results/annotated_downregulated_genes.csv",
  row.names = FALSE
)


cat("Total significant probes:",
    nrow(significant_annotated), "\n")

cat("Upregulated probes:",
    nrow(upregulated_annotated), "\n")

cat("Downregulated probes:",
    nrow(downregulated_annotated), "\n")

cat("Unique annotated gene symbols:",
    length(unique(na.omit(significant_annotated$GeneSymbol))),
    "\n")

cat("Missing gene symbols:",
    sum(is.na(significant_annotated$GeneSymbol)), "\n")
