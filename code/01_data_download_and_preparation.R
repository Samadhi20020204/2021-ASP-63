
# ============================================================
# Project: Differential Gene Expression Analysis
# Dataset: GSE15852
# Task: Data Download and Preparation
# ============================================================

# 1. Create project directories

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 2. Load required libraries

library(GEOquery)
library(Biobase)

# 3. Download GEO dataset

gse <- getGEO(
  "GSE15852",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data/raw"
)

# 4. Extract ExpressionSet

eset <- gse[[1]]

# Save a local copy of the ExpressionSet
# This file is excluded from GitHub using .gitignore

saveRDS(
  eset,
  "data/raw/GSE15852_ExpressionSet.rds"
)

# 5. Extract expression data and sample information

expression_data <- exprs(eset)
sample_info <- pData(eset)

# 6. Prepare sample metadata

sample_info$sample_id <- rownames(sample_info)

sample_info$group <- ifelse(
  grepl("^Normal", sample_info$title),
  "Normal",
  "Cancer"
)

sample_info$patient_id <- sub(
  "^(Normal|Cancer)\\s+",
  "",
  sample_info$title
)

sample_info$patient_id <- sub(
  "[NT]$",
  "",
  sample_info$patient_id
)

# 7. Validate sample groups and patient pairs

print(table(sample_info$group))

print(table(
  sample_info$patient_id,
  sample_info$group
))

# 8. Check sample alignment

stopifnot(
  all(colnames(expression_data) == rownames(sample_info))
)

# 9. Save processed data

saveRDS(
  expression_data,
  "data/processed/expression_data.rds"
)

write.csv(
  sample_info,
  "data/processed/sample_info.csv",
  row.names = FALSE
)

# 10. Display dataset dimensions

print(dim(expression_data))
print(dim(sample_info))

# ============================================================
# End of data preparation
# ============================================================
