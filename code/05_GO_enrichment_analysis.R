
# ============================================================
# 05_GO_enrichment_analysis.R
# Gene Ontology Biological Process Enrichment Analysis
# Dataset: GSE15852
# ============================================================

# 1. Load required packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# 2. Check required analysis objects

required_objects <- c(
  "up_entrez",
  "down_entrez",
  "universe_entrez"
)

missing_objects <- required_objects[
  !sapply(required_objects, exists)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects: ",
    paste(missing_objects, collapse = ", "),
    ". Run the previous analysis and ID mapping steps first."
  )
}

# 3. Create results directory

dir.create("results", showWarnings = FALSE)

# 4. GO Biological Process enrichment
# Upregulated genes

go_up <- enrichGO(
  gene = up_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# 5. GO Biological Process enrichment
# Downregulated genes

go_down <- enrichGO(
  gene = down_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# 6. Convert results into data frames

go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

# 7. Display summary

cat("Upregulated GO terms:", nrow(go_up_df), "\n")
cat("Downregulated GO terms:", nrow(go_down_df), "\n")

# 8. Save GO enrichment results

write.csv(
  go_up_df,
  "results/GO_upregulated_results.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "results/GO_downregulated_results.csv",
  row.names = FALSE
)

# 9. Generate dot plots

if (nrow(go_up_df) > 0) {
  
  up_plot <- dotplot(
    go_up,
    showCategory = 15,
    title = "GO Biological Processes - Upregulated Genes"
  )
  
  ggsave(
    "results/GO_upregulated_dotplot.png",
    plot = up_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  print(up_plot)
}

if (nrow(go_down_df) > 0) {
  
  down_plot <- dotplot(
    go_down,
    showCategory = 15,
    title = "GO Biological Processes - Downregulated Genes"
  )
  
  ggsave(
    "results/GO_downregulated_dotplot.png",
    plot = down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  print(down_plot)
}

# 10. Final output verification

cat("\nGO analysis completed.\n")

print(list.files("results", pattern = "^GO_"))
