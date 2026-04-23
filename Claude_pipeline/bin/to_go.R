library(tidyverse)

# 1. Load the data
go_res <- read_tsv("/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results_interactive/goseq_GO_results.tsv", show_col_types = FALSE)

# 2. Filter and prepare the top 15 Biological Processes
top_bp <- go_res |>
  filter(ontology == "BP", padj < 0.05) |>    # Keep only significant Biological Processes
  arrange(padj) |>                            # Sort by most significant
  slice_head(n = 15) |>                       # Take the top 15
  mutate(
    # Calculate Gene Ratio (DE genes in term / Total genes in term)
    gene_ratio = numDEInCat / numInCat,
    # Truncate very long GO terms so they fit on the plot
    term_short = str_trunc(term, width = 50),
    # Order the terms by significance so the plot is neatly sorted
    term_short = fct_reorder(term_short, -padj) 
  )

# 3. Generate the Dot Plot
ggplot(top_bp, aes(x = gene_ratio, y = term_short, size = numDEInCat, color = padj)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red",high = "blue") + # Inverse colors: lower p-value = brighter
  labs(
    title = "Top 15 Enriched Biological Processes",
    subtitle = "Length-bias corrected (goseq)",
    x = "Gene Ratio (DE genes / Total genes in term)",
    y = NULL,
    size = "DE Genes\n(Count)",
    color = "FDR (padj)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed", color = "grey80"),
    axis.text.y = element_text(color = "black")
  )

ggsave("/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results_interactive/goseq_dotplot.svg", width = 8, height = 6)
