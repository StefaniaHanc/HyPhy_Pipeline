#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(tidyverse) 

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error : Provide target genes file and GMT file as arguments.\nUsage: Rscript path/07_GOE.R <input_file.tsv> <gmt_file.gmt>", call. = FALSE)
}

INPUT_FILE <- args[1]
FILE_GMT <- args[2]

base_name <- tools::file_path_sans_ext(basename(INPUT_FILE))
output_dir <- dirname(INPUT_FILE)

################################################################

# signif genes
genes_input <- read_table(INPUT_FILE, col_names = TRUE)
target_genes <- gsub('"', '', genes_input[[1]])

################################################################

# GMT file
clean_gmt <- read.gmt(FILE_GMT) %>% 
  filter(term != "") %>% 
  separate_rows(term, sep = "\\|")

################################################################

# split GO terms to BP, MF, CC

# map to ontology 
go_ontologies <- Ontology(unique(clean_gmt$term))

go_ontology_df <- data.frame(
  term = names(go_ontologies),
  ontology = as.character(go_ontologies),
  stringsAsFactors = FALSE
)

# merge with main annot
clean_gmt_with_onto <- merge(clean_gmt, go_ontology_df, by = "term")

# GMT df for each category
gmt_BP <- clean_gmt_with_onto[clean_gmt_with_onto$ontology == "BP", c("term", "gene")]
gmt_MF <- clean_gmt_with_onto[clean_gmt_with_onto$ontology == "MF", c("term", "gene")]
gmt_CC <- clean_gmt_with_onto[clean_gmt_with_onto$ontology == "CC", c("term", "gene")]

################################################################

# GO term to name dict
unique_go_terms <- unique(clean_gmt$term)

term2name <- data.frame(
  term = unique_go_terms,
  name = Term(unique_go_terms)
)

################################################################

# ENRICHER

# keep cutoffs at 1 to fetch all terms in order to display all
results_BP <- enricher(gene = target_genes, TERM2GENE = gmt_BP, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)
results_MF <- enricher(gene = target_genes, TERM2GENE = gmt_MF, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)
results_CC <- enricher(gene = target_genes, TERM2GENE = gmt_CC, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)


################################################################

# plots

plot_all_go_terms <- function(enrich_result, title_str = "GO Enrichment", showCategory = 25) {
  
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    return(ggplot() + labs(title = title_str) + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "No terms mapped"))
  }
  
  # extract df & sort by p-value
  df <- as.data.frame(enrich_result)
  df <- df[order(df$pvalue), ]
  
  # limit to top N
  if (nrow(df) > showCategory) {
    df <- df[1:showCategory, ]
  }
  
  # numeric Gene ratio
  df <- df %>%
    mutate(
      num = as.numeric(sub("/.*", "", GeneRatio)),
      denom = as.numeric(sub(".*/", "", GeneRatio)),
      GeneRatio_num = num / denom
    ) %>%
    mutate(Significant = ifelse(pvalue < 0.05, "p < 0.05", "p >= 0.05"))
  
  # truncate too long descriptions
  #df$Description <- str_trunc(df$Description, 40)
  
  # sort by Gene ratio value for vertical ordering
  df$Description <- factor(df$Description, levels = unique(df$Description[order(df$GeneRatio_num)]))
  
  # ggplot dotplot
  p <- ggplot(df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = pvalue, shape = Significant), stroke = 1.2) +
    scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1), drop = FALSE, name = "Significance") +
    scale_color_gradient(low = "red", high = "blue") + 
    labs(
      title = title_str,
      x = "Gene ratio",
      y = "GO description",
      size = "Gene count",
      color = "p-value"
    ) +
    
    #legend order 1. significance, 2. count, 3. p-val
    guides(
      shape = guide_legend(order = 2),
      size = guide_legend(order = 3),
      color = guide_colorbar(order = 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}

################################################################

# plot & save

p_bp <- plot_all_go_terms(results_BP, title_str = paste(base_name, "- Biological process"))
p_mf <- plot_all_go_terms(results_MF, title_str = paste(base_name, "- Molecular function"))
p_cc <- plot_all_go_terms(results_CC, title_str = paste(base_name, "- Cellular component"))

ggsave(file.path(output_dir, paste0(base_name, "_BP.png")), plot = p_bp, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0(base_name, "_MF.png")), plot = p_mf, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0(base_name, "_CC.png")), plot = p_cc, width = 8, height = 6, dpi = 300)

# write.csv(as.data.frame(results_BP), file = paste0(base_name, "_ClusterProfiler_BP_results.csv"), row.names = FALSE)
# write.csv(as.data.frame(results_MF), file = paste0(base_name, "_ClusterProfiler_MF_results.csv"), row.names = FALSE)
# write.csv(as.data.frame(results_CC), file = paste0(base_name, "_ClusterProfiler_CC_results.csv"), row.names = FALSE)