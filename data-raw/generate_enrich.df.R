# R script used to generate the enrichment result
library(HPiP)


data('predicted_PPIs')

#perform enrichment
enrich.df <- enrichfindP(predicted_PPIs,threshold = 0.05,
                         sources = c("GO", "KEGG"),
                         p.corrction.method = "bonferroni",
                         org = "hsapiens")

use_data(enrich.df, overwrite = TRUE)




