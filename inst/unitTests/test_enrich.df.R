test_enrich.df <- function() {
  data('predicted_PPIs')
  enrich.df <- enrichfindP(predicted_PPIs,
                           threshold = 0.05,
                           sources = c("GO", "KEGG"),
                           p.corrction.method = "bonferroni",
                           org = "hsapiens")

  checkTrue(is.data.frame(enrich.df) == TRUE)
}

