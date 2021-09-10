# R script used to generate the viral SummerizedExperiment object
#Load the gold-standard data
data('Gold_ReferenceSet')
gd <- Gold_ReferenceSet

#retrive sequences
viralid <- unique(gd$Pathogen_Protein)
viralidseq <-getFASTA(viralid)
viralidseq <- do.call(rbind, viralidseq) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("UniprotKB")


#generate descriptors
CTDC_df <- calculateCTDC(viralidseq)
CTDC_df <- CTDC_df[order(match(CTDC_df$identifier,viralid)),]
CTDC_m <- as.matrix(CTDC_df[, -1])
row.names(CTDC_m) <- CTDC_df$identifier

CTDT_df <- calculateCTDT(viralidseq)
CTDT_df <- CTDT_df[order(match(CTDT_df$identifier,viralid)),]
CTDT_m <- as.matrix(CTDT_df[, -1])
row.names(CTDT_m) <- CTDT_df$identifier

CTDD_df <- calculateCTDD(viralidseq)
CTDD_df <- CTDD_df[order(match(CTDD_df$identifier,viralid)),]
CTDD_m <- as.matrix(CTDD_df[, -1])
row.names(CTDD_m) <- CTDD_df$identifier



#convert df to se object

ctdc_se <- SummarizedExperiment(assays = list(counts = CTDC_m),
                                colData = paste0(colnames(CTDC_df[,-1]),
                                                 "CTDC"))
ctdt_se <- SummarizedExperiment(assays = list(counts = CTDT_m),
                                colData = paste0(colnames(CTDT_df[,-1]),
                                                 "CTDT"))
ctdd_se <- SummarizedExperiment(assays = list(counts = CTDD_m),
                                colData = paste0(colnames(CTDD_df[,-1]),
                                                 "CTDD"))

viral_se <- cbind(ctdc_se,ctdd_se,ctdt_se)

use_data(viral_se,overwrite = TRUE)


