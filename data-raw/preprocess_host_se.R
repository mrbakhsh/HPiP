# R script used to generate the host SummerizedExperiment object
#Load the gold-standard data
data('Gold_ReferenceSet')
gd <- Gold_ReferenceSet

#retrive sequences
hostid <- unique(gd$Host_Protein)
hostseq <-getFASTA(hostid)
hostseq <- do.call(rbind, hostseq) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("UniprotKB")


#generate descriptors
CTDC_df <- calculateCTDC(hostseq)
CTDC_df <- CTDC_df[order(match(CTDC_df$identifier,hostid)),]
CTDC_m <- as.matrix(CTDC_df[, -1])
row.names(CTDC_m) <- CTDC_df$identifier

CTDT_df <- calculateCTDT(hostseq)
CTDT_df <- CTDT_df[order(match(CTDT_df$identifier,hostid)),]
CTDT_m <- as.matrix(CTDT_df[, -1])
row.names(CTDT_m) <- CTDT_df$identifier

CTDD_df <- calculateCTDD(hostseq)
CTDD_df <- CTDD_df[order(match(CTDD_df$identifier,hostid)),]
CTDD_m <- as.matrix(CTDD_df[, -1])
row.names(CTDD_m) <- CTDD_df$identifier

ctdc_se <- SummarizedExperiment(assays = list(counts = CTDC_m),
                                colData = paste0(colnames(CTDC_df[,-1]),
                                                 "CTDC"))
ctdt_se <- SummarizedExperiment(assays = list(counts = CTDT_m),
                                colData = paste0(colnames(CTDT_df[,-1]),
                                                 "CTDT"))
ctdd_se <- SummarizedExperiment(assays = list(counts = CTDD_m),
                                colData = paste0(colnames(CTDD_df[,-1]),
                                                 "CTDD"))

host_se <-
  cbind(ctdc_se,ctdd_se,ctdt_se)

use_data(host_se,overwrite = TRUE)

