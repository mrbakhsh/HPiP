test_getFSTA = function() {
  library(protr)
  uniprot.id <- c("P0DTC4", "P0DTC5", "P0DTC9")
  local = tempdir()
  fasta_df <- getFASTA(uniprot.id,
                              filename = 'FASTA.RData', path = local)
  checkTrue(length(fasta_df) == 3)

}


