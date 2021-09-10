# R script used to retrive SARS-CoV-2 FASTA Sequences

# UniProt identifier of SARS-CoV-2 proteins are retrived from
#https://www.uniprot.org/uniprot/?query=proteome:UP000464024

id <- c("P0DTC2","P0DTD8","P0DTD1","P0DTC6","P0DTC1","P0DTC8",
        "P0DTF1","P0DTC5","P0DTD3","P0DTC3","P0DTG0","P0DTG1","P0DTC7",
        "P0DTD2","P0DTC9","P0DTC4","A0A663DJA2")

local = "C:/Users/Matine/OneDrive/Desktop/cran_packages/HPiP/data"
getFASTA(id,filename = 'UP000464024_df.RData')
