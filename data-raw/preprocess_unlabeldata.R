#R script used to generate unlabeled HP-PPIs

#loading the required packages
library(readxl)
library(stringr)
library(dplyr)
library(HPiP)


dat <-
  read_excel("41586_2020_2286_MOESM5_ESM.xlsx")
viralp <-
  HPiP:::viralprot
colnames(viralp)[2] <- "Bait"
data("Gold_ReferenceSet")


`%nin%` <- Negate(`%in%`)
dat_p <-
  dat %>%
  mutate(Bait = str_remove_all(Bait, "SARS-CoV2\\s+")) %>%
  mutate(Bait = str_replace(Bait,"Spike", "S")) %>%
  mutate(Bait = str_replace(Bait,"orf", "ORF")) %>%
  filter(!str_detect(Bait, "nsp5_C145A")) %>%
  left_join(., viralp) %>%
  mutate(p = paste(`Bait`, `Pathogen_Protein`, sep = ":")) %>%
  mutate(PPI = paste(`p`, `Preys`, sep = "~")) %>%
  dplyr::select(13,11,2) %>%
  dplyr::rename(Host_Protein = Preys)%>%
  dplyr::filter(PPI %nin% Gold_ReferenceSet$PPI)

# retrive sequences (viral)
viralid <- unique(dat_p$Pathogen_Protein)
viralidseq <-getFASTA(viralid)
viralidseq <- do.call(rbind, viralidseq) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("UniprotKB")


#generate descriptors
CTDC_df <- calculateCTDC(viralidseq)
CTDC_df <- CTDC_df[order(match(CTDC_df$identifier,viralid)),]

CTDT_df <- calculateCTDT(viralidseq)
CTDT_df <- CTDT_df[order(match(CTDT_df$identifier,viralid)),]


CTDD_df <- calculateCTDD(viralidseq)
CTDD_df <- CTDD_df[order(match(CTDD_df$identifier,viralid)),]

desc_viral <-
  full_join(CTDC_df,CTDT_df, by = "identifier") %>%
  full_join(., CTDD_df)
cdata <- as.data.frame(colData(viral_se))
desc_viral1 <- as.matrix(desc_viral[, -1])
row.names(desc_viral1) <- desc_viral$identifier





#retrive sequences (host)
hostid <- unique(dat_p$Host_Protein)
hostseq <-getFASTA(hostid)
hostseq <- do.call(rbind, hostseq) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("UniprotKB")


#generate descriptors
CTDC_df <- calculateCTDC(hostseq)
CTDC_df <- CTDC_df[order(match(CTDC_df$identifier,hostid)),]

CTDT_df <- calculateCTDT(hostseq)
CTDT_df <- CTDT_df[order(match(CTDT_df$identifier,hostid)),]

CTDD_df <- calculateCTDD(hostseq)
CTDD_df <- CTDD_df[order(match(CTDD_df$identifier,hostid)),]

desc_host <-
  full_join(CTDC_df,CTDT_df, by = "identifier") %>%
  full_join(., CTDD_df)
cdata <- as.data.frame(colData(host_se))
desc_host1 <- as.matrix(desc_host[, -1])
row.names(desc_host1) <- desc_host$identifier

dat_p <-
  filter(dat_p, Host_Protein %in% desc_host$identifier)

#Map them to the orgininal df
x1_host <- matrix(NA, nrow = nrow(dat_p), ncol = ncol(desc_host1))
for (i in 1:nrow(dat_p))
  x1_host[i, ] <-
  desc_host1[which(dat_p$Host_Protein[i] == desc_host$identifier), ]

#Map them to the orgininal df
x1_viral <- matrix(NA, nrow = nrow(dat_p), ncol = ncol(desc_viral1))
for (i in 1:nrow(dat_p))
  x1_viral[i, ] <- desc_viral1[which(dat_p$Pathogen_Protein[i] == viralid), ]



x <- getHPI(x1_viral,x1_host, type = "combine")
x <- as.data.frame(x)
x <- cbind(dat_p$PPI,x)
colnames(x)[1] <- c("PPI")
unlabel_data <- x
unlabel_data <- sample_n(unlabel_data, 1000)

use_data(unlabel_data,overwrite = TRUE)

