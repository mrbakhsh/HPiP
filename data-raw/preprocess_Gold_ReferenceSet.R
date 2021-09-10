# R script used to generate `Gold_ReferenceSet`


# Loading required packages
library(httr)
library(dplyr)
library(tidyr)


# Sars-cov-2 taxanomy ID
organism.taxID = 2697049
access.key = '81bb3b5a6bd9a8084a7be71f0963ab1e'
# construct query to retrive all interactiosn for the requested organism
url <-
  "http://webservice.thebiogrid.org/interactions/?"
query <- paste(url,
               paste0('accessKey=',access.key),
               "interSpeciesExcluded=false",
               "selfInteractionsExcluded=true",
               "includeHeader=true",
               paste0("taxId=", organism.taxID),
               sep = "&"
)

# number of results
n <-
  httr::GET(paste(query, "format=count", sep = "&"))
httr::stop_for_status(n)


n <- as.numeric(
  httr::content(n, type = "text/csv", encoding = "UTF-8", col_types = "i")
)

# retrieve results in batches
tpINT <- lapply(seq(1, n, 10000), function(m) {
  l <- paste(query, paste0("start=", m), sep = "&")
  request <-
    httr::GET(l)
  httr::stop_for_status(request)

  # parse the biogird results
  biog_result <- httr::content(request,
                               type = "text/tab-separated-values",
                               encoding = "UTF-8",
                               col_types = "cccccccccccccccccccccccc"
  )
  q <- subset(biog_result, select = c(
    "Official Symbol Interactor A",
    "Official Symbol Interactor B",
    "Organism Interactor A",
    "Organism Interactor B",
    "Experimental System Type",
    "Author"
  ))

  return(q)
})

tpINT <- do.call("rbind", tpINT)
dataTP <-
  tpINT
dataTP <- # remove duplicate
  tpINT[!duplicated(apply(tpINT, 1, function(x)
    paste(sort(x), collapse = ""))), ]

dataTP <- #select Samavarchi-Tehrani P (2020) PPIs
  dataTP %>%
  filter(Author == "Samavarchi-Tehrani P (2020)") %>%
  select(1,2) %>%
  #sample_n(., 500) %>%
  mutate(PPI = paste(`Official Symbol Interactor A`,
                     `Official Symbol Interactor B`, sep = "~")) %>%
  select(3,1,2)
dataTP$class <- "Positive"

#generate negative PPIs
path_prot <- unique(dataTP$`Official Symbol Interactor A`)
hostprot <- unique(dataTP$`Official Symbol Interactor B`)
TP <- dataTP$PPI
TN <- get_negativePPI(path_prot,hostprot,TP)
TN <-
  separate(TN, PPI,
           c("Official Symbol Interactor A", "Official Symbol Interactor B"),
           sep = "~", remove = FALSE)
TN$class <- "Negative"

gd <- rbind(dataTP, TN)


#map gene names to UniProt

#open viral UniProt identifier (this file is directly retrieved from UniProt)
viralprot <- HPiP:::viralprot
gd1 <- left_join(gd, viralprot)
#open host UniProt identifier (this file is directly retrieved from UniProt)
hostprot <- HPiP:::hostprot
gd2 <- left_join(gd1, hostprot)
gd2 <- na.omit(gd2)

gd2 <-
  unite(gd2, PPI,
        c(`Official Symbol Interactor A`, Pathogen_Protein),
        sep = ":", remove = FALSE)
gd2 <-
  unite(gd2, PPI,  c(PPI, Host_Protein), sep = "~", remove = FALSE)

set.seed(101)
positive <-
  gd2 %>%
  dplyr::filter(class == "Positive") %>%
  dplyr::select(2,1,3,5,6,4) %>%
  sample_n(500)

set.seed(102)
negative <-
  gd2 %>%
  dplyr::filter(class == "Negative") %>%
  dplyr::select(2,1,3,5,6,4) %>%
  sample_n(500)



Gold_ReferenceSet <-
  rbind(positive, negative)
use_data(Gold_ReferenceSet, overwrite = TRUE)

