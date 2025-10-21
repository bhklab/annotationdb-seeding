# Used to remove duplicates by CIDs (if they exist)
library(dplyr)

pubchem_table <- read.csv("output_data/pubchem_table.csv")
keep_unique_cid <- pubchem_table %>% distinct(cid, .keep_all=TRUE)

keep_unique_cid$X <- NULL

write.csv(keep_unique_cid, "output_data/pubchem_table_omitted.csv")
