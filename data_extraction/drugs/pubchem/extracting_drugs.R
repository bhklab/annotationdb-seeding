#BeatAML_2018 has cooked drug names, there will be a hard time mapping those (57,279 unique drugs as well)

library(PharmacoGx)

# Data directories
data_dir <- file.path(getwd(), "../rawdata/orcestradata/pharmacosets")
output <- file.path(getwd(), "../output_data/drugs.csv")


# Initialize drugs data structure, and load rds files

all_drugs <- character(0)
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)

# Iterate RDS files, extract drugs either through treatment (new structure) or drug (old structure)
for (f in rds_files) {
	
	pset <- readRDS(f)
	if (is.null(pset)) next

	pset <- updateObject(pset)
	drugs <- rownames(pset@treatment)

	if (is.null(drugs)) {
		drugs <- pset@drug$drug.name
	}

	if (!is.null(drugs)) {
		all_drugs <- unique(c(all_drugs, drugs))
	}
}

all_drugs <- sort(all_drugs)
out_df <- data.frame(drug = all_drugs, stringsAsFactors = FALSE)
write.csv(out_df, file = output, row.names = FALSE)


