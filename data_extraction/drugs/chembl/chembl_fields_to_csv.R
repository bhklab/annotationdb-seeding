library(AnnotationGx)
library(jsonlite)

# Get all chembl object names
chembl_resources <- getChemblResources()

# Extract all object fields
field_list <- list()
for (res in chembl_resources) {
	if (res != "document_term") {
		field_list[[res]] = getChemblResourceFields(res)
	}
}


# Output to JSON
json_out <- toJSON(field_list, pretty = TRUE, auto_unbox = TRUE)
write(json_out, file = "chembl_resources_fields.json")


# Find longest field vector
max_len <- 0
for (res in names(field_list)) {
	n <- length(field_list[[res]])
	if (n > max_len) max_len <- n
}

# Add NA to to fill missing fields below max_len, utilizes R idiom
for (res in names(field_list)){
	length(field_list[[res]]) <- max_len
}

# Output to csv
fields_df <- data.frame(field_list, stringsAsFactors = FALSE)
write.csv(fields_df, "chembl_resources_fields.csv", row.names = FALSE)