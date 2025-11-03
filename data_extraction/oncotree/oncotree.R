library(AnnotationGx)

df <- as.data.frame(getOncotreeTumorTypes())
list_cols <- names(Filter(is.list, df))

df[list_cols] <- lapply(df[list_cols], function(col) {
  vapply(col, function(x) {
    if (is.null(x) || length(x) == 0) return(NA_character_)
    if (is.list(x)) paste(unlist(x, use.names = FALSE), collapse = ",")
    else paste(x, collapse = ",")
  }, character(1))
})

 # Rename columns to align schema
names(df) <- sub("mainType", "main_type", names(df))
names(df) <- sub("parent", "parent_code", names(df))
names(df) <- sub("parentCode", "parent_code", names(df))
names(df) <- sub("externalReferences", "external_references", names(df))

ver <- as.data.frame(getOncotreeVersions())

strip_html <- function(s) gsub("<[^>]+>", "", s)

latest_row <- ver[ trimws(strip_html(ver$description)) ==
                     "This is the latest approved version for public use.", ]

version_api_identifier <- latest_row$api_identifier[1]
version_release_date   <- as.character(as.Date(latest_row$release_date[1]))

df$version_api_identifier <- version_api_identifier
df$version_release_date   <- version_release_date

write.csv(df, "output_data/oncotree.csv", row.names = FALSE, fileEncoding = "UTF-8")
