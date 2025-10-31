library(AnnotationGx)

# --- 1) Pull tumor types and flatten list-columns ---
df <- as.data.frame(getOncotreeTumorTypes())
list_cols <- names(Filter(is.list, df))

df[list_cols] <- lapply(df[list_cols], function(col) {
  vapply(col, function(x) {
    if (is.null(x) || length(x) == 0) return(NA_character_)
    if (is.list(x)) paste(unlist(x, use.names = FALSE), collapse = ",")
    else paste(x, collapse = ",")
  }, character(1))
})

# --- 2) Get the current OncoTree version metadata ---
ver <- as.data.frame(getOncotreeVersions())

# Strip any HTML from descriptions and match the required sentence
strip_html <- function(s) gsub("<[^>]+>", "", s)

latest_row <- ver[ trimws(strip_html(ver$description)) ==
                     "This is the latest approved version for public use.", ]

version_api_identifier <- latest_row$api_identifier[1]
version_release_date   <- as.character(as.Date(latest_row$release_date[1]))

# --- 3) Append version fields to every tumor-type row ---
df$version_api_identifier <- version_api_identifier
df$version_release_date   <- version_release_date

# --- 4) Write CSV ---
write.csv(df, "output_data/oncotree.csv", row.names = FALSE, fileEncoding = "UTF-8")

