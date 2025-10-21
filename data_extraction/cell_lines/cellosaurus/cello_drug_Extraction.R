library(AnnotationGx)
library(data.table)
library(jsonlite)

# Import drug list
cell_lines_file <- read.csv("../output_data/cellosaurus_ids.csv")
cell_lines <- cell_lines_file$cell_line


# Map cell line anmes to cellosaurus accession id, remove NA's
# result <- mapCell2Accession(cell_lines, parsed = TRUE)
# result <- result[!is.na(result$accession)]

flatten_col <- function(x) {
  if (!is.list(x)) return(x)
  vapply(
    x,
    function(el) {
      if (is.null(el) || length(el) == 0) {
        NA_character_
      } else if (is.list(el) && !is.null(names(el))) {
        # Named/nested list: keep structure as JSON
        jsonlite::toJSON(el, auto_unbox = TRUE, null = "null")
      } else {
        # Vector/unnamed list: collapse unique values
        paste(unique(as.character(unlist(el))), collapse = " | ")
      }
    },
    character(1)
  )
}


# List wanted cellosaurus fields
fields <- c("id", "sy", "idsy", "ac", "acas", "dr", "ww", "genome-ancestry", "hla", "registration", "sequence-variation", "anecdotal", "biotechnology", "breed", "caution", "cell-type", "characteristics", "donor-info", "derived-from-site", "discontinued", "doubling-time", "from", "group", "karyotype", "knockout", "msi", "miscellaneous", "misspelling", "mab-isotype", "mab-target", "omics", "part-of", "population", "problematic", "resistance", "senescence", "integrated", "transformant", "virology", "cc", "str", "di", "din", "dio", "ox", "sx", "ag", "oi", "hi", "ch", "ca", "dt", "dtc", "dtu", "dtv")
field_data <- annotateCellAccession(cell_lines, to = fields)
field_data_flat <- as.data.frame(lapply(field_data, flatten_col), stringsAsFactors = FALSE)

write.csv(field_data_flat, "../output_data/cell_line_table.csv", row.names=TRUE)
