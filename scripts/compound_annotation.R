library(AnnotationGx)

# Define compound names
compound_names <- c(
  "2-Hydroxy-5-(2,5-dihydroxybenzylamino)benzoic acid",
  "5-Nonyloxytryptamine oxalate",
  "6-azauridine",
  "Acetaminophen",
  "AG-126",
  "AG-494",
  "Albendazole",
  "Alfacalcidol",
  "Altretamine"
)

# Step 1: Get CIDs
result <- mapCompound2CID(compound_names, first = TRUE)

# Step 2: Remove compounds that failed to map
result <- result[!is.na(result$cids)]

# Step 3: Define properties to fetch
properties <- c("Title", "MolecularFormula", "MolecularWeight", "InChIKey")

# Step 4: Query PubChem for those properties
property_data <- mapCID2Properties(ids = result$cids, properties = properties)

# Step 5: Add compound names
property_data$name <- result$name

# Step 6: Fetch ChEMBL IDs
property_data$ChEMBL_ID <- annotatePubchemCompound(result$cids, "ChEMBL ID")

# Optional: Reorder columns
property_data <- property_data[, c("name", "CID", "Title", "MolecularFormula", 
                                   "MolecularWeight", "InChIKey", "ChEMBL_ID")]

# Step 7: View or save
print(property_data)
# View(property_data)  # Optional
# write.csv(property_data, "compound_summary.csv", row.names = FALSE)  # Optional 