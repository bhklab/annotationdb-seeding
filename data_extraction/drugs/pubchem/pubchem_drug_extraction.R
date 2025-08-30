library(AnnotationGx)
library(jsonlite)

# Import drug list
drugs <- read.csv("output_data/drugs.csv")
drug_names <- drugs$drug

# Map drug names to cids, remove drugs that don't have cids
result <- mapCompound2CID(drug_names, first = TRUE)
result <- result[!is.na(result$cids)]

# List wanted pubchem properties, map cids to wanted properties
properties <- c("Title", "MolecularFormula", "MolecularWeight", "SMILES", "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass", "TPSA", "Complexity", "Charge", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount", "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount", "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount", "UndefinedBondStereoCount", "CovalentUnitCount", "Volume3D", "XStericQuadrupole3D", "YStericQuadrupole3D", "ZStericQuadrupole3D", "FeatureCount3D", "FeatureAcceptorCount3D", "FeatureDonorCount3D", "FeatureAnionCount3D", "FeatureCationCount3D", "FeatureRingCount3D", "FeatureHydrophobeCount3D", "ConformerModelRMSD3D", "EffectiveRotorCount3D", "ConformerCount3D", "Fingerprint2D", "PatentCount", "PatentFamilyCount", "LiteratureCount", "AnnotationTypes", "AnnotationTypeCount")
property_data <- mapCID2Properties(ids = result$cids, properties = properties)

# Mapping drug names and ChEMBL ID to pubchem data 
property_data$name <- result$name
property_data$ChEMBL_ID <- annotatePubchemCompound(result$cids, "ChEMBL ID")

# Gather synonyms for drug names
property_data$synonyms <- character(nrow(property_data))
for (i in seq_along(property_data$CID)) {
	cid <- property_data$CID[i]
	if (is.na(cid)) { 
		next
	}
	
	url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/TXT", cid)
  	synonyms <- tryCatch(unique(readLines(url, warn = FALSE)), error = function(e) character(0))

	property_data$synonyms[i] <- if (length(synonyms)) paste(synonyms, collapse = ", ") else "" # Split list into comma separated string
	
	Sys.sleep(0.1)  # Small API buffer to not get blocked
}

write.csv(property_data, "output_data/pubchem_table.csv", row.names=TRUE)
