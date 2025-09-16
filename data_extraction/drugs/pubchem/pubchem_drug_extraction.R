library(AnnotationGx)
library(jsonlite)

# Import drug list
drugs <- read.csv("output_data/drugs.csv")
drug_names <- drugs$drug

# Map drug names to cids, remove drugs that don't have cids
result <- mapCompound2CID(drug_names, first = TRUE)
result <- lapply(result[!is.na(result$cids)], unique)

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

# Column header renaming before creating csv to:  "","cid","molecular_formula","molecular_weight","smiles","connectivity_smiles","inchi","inchikey","iupac_name","xlogp","exact_mass","monoisotopic_mass","tpsa","complexity","charge","h_bond_donor_count","h_bond_acceptor_count","rotatable_bond_count","heavy_atom_count","isotope_atom_count","atom_stereo_count","defined_atom_stereo_count","undefined_atom_stereo_count","bond_stereo_count","defined_bond_stereo_count","undefined_bond_stereo_count","covalent_unit_count","volume_3d","x_steric_quadrupole_3d","y_steric_quadrupole_3d","z_steric_quadrupole_3d","feature_count_3d","feature_acceptor_count_3d","feature_donor_count_3d","feature_anion_count_3d","feature_cation_count_3d","feature_ring_count_3d","feature_hydrophobe_count_3d","conformer_model_rmsd_3d","effective_rotor_count_3d","conformer_count_3d","fingerprint_2d","title","patent_count","patent_family_count","literature_count","annotation_types","annotation_type_count","name","chembl_id","synonyms" 

write.csv(property_data, "output_data/pubchem_table.csv", row.names=TRUE)
