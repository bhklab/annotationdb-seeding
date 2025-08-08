import os
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages

# Extract unique drugs from csv
drugs_df = pd.read_csv(f'{os.getcwd()}/seeding/input_data/drugs.csv', encoding='utf-8')
drugs_list = drugs_df["drug"].to_list()

drugs = set(drugs_list)

# Extract unique cell lines from csv
cells_df = pd.read_csv(f'{os.getcwd()}/seeding/input_data/drugs.csv', encoding='utf-8')
cells_list = drugs_df["cells"].to_list()

cells = set(drugs_list)


# AnnotationGx Functions to get PubChem, ChEMBL, and Cellosaurus data

robjects.r('install.packages("BiocManager", repos="http://cran.r-project.org")')
robjects.r('BiocManager::install("AnnotationGx")')

annotationgx = rpackages.importr("AnnotationGx")
print(annotationgx.package_version("AnnotationGx"))

'''
	Utilizing AnnotationGx to extract compound data from Pubchem to populate SQL table
'''

robjects.r('result <- mapCompound2CID(compound_names, first = TRUE)') # Map Compounds to their PubChem CIDs
robjects.r('result <- result[!is.na(result$cids)]') # Remove failed mappings
robjects.r('properties <- c("Title", "MolecularFormula", "MolecularWeight", "SMILES", "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass", "TPSA", "Complexity", "Charge", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount", "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount", "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount", "UndefinedBondStereoCount", "CovalentUnitCount", "Volume3D", "XStericQuadrupole3D", "YStericQuadrupole3D", "ZStericQuadrupole3D", "FeatureCount3D", "FeatureAcceptorCount3D", "FeatureDonorCount3D", "FeatureAnionCount3D", "FeatureCationCount3D", "FeatureRingCount3D", "FeatureHydrophobeCount3D", "ConformerModelRMSD3D", "EffectiveRotorCount3D", "ConformerCount3D", "Fingerprint2D", "Title", "PatentCount", "PatentFamilyCount", "LiteratureCount", "AnnotationTypes", "AnnotationTypeCount")') # Choose properties to be extracted from PubChem

robjects.r('property_data <- mapCID2Properties(ids = result$cids, properties = properties)') # Get data for each property of interest

robjects.r('property_data$name <- result$name') # Add compound names to result
robjects.r('property_data$ChEMBL_ID <- annotatePubchemCompound(result$cids, "ChEMBL ID")')# Add ChEMBL ids to result

robjects.r('write.csv(property_data), "pubchem_table.csv", row.names=TRUE')

'''
	Utilizing AnnotationGx to extract cell line data from cellosaurus to populate SQL table
'''

robjects.r('acc <- mapCell2Accession(cell_names)')
robjects.r('to <- c("id", "sy", "idsy", "ac", "acas", "dr", "ww", "genome-ancestry", "hla", "registration", "sequence-variation", "anecdotal", "biotechnology", "breed", "caution", "cell-type", "characteristics", "donor-info", "derived-from-site", "discontinued", "doubling-time", "from", "group", "karyotype", "knockout", "msi", "miscellaneous", "misspelling", "mab-isotype", "mab-target", "omics", "part-of", "population", "problematic", "resistance", "senescence", "integrated", "transformant", "virology", "cc", "str", "di", "din", "dio", "ox", "sx", "ag", "oi", "hi", "ch", "ca", "dt", "dtc", "dtu", "dtv")')
robjects.r('result <- annotateCellAccession(acc$accession, to=to, query_only=FALSE, raw=FALSE)')

