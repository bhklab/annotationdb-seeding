# ---- Setup ----
library(AnnotationGx)
library(jsonlite)
library(data.table)
library(future.apply)
library(filelock)
library(httr)


# Import drug list
drugs <- read.csv("output_data/drugs.csv")
drug_names <- drugs$drug

# Map drug names to cids, remove drugs that don't have cids
result <- mapCompound2CID(drug_names, first = TRUE)
result <- result[!is.na(result$cids), , drop = FALSE]

# Desired PubChem properties
properties <- c(
  "Title","MolecularFormula","MolecularWeight","SMILES","ConnectivitySMILES",
  "InChI","InChIKey","IUPACName","XLogP","ExactMass","MonoisotopicMass","TPSA",
  "Complexity","Charge","HBondDonorCount","HBondAcceptorCount","RotatableBondCount",
  "HeavyAtomCount","IsotopeAtomCount","AtomStereoCount","DefinedAtomStereoCount",
  "UndefinedAtomStereoCount","BondStereoCount","DefinedBondStereoCount",
  "UndefinedBondStereoCount","CovalentUnitCount","Volume3D","XStericQuadrupole3D",
  "YStericQuadrupole3D","ZStericQuadrupole3D","FeatureCount3D","FeatureAcceptorCount3D",
  "FeatureDonorCount3D","FeatureAnionCount3D","FeatureCationCount3D","FeatureRingCount3D",
  "FeatureHydrophobeCount3D","ConformerModelRMSD3D","EffectiveRotorCount3D",
  "ConformerCount3D","Fingerprint2D","PatentCount","PatentFamilyCount","LiteratureCount",
  "AnnotationTypes","AnnotationTypeCount"
)

# Output files
dir.create("output_data", showWarnings = FALSE, recursive = TRUE)
prop_out <- "output_data/drug_out.csv"
syn_out  <- "output_data/drug_synonyms.csv"
lockfile <- "output_data/.csv.lock"

# Column rename map (same as yours)
rename_map <- c(
  CID = "cid",
  Title = "title",
  name = "mapped_name",
  ChEMBL_ID = "molecule_chembl_id",
  MolecularFormula = "molecular_formula",
  MolecularWeight = "molecular_weight",
  SMILES = "smiles",
  ConnectivitySMILES = "connectivity_smiles",
  InChI = "inchi",
  InChIKey = "inchikey",
  IUPACName = "iupac_name",
  XLogP = "xlogp",
  ExactMass = "exact_mass",
  MonoisotopicMass = "monoisotopic_mass",
  TPSA = "tpsa",
  Complexity = "complexity",
  Charge = "charge",
  HBondDonorCount = "h_bond_donor_count",
  HBondAcceptorCount = "h_bond_acceptor_count",
  RotatableBondCount = "rotatable_bond_count",
  HeavyAtomCount = "heavy_atom_count",
  IsotopeAtomCount = "isotope_atom_count",
  AtomStereoCount = "atom_stereo_count",
  DefinedAtomStereoCount = "defined_atom_stereo_count",
  UndefinedAtomStereoCount = "undefined_atom_stereo_count",
  BondStereoCount = "bond_stereo_count",
  DefinedBondStereoCount = "defined_bond_stereo_count",
  UndefinedBondStereoCount = "undefined_bond_stereo_count",
  CovalentUnitCount = "covalent_unit_count",
  Volume3D = "volume_3d",
  XStericQuadrupole3D = "x_steric_quadrupole_3d",
  YStericQuadrupole3D = "y_steric_quadrupole_3d",
  ZStericQuadrupole3D = "z_steric_quadrupole_3d",
  FeatureCount3D = "feature_count_3d",
  FeatureAcceptorCount3D = "feature_acceptor_count_3d",
  FeatureDonorCount3D = "feature_donor_count_3d",
  FeatureAnionCount3D = "feature_anion_count_3d",
  FeatureCationCount3D = "feature_cation_count_3d",
  FeatureRingCount3D = "feature_ring_count_3d",
  FeatureHydrophobeCount3D = "feature_hydrophobe_count_3d",
  ConformerModelRMSD3D = "conformer_model_rmsd_3d",
  EffectiveRotorCount3D = "effective_rotor_count_3d",
  ConformerCount3D = "conformer_count_3d",
  Fingerprint2D = "fingerprint_2d",
  PatentCount = "patent_count",
  PatentFamilyCount = "patent_family_count",
  LiteratureCount = "literature_count",
  AnnotationTypes = "annotation_types",
  AnnotationTypeCount = "annotation_type_count"
)

# Resume support: skip CIDs already processed into prop_out
already_done <- character(0)
if (file.exists(prop_out) && isTRUE(file.info(prop_out)$size > 0)) {
  existing <- tryCatch(data.table::fread(prop_out, nThread = 1, showProgress = FALSE),
                       error = function(e) NULL)
  if (!is.null(existing) && "cid" %in% names(existing)) {
    already_done <- as.character(unique(existing$cid))
  }
}


# Merge names with CIDs (so each worker knows the mapped_name)
cid_df <- data.frame(
  cid  = as.character(result$cids),
  name = as.character(result$name),
  stringsAsFactors = FALSE
)
cid_df <- cid_df[!is.na(cid_df$cid) & nzchar(cid_df$cid), ]
cid_df <- cid_df[!(cid_df$cid %in% already_done), ]
if (!nrow(cid_df)) {
  message("Nothing to do; all CIDs already processed.")
}

# ---- Small helpers ----
append_csv_safely <- function(dt, path) {
  lf <- filelock::lock(lockfile, timeout = 60000)
  on.exit(filelock::unlock(lf), add = TRUE)

  write_with_header <- !file.exists(path) || file.info(path)$size == 0
  data.table::fwrite(
    dt, path,
    append = !write_with_header,
    col.names = write_with_header
  )
}


fetch_synonyms <- function(cid) {
  url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/TXT", cid)
  syn <- tryCatch(unique(readLines(url, warn = FALSE)), error = function(e) character(0))
  syn <- unique(trimws(syn[nzchar(syn)]))
  if (!length(syn)) return(NULL)
  utc_now <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  data.table(
    synonym = syn,
    pubchem_cid = rep(cid, length(syn)),
    source = "PubChem",
    version = rep(utc_now, length(syn))
  )
}

# ---- Per-CID worker ----
process_one_cid <- function(cid, mapped_name) {
  # Convert to numeric for AnnotationGx
  cid_num <- suppressWarnings(as.numeric(cid))
  if (is.na(cid_num)) {
    message(sprintf("[CID %s] skipped: cannot coerce to numeric", cid))
    return(invisible(NULL))
  }

  # 1) properties
  prop <- tryCatch(
    AnnotationGx::mapCID2Properties(ids = cid_num, properties = properties),
    error = function(e) { message(sprintf("[CID %s] mapCID2Properties error: %s", cid, e$message)); NULL }
  )
  if (is.null(prop) || !nrow(prop)) return(invisible(NULL))

  prop$name <- mapped_name

  # 2) ChEMBL ID
  chembl <- tryCatch(
    AnnotationGx::annotatePubchemCompound(cid_num, "ChEMBL ID"),
    error = function(e) NA_character_
  )
  prop$ChEMBL_ID <- chembl


  # 3) rename columns to your schema
  present <- intersect(names(rename_map), names(prop))
  names(prop)[match(present, names(prop))] <- rename_map[present]

  # Ensure CID is named "cid" (after rename)
  if (!"cid" %in% names(prop) && "CID" %in% names(prop)) {
    setnames(prop, "CID", "cid")
  }

  # Normalize to data.table, single row
  prop_dt <- as.data.table(prop)

  # 4) synonyms
  syn_dt <- fetch_synonyms(cid)
  if (!is.null(syn_dt)) {
    # Deduplicate (cid, synonym) ignoring case
    syn_dt[, key := paste0(tolower(synonym), "::", pubchem_cid)]
    syn_dt <- syn_dt[!duplicated(key)][, key := NULL]
  }

  # 5) append both immediately (thread-safe)
  append_csv_safely(prop_dt, prop_out)
  if (!is.null(syn_dt) && nrow(syn_dt)) {
    append_csv_safely(syn_dt, syn_out)
  }

  # be gentle to APIs (even in parallel)
  Sys.sleep(0.1)
  invisible(NULL)
}

# ---- Parallel run (set workers to taste; 3–6 is a good, polite range) ----
if (nrow(cid_df)) {
  plan(multisession, workers = max(1, min(6, parallel::detectCores() - 1)))
  # Use future_lapply to parallelize over rows
  invisible(future_lapply(
    seq_len(nrow(cid_df)),
    function(i) process_one_cid(cid_df$cid[i], cid_df$name[i]),
    future.seed = TRUE
  ))
  plan(sequential)
}