library(AnnotationGx)
library(jsonlite)
library(data.table)
library(future.apply)
library(httr)

`%||%` <- function(a, b) if (!is.null(a)) a else b
printf <- function(...) { cat(sprintf(...), "\n"); flush.console() }
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
tic <- function() proc.time()[["elapsed"]]

# Define input/output directories/filenames
in_csv   <- "output_data/pset_drugs/pset_drugs_test.csv"
out_dir  <- "output_data/pset_drugs/complete"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
prop_out <- file.path(out_dir, "pset_drugs_out.csv")
syn_out  <- file.path(out_dir, "pset_drugs_synonyms.csv")
bio_out  <- file.path(out_dir, "pset_drugs_bioassays.csv")

# Load allowed AIDs (Homo sapiens assays)
aids_file <- "homosapien_aids.txt"
allowed_aids <- character()
aids_dt <- tryCatch(
    fread(aids_file),
    error = function(e) {
        printf("[%s] WARNING: could not read AIDs file '%s': %s",
               ts(), aids_file, e$message)
        NULL
    }
)
if (!is.null(aids_dt) && "aid" %in% names(aids_dt)) {
    allowed_aids <- unique(as.character(trimws(aids_dt$aid)))
    allowed_aids <- allowed_aids[nzchar(allowed_aids)]
    printf("[%s] Loaded %d allowed AIDs from %s",
           ts(), length(allowed_aids), aids_file)
} else {
    printf("[%s] WARNING: AIDs file missing or no 'aid' column; no AID filtering will be applied.",
           ts())
    allowed_aids <- character()
}

printf("[%s] Starting run", ts())
run_start <- tic()

# Read in drug names
printf("[%s] Loading input CSV: %s", ts(), in_csv)
drugs <- read.csv(in_csv, stringsAsFactors = FALSE)
if (!"drug" %in% names(drugs)) {
    stop(sprintf("Input must contain a 'drug' column. Columns found: %s",
                 paste(names(drugs), collapse = ", ")))
}
drug_names <- unique(na.omit(trimws(as.character(drugs$drug))))
printf("[%s] Loaded %d unique names", ts(), length(drug_names))

# Name mapping to CID, 10 names at a time with 5s delay
safe_map_names_to_cids <- function(names_vec) {
    names_vec <- unique(trimws(na.omit(as.character(names_vec))))
    names_vec <- names_vec[nzchar(names_vec)]
    if (!length(names_vec)) {
        return(data.frame(name = character(), cid = character(), stringsAsFactors = FALSE))
    }

    printf("[%s] Starting mapCompound2CID in mini-batches of 5", ts())

    # split names into groups (current code uses 10; change to 5 if you want strict size 5)
    mini_chunks <- split(seq_along(names_vec),
                         ceiling(seq_along(names_vec) / 10L))

    all_res <- vector("list", length(mini_chunks))
    mini_i  <- 0L

    for (mc in mini_chunks) {
        mini_i <- mini_i + 1L
        sub_names <- names_vec[mc]

        printf("[%s]   mapCompound2CID mini-batch %d/%d | names: %d",
               ts(), mini_i, length(mini_chunks), length(sub_names))

        agx <- tryCatch(
            AnnotationGx::mapCompound2CID(sub_names, first = TRUE),
            error = function(e) {
                printf("[%s]   mapCompound2CID error (mini-batch %d): %s",
                       ts(), mini_i, e$message)
                NULL
            }
        )

        if (!is.null(agx) && nrow(agx)) {
            stopifnot(all(c("name", "cids") %in% names(agx)))
            df <- data.frame(
                name = as.character(agx$name),
                cid  = as.character(agx$cids),
                stringsAsFactors = FALSE
            )
            df <- df[!is.na(df$cid) & nzchar(df$cid), , drop = FALSE]
            if (nrow(df)) {
                all_res[[mini_i]] <- df
            }
        } else {
            printf("[%s]   No CIDs returned in this mini-batch", ts())
        }

        # sleep between mini-batches to be gentle on PubChem
        Sys.sleep(5)
    }

    out <- do.call(rbind, all_res)
    if (is.null(out)) {
        out <- data.frame(name = character(), cid = character(), stringsAsFactors = FALSE)
    }

    out <- unique(out)
    printf("[%s] Finished mapCompound2CID: %d mapped name(s)", ts(), nrow(out))
    out
}

map_t0 <- tic()
cid_map <- safe_map_names_to_cids(drug_names)
printf("[%s] Name->CID mapping finished in %.1fs", ts(), tic() - map_t0)

cid_map <- cid_map[!is.na(cid_map$cid) & nzchar(cid_map$cid), , drop = FALSE]
stopifnot(nrow(cid_map) > 0)
printf("[%s] Mapped %d/%d names to PubChem CIDs.", ts(), nrow(cid_map), length(drug_names))

# Pubchem properties
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

# Renaming
rename_map <- c(
    CID = "cid", Title = "title", name = "mapped_name", ChEMBL_ID = "molecule_chembl_id", MolecularFormula = "molecular_formula",
    MolecularWeight = "molecular_weight", SMILES = "smiles", ConnectivitySMILES = "connectivity_smiles", InChI = "inchi",
    InChIKey = "inchikey", IUPACName = "iupac_name", XLogP = "xlogp", ExactMass = "exact_mass", MonoisotopicMass = "monoisotopic_mass",
    TPSA = "tpsa", Complexity = "complexity", Charge = "charge", HBondDonorCount = "h_bond_donor_count", HBondAcceptorCount = "h_bond_acceptor_count",
    RotatableBondCount = "rotatable_bond_count", HeavyAtomCount = "heavy_atom_count", IsotopeAtomCount = "isotope_atom_count",
    AtomStereoCount = "atom_stereo_count", DefinedAtomStereoCount = "defined_atom_stereo_count", UndefinedAtomStereoCount = "undefined_atom_stereo_count",
    BondStereoCount = "bond_stereo_count", DefinedBondStereoCount = "defined_bond_stereo_count", UndefinedBondStereoCount = "undefined_bond_stereo_count",
    CovalentUnitCount = "covalent_unit_count", Volume3D = "volume_3d", XStericQuadrupole3D = "x_steric_quadrupole_3d",
    YStericQuadrupole3D = "y_steric_quadrupole_3d", ZStericQuadrupole3D = "z_steric_quadrupole_3d", FeatureCount3D = "feature_count_3d",
    FeatureAcceptorCount3D = "feature_acceptor_count_3d", FeatureDonorCount3D = "feature_donor_count_3d", FeatureAnionCount3D = "feature_anion_count_3d",
    FeatureCationCount3D = "feature_cation_count_3d", FeatureRingCount3D = "feature_ring_count_3d", FeatureHydrophobeCount3D = "feature_hydrophobe_count_3d",
    ConformerModelRMSD3D = "conformer_model_rmsd_3d", EffectiveRotorCount3D = "effective_rotor_count_3d", ConformerCount3D = "conformer_count_3d",
    Fingerprint2D = "fingerprint_2d", PatentCount = "patent_count", PatentFamilyCount = "patent_family_count", LiteratureCount = "literature_count",
    AnnotationTypes = "annotation_types", AnnotationTypeCount = "annotation_type_count"
)

.syn_handle <- httr::handle("https://pubchem.ncbi.nlm.nih.gov")

retry <- function(fun, times = 2, base_wait = 0.25) {
    for (i in seq_len(times)) {
        res <- try(fun(), silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
        if (i < times) Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1, 0, 0.1))
    }
    stop("retry failed")
}

fetch_synonyms_batch <- function(cids_chr) {
    if (!length(cids_chr)) return(data.table())
    path <- sprintf("/rest/pug/compound/cid/%s/synonyms/JSON", paste(cids_chr, collapse = ","))
    req_fun <- function() {
        resp <- httr::GET(
            handle = .syn_handle, path = path,
            httr::timeout(30),
            httr::user_agent("AnnotationDB/pset_drugs (synonyms_batch)")
        )
        sc <- httr::status_code(resp)
        if (sc >= 500 || sc == 429) stop("server busy")
        resp
    }
    res <- tryCatch(retry(req_fun), error = function(e) NULL)
    out <- data.table()

    if (!is.null(res) && httr::status_code(res) == 200) {
        txt <- httr::content(res, as = "text", encoding = "UTF-8")
        j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)
        if (!inherits(j, "try-error")) {
            info <- j$InformationList$Information
            if (!is.null(info) && length(info)) {
                out <- rbindlist(lapply(info, function(x) {
                    if (is.null(x$Synonym)) return(NULL)
                    data.table(
                        synonym     = unique(trimws(x$Synonym)),
                        pubchem_cid = as.character(x$CID),
                        source      = "PubChem",
                        version     = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
                    )
                }), fill = TRUE)
            }
        }
    }

    have_cid <- unique(out$pubchem_cid)
    missing  <- setdiff(as.character(cids_chr), have_cid)
    if (length(missing)) {
        txt_list <- lapply(missing, function(cid) {
            url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/TXT", cid)
            syn <- tryCatch(unique(readLines(url, warn = FALSE)), error = function(e) character(0))
            syn <- unique(trimws(syn[nzchar(syn)]))
            if (!length(syn)) return(NULL)
            data.table(
                synonym     = syn,
                pubchem_cid = rep(as.character(cid), length(syn)),
                source      = "PubChem",
                version     = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
            )
        })
        out <- rbindlist(c(list(out), txt_list), fill = TRUE)
    }

    if (nrow(out)) {
        out[, key := paste0(tolower(synonym), "::", pubchem_cid)]
        out <- out[!duplicated(key)][, key := NULL]
    }
    out
}

# ChEMBL max_phase (for reporting)
.fetch_chembl_phase_cache <- new.env(parent = emptyenv())

fetch_chembl_phase <- function(chembl_id) {
    if (is.na(chembl_id) || !nzchar(chembl_id)) return(NA_integer_)
    if (exists(chembl_id, envir = .fetch_chembl_phase_cache, inherits = FALSE)) {
        return(get(chembl_id, envir = .fetch_chembl_phase_cache))
    }
    url <- sprintf("https://www.ebi.ac.uk/chembl/api/data/molecule/%s?format=json", chembl_id)
    resp <- try(httr::GET(
        url, httr::timeout(20),
        httr::user_agent("AnnotationDB/pset_drugs (fetch_chembl_phase)")
    ), silent = TRUE)
    if (inherits(resp, "try-error") || httr::status_code(resp) != 200) {
        assign(chembl_id, NA_integer_, envir = .fetch_chembl_phase_cache); return(NA_integer_)
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)
    if (inherits(j, "try-error") || is.null(j$max_phase)) {
        assign(chembl_id, NA_integer_, envir = .fetch_chembl_phase_cache); return(NA_integer_)
    }
    phase <- suppressWarnings(as.integer(j$max_phase))
    assign(chembl_id, phase, envir = .fetch_chembl_phase_cache)
    phase
}

# Bio Assay Description lookup
fetch_assay_description <- function(aid, tries = 2, base_wait = 0.25) {
    # Returns one row with: aid, protein_target_name, protein_accession, source_name, source_id
    aid_chr <- as.character(aid)
    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/description/JSON", aid_chr)

    for (i in seq_len(tries)) {
        resp <- try(httr::GET(
            url, httr::timeout(20),
            httr::user_agent("AnnotationDB/pset_drugs (assay_description)")
        ), silent = TRUE)
        if (inherits(resp, "try-error")) {
            if (i < tries) { Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.1)); next }
            return(data.table(
                aid                 = aid_chr,
                protein_target_name = NA_character_,
                protein_accession   = NA_character_,
                source_name         = NA_character_,
                source_id           = NA_character_
            ))
        }

        if (httr::status_code(resp) == 200) {
            txt <- httr::content(resp, as = "text", encoding = "UTF-8")
            j <- try(jsonlite::fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
            if (inherits(j, "try-error")) break

            # Navigate: PC_AssayContainer[[1]]$assay$descr
            descr <- try(j$PC_AssayContainer[[1]]$assay$descr, silent = TRUE)
            if (inherits(descr, "try-error") || is.null(descr)) break

            # Source name + id (robust to differing shapes)
            src <- descr$aid_source
            source_name <- try(src$db$name, silent = TRUE)
            if (inherits(source_name, "try-error") || is.null(source_name)) source_name <- NA_character_
            source_id <- NA_character_
            if (!is.null(src)) {
                # observed variants: source_id$str, source_id$id, or scalar
                sid1 <- try(src$source_id$str, silent = TRUE)
                sid2 <- try(src$source_id$id,  silent = TRUE)
                sid3 <- try(src$source_id,     silent = TRUE)
                cand <- c(
                    if (!inherits(sid1, "try-error")) sid1 else NULL,
                    if (!inherits(sid2, "try-error")) sid2 else NULL,
                    if (!inherits(sid3, "try-error")) sid3 else NULL
                )
                if (length(cand)) source_id <- as.character(cand[[1]])
                if (is.null(source_id)) source_id <- NA_character_
            }

            # Target protein fields
            target <- descr$target
            protein_accession   <- NA_character_
            protein_target_name <- NA_character_
            if (!is.null(target) && length(target) >= 1) {
                protein_accession   <- try(target[[1]]$protein_accession, silent = TRUE)
                if (inherits(protein_accession, "try-error") || is.null(protein_accession)) protein_accession <- NA_character_
                protein_target_name <- try(target[[1]]$name, silent = TRUE)
                if (inherits(protein_target_name, "try-error") || is.null(protein_target_name)) protein_target_name <- NA_character_
            }

            return(data.table(
                aid                 = aid_chr,
                protein_target_name = as.character(protein_target_name %||% NA_character_),
                protein_accession   = as.character(protein_accession %||% NA_character_),
                source_name         = as.character(source_name %||% NA_character_),
                source_id           = as.character(source_id %||% NA_character_)
            ))
        }

        sc <- httr::status_code(resp)
        if (sc %in% c(429) || sc >= 500) {
            if (i < tries) { Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.2)); next }
        }
        break
    }

    data.table(
        aid                 = aid_chr,
        protein_target_name = NA_character_,
        protein_accession   = NA_character_,
        source_name         = NA_character_,
        source_id           = NA_character_
    )
}

fetch_assay_descriptions_for_aids <- function(aids_chr) {
    if (!length(aids_chr)) {
        return(data.table(
            aid                 = character(),
            protein_target_name = character(),
            protein_accession   = character(),
            source_name         = character(),
            source_id           = character()
        ))
    }
    out <- lapply(aids_chr, fetch_assay_description)
    rbindlist(out, fill = TRUE)
}

fetch_bioassay_summary_one <- function(cid, tries = 2, base_wait = 0.25) {
    cid_chr <- as.character(cid)
    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON", cid_chr)

    for (i in seq_len(tries)) {
        resp <- try(httr::GET(
            url, httr::timeout(30),
            httr::user_agent("AnnotationDB/pset_drugs (assaysummary)")
        ), silent = TRUE)

        if (inherits(resp, "try-error")) {
            if (i < tries) { Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.1)); next }
            return(data.table())
        }

        if (httr::status_code(resp) == 200) {
            txt <- httr::content(resp, as = "text", encoding = "UTF-8")
            j <- try(jsonlite::fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
            if (inherits(j, "try-error")) return(data.table())

            tab <- j$Table
            if (is.null(tab)) return(data.table())

            col_names <- try(unlist(tab$Columns$Column, use.names = FALSE), silent = TRUE)
            if (inherits(col_names, "try-error") || is.null(col_names)) return(data.table())

            idx <- function(needles) which(tolower(col_names) %in% tolower(needles))[1] %||% NA_integer_
            ix_aid    <- idx(c("AID","Bioassay AID"))
            ix_tname  <- idx(c("Target Name","Bioassay Target Name"))
            ix_tacc   <- idx(c("Target Accession"))
            ix_gene   <- idx(c("Target GeneID","GeneID"))
            ix_aname  <- idx(c("Assay Name"))

            rows <- tab$Row
            if (is.null(rows) || !length(rows)) return(data.table())

            # ⛔ OLD: rows <- rows[1:50]
            # ✅ NEW: use ALL rows without a cap
            # (do nothing)

            out <- lapply(rows, function(r) {
                vals <- try(unlist(r$Cell, use.names = FALSE), silent = TRUE)
                if (inherits(vals, "try-error")) return(NULL)
                pick <- function(ix) if (is.na(ix) || ix < 1 || ix > length(vals)) NA_character_ else as.character(vals[ix])

                aid <- pick(ix_aid)

                tname <- pick(ix_tname)
                if (is.na(tname) || !nzchar(tname)) tname <- pick(ix_tacc)
                if (is.na(tname) || !nzchar(tname)) {
                    g <- pick(ix_gene)
                    if (!is.na(g) && nzchar(g)) tname <- paste0("GeneID:", g)
                }
                if (is.na(tname) || !nzchar(tname)) tname <- pick(ix_aname)

                data.table(
                    pubchem_cid          = cid_chr,
                    bioassay_aid         = aid,
                    bioassay_target_name = tname
                )
            })

            dt <- rbindlist(out, fill = TRUE)
            if (!nrow(dt)) return(data.table())
            dt <- dt[!is.na(bioassay_aid) & nzchar(bioassay_aid)]
            dt <- unique(dt)

            # augment with assay descriptions
            aids <- unique(dt$bioassay_aid)
            desc_dt <- fetch_assay_descriptions_for_aids(aids)
            setnames(desc_dt, "aid", "bioassay_aid")

            dt <- merge(dt, desc_dt, by = "bioassay_aid", all.x = TRUE, sort = FALSE)

            return(dt)
        }

        sc <- httr::status_code(resp)
        if (sc %in% c(429) || sc >= 500) {
            if (i < tries) { Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.2)); next }
        }

        return(data.table())
    }
    data.table()
}


fetch_bioassay_for_cids <- function(cids_chr) {
    if (!length(cids_chr)) return(data.table())
    rbindlist(lapply(cids_chr, fetch_bioassay_summary_one), fill = TRUE)
}

# Chunking
CHUNK_SIZE <- 100L
n_total <- nrow(cid_map)
chunks <- split(seq_len(n_total), ceiling(seq_along(seq_len(n_total)) / CHUNK_SIZE))
printf("[%s] Processing %d CIDs across %d chunk(s) (chunk size %d)", ts(), n_total, length(chunks), CHUNK_SIZE)

props_all <- list()
syn_all   <- list()
bio_all   <- list()

chunk_i <- 0L
for (idx in chunks) {
    chunk_i <- chunk_i + 1L
    chunk <- cid_map[idx, , drop = FALSE]
    printf("[%s] Chunk %d/%d | CIDs: %d", ts(), chunk_i, length(chunks), nrow(chunk))

    # --- Properties ---
    prop_rows <- 0L
    local_props_list <- list()

    cids_chr_chunk <- as.character(chunk$cid)
    mini_chunks <- split(seq_along(cids_chr_chunk),
                         ceiling(seq_along(cids_chr_chunk) / 10L))

    mini_i <- 0L
    for (mc in mini_chunks) {
        mini_i <- mini_i + 1L
        sub_cids_chr <- cids_chr_chunk[mc]
        sub_cids_num <- suppressWarnings(as.numeric(sub_cids_chr))

        names_map <- setNames(chunk$name[match(sub_cids_chr, chunk$cid)],
                              nm = as.character(sub_cids_chr))

        printf("[%s]   Properties mini-batch %d/%d | CIDs: %d",
               ts(), mini_i, length(mini_chunks), length(sub_cids_chr))

        prop <- tryCatch(
            AnnotationGx::mapCID2Properties(ids = sub_cids_num, properties = properties),
            error = function(e) {
                printf("[%s]   mapCID2Properties error (mini-batch): %s", ts(), e$message)
                NULL
            }
        )

        if (!is.null(prop) && nrow(prop)) {
            prop$CID  <- as.character(prop$CID)
            prop$name <- names_map[prop$CID]

            chembl <- tryCatch(
                AnnotationGx::annotatePubchemCompound(as.numeric(prop$CID), "ChEMBL ID"),
                error = function(e) {
                    printf("[%s]   annotatePubchemCompound error (mini-batch): %s", ts(), e$message)
                    NA_character_
                }
            )

            if (is.atomic(chembl) && length(chembl) %in% c(1, nrow(prop))) {
                prop$ChEMBL_ID <- chembl
            } else {
                prop$ChEMBL_ID <- vapply(prop$CID, function(x)
                    tryCatch(AnnotationGx::annotatePubchemCompound(as.numeric(x), "ChEMBL ID"),
                             error = function(e) NA_character_), character(1))
            }

            # ChEMBL max_phase
            prop$chembl_max_phase <- vapply(as.character(prop$ChEMBL_ID), fetch_chembl_phase, integer(1))

            local_props_list[[length(local_props_list) + 1L]] <- as.data.table(prop)
        } else {
            printf("[%s]   No properties returned for this mini-batch", ts())
        }

        # Respectful pause between properties mini-batches
        Sys.sleep(5)
    }

    if (length(local_props_list)) {
        prop_dt <- rbindlist(local_props_list, fill = TRUE)

        # Rename to target schema
        present <- intersect(names(rename_map), names(prop_dt))
        names(prop_dt)[match(present, names(prop_dt))] <- rename_map[present]
        if (!"cid" %in% names(prop_dt) && "CID" %in% names(prop_dt)) {
            data.table::setnames(prop_dt, "CID", "cid")
        }

        # drug_like derived from annotation_types
        prop_dt[, drug_like := !is.na(annotation_types) &
                            grepl("(^|\\|)\\s*Drug and Medication Information\\s*(\\||$)",
                                  annotation_types, ignore.case = TRUE)]

        prop_rows <- nrow(prop_dt)
        props_all[[length(props_all) + 1L]] <- prop_dt
        printf("[%s]   Properties rows (chunk): %d", ts(), prop_rows)
    } else {
        printf("[%s]   No properties returned for this chunk", ts())
    }

    # --- Synonyms  ---
    syn_rows <- 0L
    local_syn_list <- list()
    syn_cids_chunk <- unique(chunk$cid)
    syn_mini_chunks <- split(seq_along(syn_cids_chunk),
                             ceiling(seq_along(syn_cids_chunk) / 10L))

    syn_i <- 0L
    for (sc_idx in syn_mini_chunks) {
        syn_i <- syn_i + 1L
        sub_syn_cids <- syn_cids_chunk[sc_idx]

        printf("[%s]   Synonyms mini-batch %d/%d | CIDs: %d",
               ts(), syn_i, length(syn_mini_chunks), length(sub_syn_cids))

        syn_dt_sub <- fetch_synonyms_batch(sub_syn_cids)
        if (nrow(syn_dt_sub)) {
            local_syn_list[[length(local_syn_list) + 1L]] <- syn_dt_sub
            syn_rows <- syn_rows + nrow(syn_dt_sub)
        }

        Sys.sleep(5)
    }

    if (length(local_syn_list)) {
        syn_dt <- rbindlist(local_syn_list, fill = TRUE)
        syn_all[[length(syn_all) + 1L]] <- syn_dt
    } else {
        syn_dt <- data.table()
    }
    printf("[%s]   Synonyms rows: %d", ts(), syn_rows)

    # --- Bioassays  ---
    bio_rows <- 0L
    local_bio_list <- list()
    bio_cids_chunk <- unique(chunk$cid)
    bio_mini_chunks <- split(seq_along(bio_cids_chunk),
                             ceiling(seq_along(bio_cids_chunk) / 10L))

    bio_i <- 0L
    for (bc_idx in bio_mini_chunks) {
        bio_i <- bio_i + 1L
        sub_bio_cids <- bio_cids_chunk[bc_idx]

        printf("[%s]   Bioassays mini-batch %d/%d | CIDs: %d",
               ts(), bio_i, length(bio_mini_chunks), length(sub_bio_cids))

        bio_dt_sub <- fetch_bioassay_for_cids(sub_bio_cids)
        if (nrow(bio_dt_sub)) {
            local_bio_list[[length(local_bio_list) + 1L]] <- bio_dt_sub
            bio_rows <- bio_rows + nrow(bio_dt_sub)
        }

        Sys.sleep(5)
    }

    if (length(local_bio_list)) {
        bio_dt <- rbindlist(local_bio_list, fill = TRUE)
        bio_all[[length(bio_all) + 1L]] <- bio_dt
    } else {
        bio_dt <- data.table()
    }
    printf("[%s]   Bioassays rows: %d", ts(), bio_rows)

    Sys.sleep(0.15)
}

# Combine
props_all <- if (length(props_all)) rbindlist(props_all, fill = TRUE) else data.table()
syn_all   <- if (length(syn_all))   rbindlist(syn_all,   fill = TRUE) else data.table()
bio_all   <- if (length(bio_all))   rbindlist(bio_all,   fill = TRUE) else data.table()
printf("[%s] Combined totals — props: %d | syns: %d | bio: %d", ts(), nrow(props_all), nrow(syn_all), nrow(bio_all))

# Filter bioassays to only those whose AID is in the allowed list
if (length(allowed_aids) && nrow(bio_all)) {
    before_n <- nrow(bio_all)
    bio_all <- bio_all[bioassay_aid %in% allowed_aids]
    printf("[%s] Filtered bio by AIDs list: %d -> %d rows",
           ts(), before_n, nrow(bio_all))
} else {
    printf("[%s] No AID filtering applied (either no allowed_aids or no bio_all rows).", ts())
}

# Deduplicate conservatively
if (nrow(props_all)) {
    props_all <- unique(props_all)
    printf("[%s] Dedup props -> %d", ts(), nrow(props_all))
}
if (nrow(syn_all)) {
    syn_all[, key := paste0(tolower(synonym), "::", pubchem_cid)]
    syn_all <- syn_all[!duplicated(key)][, key := NULL]
    printf("[%s] Dedup syns -> %d", ts(), nrow(syn_all))
}
if (nrow(bio_all)) {
    bio_all[, key := paste0(pubchem_cid, "::", bioassay_aid)]
    bio_all <- bio_all[!duplicated(key)][, key := NULL]
    printf("[%s] Dedup bio -> %d", ts(), nrow(bio_all))
}

# Writing Outputs
printf("[%s] Writing outputs ...", ts())
fwrite(props_all, prop_out)
fwrite(syn_all,   syn_out)
fwrite(bio_all,   bio_out)
printf("[%s] Wrote: %s (%d rows)", ts(), prop_out, nrow(props_all))
printf("[%s] Wrote: %s (%d rows)", ts(), syn_out,  nrow(syn_all))
printf("[%s] Wrote: %s (%d rows)", ts(), bio_out,  nrow(bio_all))

printf("[%s] Done in %.1fs", ts(), tic() - run_start)
