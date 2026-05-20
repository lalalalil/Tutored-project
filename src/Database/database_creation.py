import sqlite3
import json
import re
import requests
import pathlib
import pandas as pd

# =====================================================================
# Schema of the database
# =====================================================================

def create_database(db_name="infectio_git.db"):
    conn = sqlite3.connect(db_name)
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()
    cur.executescript("""

    -- ══════════════════════════════════════════
    -- MAIN TABLES
    -- ══════════════════════════════════════════

    CREATE TABLE IF NOT EXISTS Model (
        doi         TEXT PRIMARY KEY,  -- DOI 
        title       TEXT,
        formalism   TEXT,              -- ODE, boolean, ABM, SBML, stochastic...
        scale       TEXT,              -- molecular / cellular / tissue / host / epidemiological
        year        INTEGER,
        software    TEXT,              -- COPASI, MATLAB, Python, R...
        source_db   TEXT               -- 'biomodels', 'zenodo', 'ncbi'
    );

    CREATE TABLE IF NOT EXISTS Disease (
        doid_mondo  TEXT PRIMARY KEY,  
        name        TEXT
    );

    CREATE TABLE IF NOT EXISTS Pathogen (
        taxonomic_id    TEXT PRIMARY KEY,  -- NCBI Taxonomy ID
        species         TEXT,
        class           TEXT,              -- virus / bacterium / parasite / fungus
        life_cycle_stage TEXT              -- replicative / latent / intracellular...
    );

    CREATE TABLE IF NOT EXISTS Host (
        taxonomic_id    TEXT PRIMARY KEY,  -- NCBI Taxonomy ID
        species         TEXT,
        in_vitro_vivo   TEXT               -- 'in vitro' / 'in vivo' / 'in silico'
    );

    CREATE TABLE IF NOT EXISTS Tissue (
        ontology_id     TEXT PRIMARY KEY,  
        name            TEXT               -- lung, blood, intestine...
    );

    CREATE TABLE IF NOT EXISTS Cell_types (
        cell_ontology   TEXT PRIMARY KEY,  
        name            TEXT               -- macrophage, T lymphocyte...
    );
                      
    CREATE TABLE IF NOT EXISTS Dataset (
        identifier  TEXT PRIMARY KEY,  -- ex: GSE12345, clinical accession...
        source      TEXT,              -- GEO / clinical / in vitro / in vivo
        description TEXT
    );

    -- ══════════════════════════════════════════
    -- JUNCTION TABLES
    -- ══════════════════════════════════════════

    CREATE TABLE IF NOT EXISTS model_disease (
        model_doi       TEXT,
        disease_doid    TEXT,
        PRIMARY KEY (model_doi, disease_doid),
        FOREIGN KEY (model_doi)     REFERENCES Model(doi),
        FOREIGN KEY (disease_doid)  REFERENCES Disease(doid_mondo)
    );

    CREATE TABLE IF NOT EXISTS model_pathogen (
        model_doi           TEXT,
        pathogen_taxonomy   TEXT,
        PRIMARY KEY (model_doi, pathogen_taxonomy),
        FOREIGN KEY (model_doi)             REFERENCES Model(doi),
        FOREIGN KEY (pathogen_taxonomy)     REFERENCES Pathogen(taxonomic_id)
    );

    CREATE TABLE IF NOT EXISTS model_host (
        model_doi       TEXT,
        host_taxonomy   TEXT,
        PRIMARY KEY (model_doi, host_taxonomy),
        FOREIGN KEY (model_doi)     REFERENCES Model(doi),
        FOREIGN KEY (host_taxonomy) REFERENCES Host(taxonomic_id)
    );

    CREATE TABLE IF NOT EXISTS model_tissue (
        model_doi       TEXT,
        tissue_ontology TEXT,
        PRIMARY KEY (model_doi, tissue_ontology),
        FOREIGN KEY (model_doi)         REFERENCES Model(doi),
        FOREIGN KEY (tissue_ontology)   REFERENCES Tissue(ontology_id)
    );

    CREATE TABLE IF NOT EXISTS model_cell_types (
        model_doi           TEXT,
        cell_ontology       TEXT,
        PRIMARY KEY (model_doi, cell_ontology),
        FOREIGN KEY (model_doi)         REFERENCES Model(doi),
        FOREIGN KEY (cell_ontology)     REFERENCES Cell_types(cell_ontology)
    );

    CREATE TABLE IF NOT EXISTS model_dataset (
        model_doi           TEXT,
        dataset_identifier  TEXT,
        PRIMARY KEY (model_doi, dataset_identifier),
        FOREIGN KEY (model_doi)             REFERENCES Model(doi),
        FOREIGN KEY (dataset_identifier)    REFERENCES Dataset(identifier)
    );

    """)
    conn.commit()
    print(f"Database created : {db_name}")
    return conn

conn = create_database()

# =====================================================================
# Automatic detection of the JSON source metadata
# =====================================================================

def detect_source(data):
    # Detects the source (biomodels / zenodo / ncbi) from the JSON structure
    if isinstance(data, dict):
        # BioModels format : model_id at the root and web_metadata as a sub-key
        if "model_id" in data and "web_metadata" in data:
            return "biomodels"
        # Other BioModels format : keys directly at the root
        if "modelLevelAnnotations" in data or "publicationId" in data or "submissionId" in data:
            return "biomodels"
        if "conceptrecid" in data or (
            "links" in data and "zenodo.org" in str(data.get("links", {}))
        ):
            return "zenodo"
    if isinstance(data, list) and len(data) > 0:
        item = data[0]
        if isinstance(item, dict) and (
            "NlmUniqueID" in item or "PmcRefCount" in item or "ArticleIds" in item
        ):
            return "ncbi"
    return "unknown"

# =====================================================================
# Utilitary function
# =====================================================================

# Known taxons for Homo sapiens
HUMAN_TAXON_IDS = {"9606"}

# Keywords : pathogen class
PATHOGEN_CLASS_KEYWORDS = {
    "virus": [
        "virus", "viral", "sars", "covid", "influenza", "dengue",
        "hiv", "coronavirus", "zika", "west nile", "chikungunya",
    ]
}

# Keywords → known DOID (for inference from free text)
# IMPORTANT : only inference on the model title/name is used,
# not on the full synopsis, to avoid context false positives.
DISEASE_KEYWORDS = {
    "covid-19": ("DOID:0080600", "COVID-19"),
    "covid19": ("DOID:0080600", "COVID-19"),
    "sars-cov-2": ("DOID:0080600", "COVID-19"),
    "influenza": ("DOID:8469", "influenza"),
    "tuberculosis": ("DOID:399", "tuberculosis"),
    "dengue": ("DOID:12205", "dengue fever"),
    "hiv": ("DOID:526", "HIV infectious disease"),
    "zika": ("DOID:0060478", "Zika fever"),
    "chikungunya": ("DOID:0050012", "chikungunya"),
}

def infer_pathogen_class(name: str) -> str:
    # Gets the class of a pathogen from its name (NCBI Taxonomy)
    name_l = name.lower()
    for cls, keywords in PATHOGEN_CLASS_KEYWORDS.items():
        if any(k in name_l for k in keywords):
            return cls
    return ""


def infer_diseases_from_text(text: str) -> list:
    """
    Infers a list of diseases (doid_mondo, name) from free text.
    Deduplicates by doid_mondo.
    """
    text_l = text.lower()
    found, seen = [], set()
    for keyword, (doid, name) in DISEASE_KEYWORDS.items():
        if keyword in text_l and doid not in seen:
            seen.add(doid)
            found.append({"doid_mondo": doid, "name": name})
    return found


def detect_software_from_biomodels(data: dict) -> str:
    # Deduces the software used from the BioModels additional files
    for f in data.get("files", {}).get("additional", []):
        name = f.get("name", "").lower()
        if ".cps" in name or "copasi" in name:
            return "COPASI"
        if "matlab" in name:
            return "MATLAB"
        if "octave" in name:
            return "Octave"
        if ".ode" in name:
            return "XPP/ODE"
    fmt = data.get("format", {}).get("name", "")
    return fmt if fmt else ""

# =====================================================================
# Parsers by source
# =====================================================================

def _empty_result():
    return {
        "model": {},
        "diseases": [],
        "pathogens": [],
        "hosts": [],
        "tissues": [],
        "cell_types": [],
        "datasets": [],
    }

# ────────────────────────────────────────────
# 1. BIOMODELS PARSER
# ────────────────────────────────────────────
def parse_biomodels(data: dict) -> dict:
    """
    Parses a BioModels metadata file.
    Handles two formats :
      - model_id + year at the root, everything else under web_metadata
      - keys directly at the root (publicationId, modelLevelAnnotations…)
    """
    result = _empty_result()

    # Format detection
    if "web_metadata" in data:
        wm = data["web_metadata"]
        model_id  = data.get("model_id", "")
        year_root = data.get("year")  # reliable year from the root
    else:
        wm = data
        model_id  = None
        year_root = None

    pub = wm.get("publication", {})

    # Main identifier
    identifier = (
        model_id
        or wm.get("publicationId")
        or wm.get("submissionId", "")
    )

    # Formalism
    fmt_name   = wm.get("format", {}).get("name", "")
    fmt_ver    = wm.get("format", {}).get("version", "")
    approach   = wm.get("modellingApproach", {}).get("name", "")
    formalism  = f"{fmt_name} {fmt_ver}".strip() if fmt_name else approach

    # Year : priority to year_root (root), otherwise publication
    year = year_root or pub.get("year")

    result["model"] = {
        "doi":       identifier,
        "title":     wm.get("name", ""),
        "formalism": formalism,
        "scale":     _infer_scale_biomodels(wm),
        "year":      year,
        "software":  detect_software_from_biomodels(wm),
        "source_db": "biomodels",
    }

    # Annotations 
    for annot in wm.get("modelLevelAnnotations", []):
        qualifier  = annot.get("qualifier", "")
        resource   = annot.get("resource", "")
        accession  = str(annot.get("accession", ""))
        name       = annot.get("name", "").strip()
        uri        = annot.get("uri", "")

        # Disease : several possible qualifiers depending on the format 
        is_disease_resource = (
            "doid" in resource.lower()
            or "disease" in resource.lower()
            or "doid" in uri.lower()
            or accession.startswith("DOID:")
        )
        if is_disease_resource and qualifier in (
            "bqbiol:isVersionOf", "bqbiol:is", "bqbiol:isDescribedBy",
            "bqmodel:is", "bqbiol:hasProperty",
        ):
            doid = accession if accession.startswith("DOID:") else f"DOID:{accession}"
            if not any(d["doid_mondo"] == doid for d in result["diseases"]):
                result["diseases"].append({"doid_mondo": doid, "name": name})

        # Taxon : split host (Homo sapiens) and pathogen
        if qualifier == "bqbiol:hasTaxon":
            if accession in HUMAN_TAXON_IDS or "homo sapiens" in name.lower():
                result["hosts"].append({
                    "taxonomic_id": accession,
                    "species": name,
                    "in_vitro_vivo": "",
                })
            elif accession:  # any other non-human taxon : pathogen
                result["pathogens"].append({
                    "taxonomic_id": accession,
                    "species": name,
                    "class": infer_pathogen_class(name),
                    "life_cycle_stage": "",
                })

        # Cell types : CL ontology
        if "CL:" in accession or "cell" in resource.lower():
            if qualifier in ("bqbiol:is", "bqbiol:isVersionOf", "bqbiol:hasProperty"):
                result["cell_types"].append({
                    "cell_ontology": accession,
                    "name": name,
                })

        # Tissues : UBERON ontology
        if "UBERON:" in accession or "uberon" in uri.lower():
            if qualifier in ("bqbiol:is", "bqbiol:isVersionOf", "bqbiol:hasProperty"):
                result["tissues"].append({
                    "ontology_id": accession,
                    "name": name,
                })

    # Inference from free text if annotations are insufficient
    # Short text (title only) :
    text_free = wm.get("name", "") + " " + pub.get("title", "")

    # Diseases from text if none found via annotations
    # Only name + title used (not the synopsis) to avoid false positives
    if not result["diseases"]:
        text_short = wm.get("name", "") + " " + pub.get("title", "")
        result["diseases"] = infer_diseases_from_text(text_short)

    # Pathogens from text if none found via annotations
    # Uses name + title + synopsis, but avoids duplicates with annotations
    if not result["pathogens"]:
        result["pathogens"].extend(_infer_pathogens_from_text(text_free))
    else:
        # Even if we have pathogens via annotations, we complete from the title only
        text_short_patho = wm.get("name", "") + " " + pub.get("title", "")
        existing_ids = {p["taxonomic_id"] for p in result["pathogens"]}
        for p in _infer_pathogens_from_text(text_short_patho):
            if p["taxonomic_id"] not in existing_ids:
                result["pathogens"].append(p)
                existing_ids.add(p["taxonomic_id"])

    return result


def _infer_scale_biomodels(data: dict) -> str:
    # Infers the scale from the model name or the modelling approach
    approach = data.get("modellingApproach", {}).get("name", "").lower()
    name     = data.get("name", "").lower()
    text     = approach + " " + name
    if "population" in text or "epidemi" in text or "seir" in text or "sir" in text:
        return "epidemiological"
    if "intracellular" in text or "molecular" in text or "genome" in text or "metaboli" in text:
        return "molecular"
    if "tissue" in text or "organ" in text:
        return "tissue"
    if "cellular" in text or "cell" in text:
        return "cellular"
    return ""


# Known species dictionary → (taxon_id, class)
KNOWN_PATHOGENS_TEXT = [
    ("sars-cov-2",            ("2697049", "Severe acute respiratory syndrome coronavirus 2", "virus")),
    ("sars-cov",              ("2697049", "Severe acute respiratory syndrome coronavirus 2", "virus")),
    ("covid",                 ("2697049", "Severe acute respiratory syndrome coronavirus 2", "virus")),
    ("coronavirus",           ("2697049", "Severe acute respiratory syndrome coronavirus 2", "virus")),
    ("influenza",             ("11520",   "Influenza A virus", "virus")),
    ("mycobacterium tuberculosis", ("1773", "Mycobacterium tuberculosis", "bacterie")),
    ("tuberculosis",          ("1773",   "Mycobacterium tuberculosis", "bacterie")),
    ("dengue",                ("12637",  "Dengue virus", "virus")),
    ("hiv",                   ("11676",  "Human immunodeficiency virus", "virus")),
    ("chikungunya",           ("37124",  "Chikungunya virus", "virus")),
    ("zika",                  ("64320",  "Zika virus", "virus")),
]


def _infer_pathogens_from_text(text: str) -> list:
    # Infers pathogens from free text via a keyword list
    text_l = text.lower()
    found, seen_ids = [], set()
    for keyword, (taxid, species, cls) in KNOWN_PATHOGENS_TEXT:
        if keyword in text_l and taxid not in seen_ids:
            seen_ids.add(taxid)
            found.append({
                "taxonomic_id": taxid,
                "species": species,
                "class": cls,
                "life_cycle_stage": "",
            })
    return found

# ────────────────────────────────────────────
# 2. ZENODO PARSER
# ────────────────────────────────────────────
def parse_zenodo(data: dict) -> dict:
    """
    Parses a Zenodo metadata file.
    Structured biological metadata is rare here ;
    disease/pathogen is inferred from the title/description.
    """
    result = _empty_result()
    meta = data.get("metadata", data)  # fields are sometimes at the root

    doi   = data.get("doi") or meta.get("doi", "")
    title = data.get("title") or meta.get("title", "")

    # Programming language : software
    lang_list = meta.get("custom", {}).get("code:programmingLanguage", [])
    software  = ", ".join(
        lang.get("title", {}).get("en", "") for lang in lang_list
    ) if lang_list else ""

    # Publication date : year
    pub_date = meta.get("publication_date", "")
    year = int(pub_date[:4]) if pub_date and pub_date[:4].isdigit() else None

    # Resource type : formalism
    resource_type = meta.get("resource_type", {})
    formalism = resource_type.get("title", resource_type.get("type", ""))

    result["model"] = {
        "doi":       doi,
        "title":     title,
        "formalism": formalism,
        "scale":     "",  # not present in Zenodo
        "year":      year,
        "software":  software,
        "source_db": "zenodo",
    }

    # Inference from free text (title + description)
    desc = meta.get("description", "")
    text = title + " " + desc
    result["diseases"] = infer_diseases_from_text(text)

    # Linked datasets (related_identifiers)
    for rel in meta.get("related_identifiers", []):
        ident = rel.get("identifier", "")
        if ident:
            result["datasets"].append({
                "identifier":  ident,
                "source":      rel.get("resource_type", ""),
                "description": f"relation: {rel.get('relation', '')}",
            })

    return result

# ────────────────────────────────────────────
# 3. NCBI PARSER
# ────────────────────────────────────────────
def parse_ncbi(data) -> dict:
    """
    Parses an NCBI/PubMed metadata file (list format).
    Mainly publication metadata, little biological data.
    Disease/pathogen inferred from the article title.
    """
    result = _empty_result()
    item = data[0] if isinstance(data, list) else data

    doi = (
        item.get("DOI")
        or item.get("ArticleIds", {}).get("doi", "")
        or item.get("ELocationID", "").replace("doi: ", "")
    )
    title = item.get("Title", "")

    # Year
    pub_date = item.get("PubDate", "")
    year = int(pub_date[:4]) if pub_date and pub_date[:4].isdigit() else None

    result["model"] = {
        "doi":       doi,
        "title":     title,
        "formalism": "",  # not structured in NCBI
        "scale":     "",
        "year":      year,
        "software":  "",
        "source_db": "ncbi",
    }

    # Disease inference from title
    result["diseases"] = infer_diseases_from_text(title)

    return result

# ────────────────────────────────────────────
# Main dispatcher
# ────────────────────────────────────────────
def parse_metadata(data, path="") -> dict:
    # Detects the source and calls the right parser
    source = detect_source(data)
    if source == "biomodels":
        return parse_biomodels(data)
    elif source == "zenodo":
        return parse_zenodo(data)
    elif source == "ncbi":
        return parse_ncbi(data)
    else:
        print(f"Unknown source, file ignored. ({path})")
        return _empty_result()

# =====================================================================
# Insertion in the database
# =====================================================================

def insert_parsed(conn: sqlite3.Connection, parsed: dict, path: str = "", verbose: bool = True):
    # Inserts scanned data without duplicates
    m = parsed["model"]
    doi = m.get("doi", "").strip()

    if not doi:
        print(f"No identifier (DOI), record ignored. ({path})")
        return False

    cur = conn.cursor()

    # Model
    cur.execute(
        "INSERT OR IGNORE INTO Model (doi, title, formalism, scale, year, software, source_db) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (doi, m.get("title"), m.get("formalism"), m.get("scale"),
         m.get("year"), m.get("software"), m.get("source_db"))
    )

    # Diseases
    for d in parsed["diseases"]:
        cur.execute("INSERT OR IGNORE INTO Disease VALUES (?, ?)",
                    (d["doid_mondo"], d["name"]))
        cur.execute("INSERT OR IGNORE INTO model_disease VALUES (?, ?)",
                    (doi, d["doid_mondo"]))

    # Pathogens
    for p in parsed["pathogens"]:
        cur.execute("INSERT OR IGNORE INTO Pathogen VALUES (?, ?, ?, ?)",
                    (p["taxonomic_id"], p["species"],
                     p.get("class", ""), p.get("life_cycle_stage", "")))
        cur.execute("INSERT OR IGNORE INTO model_pathogen VALUES (?, ?)",
                    (doi, p["taxonomic_id"]))

    # Hosts
    for h in parsed["hosts"]:
        cur.execute("INSERT OR IGNORE INTO Host VALUES (?, ?, ?)",
                    (h["taxonomic_id"], h["species"],
                     h.get("in_vitro_vivo", "")))
        cur.execute("INSERT OR IGNORE INTO model_host VALUES (?, ?)",
                    (doi, h["taxonomic_id"]))

    # Tissues
    for t in parsed["tissues"]:
        cur.execute("INSERT OR IGNORE INTO Tissue VALUES (?, ?)",
                    (t["ontology_id"], t["name"]))
        cur.execute("INSERT OR IGNORE INTO model_tissue VALUES (?, ?)",
                    (doi, t["ontology_id"]))

    # Cell types
    for tc in parsed["cell_types"]:
        cur.execute("INSERT OR IGNORE INTO Cell_types VALUES (?, ?)",
                    (tc["cell_ontology"], tc.get("name", "")))
        cur.execute("INSERT OR IGNORE INTO model_cell_types VALUES (?, ?)",
                    (doi, tc["cell_ontology"]))

    # Datasets
    for ds in parsed["datasets"]:
        cur.execute("INSERT OR IGNORE INTO Dataset VALUES (?, ?, ?)",
                    (ds["identifier"], ds.get("source", ""), ds.get("description", "")))
        cur.execute("INSERT OR IGNORE INTO model_dataset VALUES (?, ?)",
                    (doi, ds["identifier"]))

    conn.commit()
    return True

# =====================================================================
# Scan the GitHub repo
# =====================================================================

GITHUB_REPO   = "thalieh/InfectioGIT"
GITHUB_BRANCH = "main"
METADATA_FOLDER = "Results"  # only folder to scan in the repo

# Optional : GitHub token to avoid rate-limit 
GITHUB_TOKEN = ""  # fill in if needed

def _github_headers():
    h = {"Accept": "application/vnd.github+json"}
    if GITHUB_TOKEN:
        h["Authorization"] = f"token {GITHUB_TOKEN}"
    return h

def list_json_files_github(
    repo=GITHUB_REPO,
    branch=GITHUB_BRANCH,
    folder=METADATA_FOLDER,
) -> list:
    """
    Lists all .json files in `folder/` via the GitHub API.
    Returns a list of relative paths (ex: 'metadata/biomodels/BIOMD0000000958.json').
    """
    url = f"https://api.github.com/repos/{repo}/git/trees/{branch}?recursive=1"
    resp = requests.get(url, headers=_github_headers(), timeout=30)

    tree = resp.json().get("tree", [])
    json_files = [
        item["path"]
        for item in tree
        if item["type"] == "blob"
        and item["path"].startswith(folder)
        and item["path"].endswith(".json")
    ]
    return json_files


def fetch_json_github(path: str, repo=GITHUB_REPO, branch=GITHUB_BRANCH):
    # Downloads the JSON
    # Using the raw github.com URL (not the API) to let GitHub handle LFS
    url = f"https://github.com/{repo}/raw/{branch}/{path}"
    resp = requests.get(url, timeout=30)
    return resp.json()

def scan_github_and_populate(
    conn: sqlite3.Connection,
    repo=GITHUB_REPO,
    branch=GITHUB_BRANCH,
    folder=METADATA_FOLDER,
):
    json_files = list_json_files_github(repo, branch, folder)
    print(f"{len(json_files)} JSON files found\n")

    success, errors = 0, 0

    for path in json_files:
        #print(f"{path}")
        data = fetch_json_github(path, repo, branch)
        if data is None:
            errors += 1
            continue

        source = detect_source(data)

        if source == "unknown":
            print(f"Unknown source, ignored. ({path})")
            errors += 1
            continue

        parsed = parse_metadata(data, path)
        ok = insert_parsed(conn, parsed, path)
        if ok:
            success += 1
        else:
            errors += 1

    print(f"\nDone : {success} insertions, {errors} errors/ignored")

# =====================================================================
# Launch the scan of the GitHub InfectioGIT repo
# =====================================================================

scan_github_and_populate(conn)

# =====================================================================
# Check and visualize the database
# =====================================================================

# Summary of the database content
tables = [
    "Model", "Disease", "Pathogen", "Host",
    "Tissue", "Cell_types", "Dataset",
    "model_disease", "model_pathogen", "model_host",
    "model_tissue", "model_cell_types", "model_dataset",
]
for t in tables:
    n = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
    print(f"  {t:<35} : {n:>5} row(s)")

pd.read_sql_query("SELECT * FROM Model", conn)

# =====================================================================
# Joins
# =====================================================================

# models × diseases
pd.read_sql_query("""
    SELECT m.title, m.formalism, m.year, m.source_db,
           ma.name AS disease, ma.doid_mondo
    FROM Model m
    JOIN model_disease mm ON m.doi = mm.model_doi
    JOIN Disease ma        ON mm.disease_doid = ma.doid_mondo
""", conn)

# models × pathogens × hosts
pd.read_sql_query("""
    SELECT m.title, m.year,
           p.species AS pathogen, p.class,
           h.species AS host
    FROM Model m
    LEFT JOIN model_pathogen mp ON m.doi = mp.model_doi
    LEFT JOIN Pathogen p         ON mp.pathogen_taxonomy = p.taxonomic_id
    LEFT JOIN model_host mh      ON m.doi = mh.model_doi
    LEFT JOIN Host h             ON mh.host_taxonomy = h.taxonomic_id
""", conn)

# =====================================================================
# Close the connection
# =====================================================================

conn.close()
print("Connection closed")
