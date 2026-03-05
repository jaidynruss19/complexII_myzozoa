import os
from pathlib import Path

# Detect if running inside WSL
def running_in_wsl():
    if "WSL_DISTRO_NAME" in os.environ:
        return True
    try:
        return "microsoft" in os.uname().release.lower()
    except AttributeError:
        return False

IN_WSL = running_in_wsl()

# Project paths 
PROJECT_ROOT_WIN = Path(r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MMseqs workspace 11.25")
WSL_ROOT = "/mnt/c/Users/jaidy/OneDrive - University of Glasgow/Documents/Honours Project/MMseqs workspace 11.25"

DATA_DIR = PROJECT_ROOT_WIN / "data"
RESULTS_DIR = PROJECT_ROOT_WIN / "results"

# Active roots depending on environment
PROJECT_ROOT = WSL_ROOT if IN_WSL else PROJECT_ROOT_WIN
RESULTS_ROOT = f"{WSL_ROOT}/results" if IN_WSL else RESULTS_DIR
DATA_ROOT    = f"{WSL_ROOT}/data" if IN_WSL else DATA_DIR

# Tree file
TREE_FILE = DATA_DIR / "myzozoaDB" / "myzozoaDB_2511.clade.tree"

# MMseqs paths (WSL execution)
QUERY_DB  = f"{WSL_ROOT}/data/ComplexIIDB/QueryDB"
TARGET_DB = f"{WSL_ROOT}/data/myzozoaDB/myzozoaDB_2511"
RESULT_DB = f"{WSL_ROOT}/results/search_hits"
TMP_DIR   = f"{WSL_ROOT}/tmp/mmseqs_search"

# Output TSV
TSV_ANNOT = RESULTS_DIR / "search_hits.annot.tsv"
TSV_ANNOT_WSL = f"{WSL_ROOT}/results/search_hits.annot.tsv"

# MMseqs workspace
MMSEQS_ROOT_WIN = PROJECT_ROOT_WIN
MMSEQS_ROOT_WSL = WSL_ROOT

HITS_TSV_WIN = MMSEQS_ROOT_WIN / "results" / "search_hits.annot.tsv"
TARGET_DB_WSL = f"{MMSEQS_ROOT_WSL}/data/myzozoaDB/myzozoaDB_2511"

# MUSCLE workspace
MUSCLE_ROOT_WIN = Path(r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MUSCLE5 workspace 01.26")
MUSCLE_ROOT_WSL = "/mnt/c/Users/jaidy/OneDrive - University of Glasgow/Documents/Honours Project/MUSCLE5 workspace 01.26"

OUT_WIN = MUSCLE_ROOT_WIN / "results"
OUT_WSL = f"{MUSCLE_ROOT_WSL}/results"

# MUSCLE executable (WSL)
MUSCLE_BIN = "/home/jaidy/.local/share/mamba/envs/muscle5/bin/muscle"

# Filtering thresholds
MIN_QCOV   = 0.15
MIN_TCOV   = 0.15
MAX_EVALUE = 1e-20

# Query proteins
TOXO_QUERIES = [
    "SDHA-Tgondii_",
    "SDHB-Tgondii_",
    "SDH10_Tgondii_GT1",
    "SDH11_Tgondii_GT1",
    "SDH15_Tgondii_GT1",
    "SDH18_Tgondii_GT1",
    "SDH23_Tgondii_GT1",
    "SDH31_Tgondii_GT1",
    "MPODD_Tgondii_GT1",
    "SDHA-Tthermophila_",
    "SDHB-Tthermophila_",
]

# Utility
def ensure_dirs():
    """Create result directories safely."""
    if IN_WSL:
        Path(OUT_WSL).mkdir(parents=True, exist_ok=True)
    else:
        OUT_WIN.mkdir(parents=True, exist_ok=True)



