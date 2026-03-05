import pandas as pd
import re

COLS = [
    "query","target","fident", "alnlen",
    "qstart","qend","tstart","tend",
    "evalue","bits","qcov","tcov",
    "taxid","taxname","taxlineage"
]


def load_hits(tsv_path):
    """Load annotated MMseqs TSV and name columns"""
    hits = pd.read_csv(tsv_path, sep="\t", header=None)
    if hits.shape[1] != 15:
        raise ValueError(f"Expected 15 columns but got {hits.shape[1]}. Check your --format-output string.")
    hits.columns = [ "qId","tId","fident","alnLen", "qStart","qEnd","tStart","tEnd", "evalue","bitScore","qcov","tcov", "taxid","taxname","taxlineage" ]
    
    return hits

def process_hits(hits):
    """Convert types + compute derived metrics (fident_pct, subunit)"""
    numeric_cols = ["fident","alnLen","qStart","qEnd","tStart","tEnd","evalue","bitScore","qcov","tcov"]
    hits[numeric_cols] = hits[numeric_cols].apply(pd.to_numeric, errors="coerce")

    hits = hits.dropna(subset=["fident","alnLen","evalue","qStart","qEnd","qcov","tcov"]).copy()

    qid = hits["qId"].astype(str)

    # If it already has "-", keep it (e.g. SDHA-Tgondii)
    # If it uses "_" format, keep the first TWO tokens (e.g. SDHTT1_Tthermophila, SDH10_Tgondii)
    parts = qid.str.split("_")
    hits["subunit"] = np.where(
        qid.str.contains("-"),
        qid.str.split(r"\s|\|").str[0],  # take first token before spaces/|
        parts.str[0] + "_" + parts.str[1].fillna("")
    ).astype(str).str.rstrip("_")

    hits["fident_pct"] = hits["fident"] * 100 if hits["fident"].max() <= 1.5 else hits["fident"]

    hits["taxid"] = hits["taxid"].astype(str)
    hits["taxname"] = hits["taxname"].astype(str)

    return hits
