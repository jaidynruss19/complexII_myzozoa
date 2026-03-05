import pandas as pd
import numpy as np

COLS = [
    "query","target","fident","alnlen",
    "qstart","qend","tstart","tend",
    "evalue","bits","qcov","tcov",
    "taxid","taxname","taxlineage"
]

def load_hits(tsv_path):
    """Load annotated MMseqs TSV and name columns"""
    hits = pd.read_csv(tsv_path, sep="\t", header=None)
    if hits.shape[1] != len(COLS):
        raise ValueError(
            f"Expected {len(COLS)} columns but got {hits.shape[1]}. "
            f"Check your --format-output string."
        )
    hits.columns = COLS
    return hits

def process_hits(hits):
    """Convert types + compute derived metrics (fident_pct, subunit)"""
    numeric_cols = ["fident","alnlen","qstart","qend","tstart","tend","evalue","bits","qcov","tcov"]
    hits[numeric_cols] = hits[numeric_cols].apply(pd.to_numeric, errors="coerce")

    hits = hits.dropna(subset=["fident","alnlen","evalue","qstart","qend","qcov","tcov"]).copy()

    qid = hits["query"].astype(str)

    base = qid.str.split(r"[\s|]+").str[0]
    parts = base.str.split("_")

    hits["subunit"] = np.where(
        base.str.contains("-"),
        base,
        np.where(parts.str.len() >= 2, parts.str[0] + "_" + parts.str[1], parts.str[0])
    ).astype(str)

    hits["fident_pct"] = hits["fident"] * 100 if hits["fident"].max() <= 1.5 else hits["fident"]

    hits["taxid"] = hits["taxid"].astype(str)
    hits["taxname"] = hits["taxname"].astype(str)

    return hits
