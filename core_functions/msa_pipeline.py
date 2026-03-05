from pathlib import Path
import pandas as pd

from paths_and_parameters import (
    HITS_TSV_WIN,
    TARGET_DB_WSL,
    OUT_WIN,
    OUT_WSL,
    MUSCLE_BIN,
    MIN_QCOV,
    MIN_TCOV,
    MAX_EVALUE,
    TOXO_QUERIES,
)

from core_functions.mmseqs_search import run_wsl

COLS = [
    "query", "target", "fident", "alnlen",
    "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "qcov", "tcov",
    "taxid", "taxname", "taxlineage"
]


def load_hits_for_msa(tsv_path: Path = HITS_TSV_WIN):
    """Load headerless MMseqs TSV using known column order."""
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        raise FileNotFoundError(f"Missing MMseqs TSV: {tsv_path}")
    return pd.read_csv(tsv_path, sep="\t", header=None, names=COLS)


def filter_and_represent(df: pd.DataFrame):
    """Filter low-quality hits then collapse to one representative per taxid."""
    good = df[
        (df["qcov"] >= MIN_QCOV) &
        (df["tcov"] >= MIN_TCOV) &
        (df["evalue"] <= MAX_EVALUE)
    ].copy()

    if good.empty:
        return good

    # Ensure numeric types for filtering/ranking
    for c in ["qcov", "tcov", "evalue", "bits"]:
        good[c] = pd.to_numeric(good[c], errors="coerce")

    # Ensure taxid is consistent for grouping
    good["taxid"] = good["taxid"].astype(str)

    good = good.dropna(subset=["taxid", "target", "bits"])

    # Prefer best bitscore per taxid
    good = good.sort_values("bits", ascending=False)

    # One rep per taxid
    rep = good.groupby("taxid", as_index=False).head(1).copy()
    return rep


def write_rep_ids(rep: pd.DataFrame, out_ids_win: Path):
    """Write representative target IDs to a file (one per line). Return (win_path, wsl_path)."""
    out_ids_win = Path(out_ids_win)
    out_ids_win.parent.mkdir(parents=True, exist_ok=True)

    rep["target"].astype(str).drop_duplicates().to_csv(out_ids_win, index=False, header=False)

    out_ids_wsl = f"{OUT_WSL}/{out_ids_win.name}"
    return out_ids_win, out_ids_wsl


def gap_fraction_from_fasta(aln_fasta_win: Path):
    """Compute fraction of gaps in an alignment FASTA."""
    aln_fasta_win = Path(aln_fasta_win)
    if not aln_fasta_win.exists():
        return float("nan")

    seqs = []
    current = []
    for line in aln_fasta_win.read_text().splitlines():
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        seqs.append("".join(current))

    if not seqs:
        return float("nan")

    total = sum(len(s) for s in seqs)
    gaps = sum(s.count("-") for s in seqs)
    return gaps / total if total else float("nan")


def run_query_msa(hits: pd.DataFrame, query_name: str):
    """Run representative extraction + MUSCLE alignment for a single query."""
    sub = hits[hits["query"] == query_name].copy()
    n_raw = len(sub)

    rep = filter_and_represent(sub)
    n_rep = len(rep)

    out_dir_win = OUT_WIN / query_name
    out_dir_win.mkdir(parents=True, exist_ok=True)
    out_dir_wsl = f"{OUT_WSL}/{query_name}"

    ids_win = out_dir_win / f"{query_name}_rep_ids.txt"
    _, ids_wsl = write_rep_ids(rep, ids_win)

    subdb_wsl = f"{out_dir_wsl}/{query_name}_rep_subdb"
    fasta_wsl = f"{out_dir_wsl}/{query_name}_rep_seqs.fasta"
    aln_wsl = f"{out_dir_wsl}/{query_name}_rep_seqs.muscle5.fasta"

    aln_win = out_dir_win / f"{query_name}_rep_seqs.muscle5.fasta"

    if n_rep == 0:
        return {
            "query": query_name,
            "raw_hits": n_raw,
            "representatives": 0,
            "n_sequences": 0,
            "gap_fraction": float("nan"),
        }

    run_wsl(f'''
    mkdir -p "{out_dir_wsl}"
    mmseqs createsubdb "{TARGET_DB_WSL}" "{ids_wsl}" "{subdb_wsl}"
    mmseqs convert2fasta "{subdb_wsl}" "{fasta_wsl}"
    "{MUSCLE_BIN}" -align "{fasta_wsl}" -output "{aln_wsl}"
    ''')

    gap_frac = gap_fraction_from_fasta(aln_win)

    return {
        "query": query_name,
        "raw_hits": n_raw,
        "representatives": n_rep,
        "n_sequences": n_rep,
        "gap_fraction": gap_frac,
    }


def run_msa_pipeline(tsv_path: Path = HITS_TSV_WIN, queries=None) -> pd.DataFrame:
    """Run MSA pipeline over all queries and save a summary CSV in OUT_WIN."""
    if queries is None:
        queries = TOXO_QUERIES

    OUT_WIN.mkdir(parents=True, exist_ok=True)

    hits = load_hits_for_msa(tsv_path)

    rows = []
    for q in queries:
        print("Running:", q)
        rows.append(run_query_msa(hits, q))

    summary_df = pd.DataFrame(rows)
    summary_df = summary_df[["query", "raw_hits", "representatives", "n_sequences", "gap_fraction"]]
    summary_df.to_csv(OUT_WIN / "msa_summary_muscle5.csv", index=False)
    return summary_df
