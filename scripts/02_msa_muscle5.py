#!/usr/bin/env python
# coding: utf-8

# In[8]:


from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
from IPython.display import display


# In[9]:


# MMseqs workspace 
MMSEQS_ROOT_WIN = Path(r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MMseqs workspace 11.25")
MMSEQS_ROOT_WSL = "/mnt/c/Users/jaidy/OneDrive - University of Glasgow/Documents/Honours Project/MMseqs workspace 11.25"

HITS_TSV_WIN = MMSEQS_ROOT_WIN / "results" / "search_hits.annot.tsv"
TARGET_DB_WSL = f"{MMSEQS_ROOT_WSL}/data/myzozoaDB/myzozoaDB_2511"

# MUSCLE workspace
MUSCLE_ROOT_WIN = Path(r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MUSCLE5 workspace 01.26")
MUSCLE_ROOT_WSL = "/mnt/c/Users/jaidy/OneDrive - University of Glasgow/Documents/Honours Project/MUSCLE5 workspace 01.26"

OUT_WIN = MUSCLE_ROOT_WIN / "results"
OUT_WIN.mkdir(parents=True, exist_ok=True)
OUT_WSL = f"{MUSCLE_ROOT_WSL}/results"

# MUSCLE executable inside micromamba env (WSL path)
MUSCLE_BIN = "/home/jaidy/.local/share/mamba/envs/muscle5/bin/muscle"

# Pre-MSA filtering thresholds (MMseqs-derived)
MIN_QCOV   = 0.15
MIN_TCOV   = 0.15
MAX_EVALUE = 1e-20

# Which subunit families to align (edit list as needed)
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


# In[14]:


def run_wsl(cmd):
    out = subprocess.run(["wsl", "bash", "-lc", cmd], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"WSL failed:\n{cmd}\n\nSTDERR:\n{out.stderr}")
    return out


COLS = [
    "query","target","fident", "alnlen",
    "qstart","qend","tstart","tend",
    "evalue","bits","qcov","tcov",
    "taxid","taxname","taxlineage"
]

# Filter MMseqs hits to remove low-quality matches, then collapse to one representative per taxid.
def filter_and_represent(df):
    good = df[
        (df["qcov"] >= MIN_QCOV) &
        (df["tcov"] >= MIN_TCOV) &
        (df["evalue"] <= MAX_EVALUE)
    ].copy()

    # Best hit per taxon by bitscore
    good = good.sort_values(["taxid", "bits"], ascending=[True, False])
    rep = good.drop_duplicates(subset=["taxid"], keep="first").copy()
    return rep

# Defining metrics for summary table 
def alignment_summary_windows(fasta_path):
    seqs = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        curr = ""
        for line in f:
            if line.startswith(">"):
                if curr:
                    seqs.append(curr)
                curr = ""
            else:
                curr += line.strip()
        if curr:
            seqs.append(curr)

    if not seqs:
        return {"n_sequences": 0, "gap_fraction": np.nan}

    aln = np.array([list(s) for s in seqs])
    gap_fraction = float(np.mean(aln == "-"))

    return {
        "n_sequences": int(len(seqs)),
        "gap_fraction": round(gap_fraction, 3),
    }

    # Alignment matrix
    aln = np.array([list(s) for s in seqs])
    n_seq, aln_len = aln.shape

    # Overall gap fraction (your old metric, kept)
    gap_fraction = float(np.mean(aln == "-"))

    # Per-sequence ungapped lengths (much more interpretable than alignment_length)
    ungapped_lengths = np.array([len(s.replace("-", "")) for s in seqs], dtype=float)
    median_ungapped_length = float(np.median(ungapped_lengths))
    iqr_ungapped_length = float(np.percentile(ungapped_lengths, 75) - np.percentile(ungapped_lengths, 25))

    # How many columns contain at least one real residue (not all gaps)
    non_gap_col_fraction = float(np.mean(np.any(aln != "-", axis=0)))

    return {
        "n_sequences": int(n_seq),
        "gap_fraction": round(gap_fraction, 3),
        "median_ungapped_length": round(median_ungapped_length, 1),
        "iqr_ungapped_length": round(iqr_ungapped_length, 1),
        "non_gap_col_fraction": round(non_gap_col_fraction, 3),
    }
def run_query_msa(hits, query_name):
    sub = hits[hits["query"] == query_name].copy()
    n_raw = len(sub)

    rep = filter_and_represent(sub)
    n_rep = len(rep)

    # Folder per query 
    out_dir_win = OUT_WIN / query_name
    out_dir_win.mkdir(parents=True, exist_ok=True)
    out_dir_wsl = f"{OUT_WSL}/{query_name}"

    # Output files
    ids_win  = out_dir_win / f"{query_name}_rep_ids.txt"
    aln_win  = out_dir_win / f"{query_name}_rep_seqs.muscle5.fasta"

    # Store path in a portable way 
    aln_rel = str(aln_win.relative_to(OUT_WIN)) if aln_win.exists() else ""

    # If already aligned, skip recomputation (MMseqs-style)
    if aln_win.exists():
        stats = alignment_summary_windows(aln_win)
        return {
            "query": query_name,
            "raw_hits": int(n_raw),
            "representatives": int(n_rep),
            **stats
        }

    # If nothing to align after filtering, return empty stats
    if n_rep == 0:
        return {
            "query": query_name,
            "raw_hits": int(n_raw),
            "representatives": 0,
            "n_sequences": 0,
            "gap_fraction": None,
        }

    # Write representative target IDs 
    rep["target"].to_csv(ids_win, index=False, header=False)

    # WSL file paths
    ids_wsl   = f"{out_dir_wsl}/{query_name}_rep_ids.txt"
    subdb_wsl = f"{out_dir_wsl}/{query_name}_rep_subdb"
    fasta_wsl = f"{out_dir_wsl}/{query_name}_rep_seqs.fasta"
    aln_wsl   = f"{out_dir_wsl}/{query_name}_rep_seqs.muscle5.fasta"

    # Extract + align
    run_wsl(f'''
    mkdir -p "{out_dir_wsl}"
    mmseqs createsubdb "{ids_wsl}" "{TARGET_DB_WSL}" "{subdb_wsl}"
    mmseqs convert2fasta "{subdb_wsl}" "{fasta_wsl}"
    "{MUSCLE_BIN}" -align "{fasta_wsl}" -output "{aln_wsl}"
    ''')

    # Summarise alignment 
    stats = alignment_summary_windows(aln_win)
    return {
        "query": query_name,
        "raw_hits": int(n_raw),
        "representatives":(n_rep),
        **stats
    }


# In[16]:


assert HITS_TSV_WIN.exists(), f"Missing MMseqs TSV: {HITS_TSV_WIN}"

hits = pd.read_csv(HITS_TSV_WIN, sep="\t", header=None, names=COLS)

rows = []
for q in TOXO_QUERIES:
    print("Running:", q)
    rows.append(run_query_msa(hits, q))

summary_df = pd.DataFrame(rows)
cols = ["query", "raw_hits", "representatives", "n_sequences", "gap_fraction"]
summary_df = summary_df[cols].sort_values(["query"])

summary_csv = OUT_WIN / "msa_summary_muscle5.csv"
summary_df.to_csv(summary_csv, index=False)

print("Saved:", summary_csv)
summary_df


# In[1]:


import os
os.getcwd()


# In[ ]:




