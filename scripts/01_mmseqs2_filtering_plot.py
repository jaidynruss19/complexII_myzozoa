#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Imports & Paths
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import altair as alt
import re
from ete3 import Tree

# Project paths - edit as needed
PROJECT_ROOT_WIN = Path(r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MMseqs workspace 11.25")
WSL_ROOT = "/mnt/c/Users/jaidy/OneDrive - University of Glasgow/Documents/Honours Project/MMseqs workspace 11.25"

DATA_DIR = PROJECT_ROOT_WIN / "data"
RESULTS_DIR = PROJECT_ROOT_WIN / "results"

TREE_FILE = DATA_DIR / "myzozoaDB" / "myzozoaDB_2511.clade.tree"

# MMseqs DB paths (WSL)
QUERY_DB  = f"{WSL_ROOT}/data/ComplexIIDB/QueryDB"
TARGET_DB = f"{WSL_ROOT}/data/myzozoaDB/myzozoaDB_2511"
RESULT_DB = f"{WSL_ROOT}/results/search_hits"
TMP_DIR = f"{WSL_ROOT}/tmp/mmseqs_search"

# Output TSV
TSV_ANNOT = RESULTS_DIR / "search_hits.annot.tsv"
TSV_ANNOT_WSL = f"{WSL_ROOT}/results/search_hits.annot.tsv"


# In[9]:


# Define Functions
def run_wsl(cmd, check = True):
    """Run a command in WSL bash and return the CompletedProcess"""
    out = subprocess.run(["wsl", "bash", "-lc", cmd], capture_output=True, text=True)
    if check and out.returncode != 0:
        raise RuntimeError(
            f"WSL command failed (code {out.returncode}).\n"
            f"CMD:\n{cmd}\n\nSTDOUT:\n{out.stdout}\n\nSTDERR:\n{out.stderr}"
        )
    return out

def mmseqs_run_search(threads=8, search_params="-s 7.5"):
    """Run mmseqs search and create a binary result DB."""
    run_wsl(f'mkdir -p "{TMP_DIR}"')
    cmd = (
        f'mmseqs search "{QUERY_DB}" "{TARGET_DB}" "{RESULT_DB}" "{TMP_DIR}" '
        f'{search_params} --threads {threads}'
    )
    return run_wsl(cmd)


def convertalis_with_taxonomy():
    """Convert MMseqs alignment DB to TSV including taxonomy + coverage"""
    cmd = f'''
    mmseqs convertalis "{QUERY_DB}" "{TARGET_DB}" "{RESULT_DB}" "{TSV_ANNOT_WSL}" \
      --format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,taxid,taxname,taxlineage"
    '''
    return run_wsl(cmd)


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


def filter_good_hits(hits, min_ident_pct=15, min_qcov=0.15, min_tcov=0.15, max_evalue=1e-5):
    """Filter to significant hits (tune thresholds as needed)"""
    return hits[
        (hits["fident_pct"] >= min_ident_pct) &
        (hits["qcov"] >= min_qcov) &
        (hits["tcov"] >= min_tcov) &
        (hits["evalue"] <= max_evalue)
    ].copy()


def get_tree_leaf_order(tree_file):
    """Return leaf order (taxids) from the clade tree"""
    tree = Tree(str(tree_file))
    return [str(x) for x in tree.get_leaf_names()]


def class_from_lineage(lineage):
    if pd.isna(lineage):
        return "Unknown"

    s = lineage.lower()

    if "apicomplexa" in s:
        return "Apicomplexa"
    if "ciliophora" in s:
        return "Ciliophora"
    if "dinophyceae" in s:
        return "Dinophyceae"
    if "bacillariophyceae" in s:
        return "Bacillariophyceae"
    if "oomycetes" in s:
        return "Oomycetes"
    if "chlorophyta" in s:
        return "Chlorophyta"
    if "fungi" in s:
        return "Fungi"
    if "metazoa" in s:
        return "Metazoa"

    return "Other Eukaryotes"


def identity_bin(max_ident):
    if max_ident < 25:
        return "0–25% (low)"
    elif max_ident < 60:
        return "25–60% (mid)"
    else:
        return "60–100% (high)"


def summarise_for_plot(good_hits):
    """Summarise per taxon × subunit for dotplot"""
    plot_data = (
        good_hits
        .groupby(["taxid","taxname","subunit"], as_index=False)
        .agg(
            max_ident=("fident_pct","max"),
            max_qcov=("qcov","max"),
            max_tcov=("tcov","max"),
            taxlineage=("taxlineage","first")
        )
    )
    plot_data["class_group"] = plot_data["taxlineage"].apply(class_from_lineage)
    plot_data["identity_bin"] = plot_data["max_ident"].apply(identity_bin)
    return plot_data


def build_taxname_order(good_hits, leaf_order):
    """Tree-ordered list of species names for x-axis"""
    lookup = (
        good_hits[["taxid","taxname"]]
        .drop_duplicates()
        .set_index("taxid")["taxname"]
        .to_dict()
    )
    return [lookup.get(tid, tid) for tid in leaf_order]

def make_dotplot(plot_data, taxname_order, leaf_order):
    # enforce tree order + consistent categories
    plot_data = plot_data[plot_data["taxid"].isin(leaf_order)].copy()
    plot_data["taxname"] = pd.Categorical(plot_data["taxname"], categories=taxname_order, ordered=True)

    class_order = ["Apicomplexa", "Ciliophora", "Dinophyceae", "Bacillariophyceae", "Oomycetes", "Chlorophyta", "Fungi", "Metazoa", "Other Eukaryotes", "Unknown"]
    plot_data["class_group"] = pd.Categorical(plot_data["class_group"], categories=class_order, ordered=True)


    fig_width = min(1600, 18 * len(taxname_order))


    def species_group_from_label(label):
        if "-Tgondii" in label or "_Tgondii" in label:
            return "1_Tgondii"
        if "-Tthermophila" in label or "_Tthermophila" in label:
            return "2_Tthermophila"
        return "3_Other"

    # Rank subunits within each species block (core first, then lineage-specific)
    order = [
        "SDHA", "SDHB", "SDHC", "SDHD",
        "SDH10", "SDH11", "SDH15", "SDH18", "SDH23", "SDH31", "MPODD",
        "SDHTT1", "SDHTT2", "SDHTT3", "SDHTT4", "SDHTT5",
        "SDHTT6", "SDHTT7", "SDHTT8", "SDHTT9", "SDHTT10", "SDHTT11",
    ]

    def subunit_rank(label):
        for i, s in enumerate(order):
            if label.startswith(s):
                return i
        return len(order) + 1

    ycol = "subunit"   

    plot_data["_species_group"] = plot_data[ycol].astype(str).apply(species_group_from_label)
    plot_data["_subunit_rank"] = plot_data[ycol].astype(str).apply(subunit_rank)

    # stable, explicit order for y-axis labels
    y_order = (
        plot_data[[ycol, "_species_group", "_subunit_rank"]]
        .drop_duplicates()
        .sort_values(by=["_species_group", "_subunit_rank", ycol])
        [ycol]
        .tolist()
    )

    plot_data[ycol] = plot_data[ycol].astype("category").cat.set_categories(y_order, ordered=True)
    plot_data = plot_data.sort_values(by=[ycol])


    chart = (
        alt.Chart(plot_data)
        .mark_circle(size=150)
        .encode(
            x=alt.X("taxname:N", sort=taxname_order, axis=alt.Axis(labelAngle=-90, labelFontSize=14, titleFontSize=14), title="Species (tree order)"),
            y=alt.Y("subunit:N", axis=alt.Axis(labelFontSize=13, titleFontSize=14), sort=y_order, title="Complex II subunit"),
            color=alt.Color("class_group:N", title="Taxonomic Class", scale=alt.Scale(scheme="tableau10"), sort=class_order),
            opacity=alt.Opacity(
                "identity_bin:N",
                scale=alt.Scale(
                    domain=["0–25% (low)", "25–60% (mid)", "60–100% (high)"],
                    range=[0.25, 0.7, 1.0]
                ),
                title="% Id"
            ),
            tooltip=["taxid", "taxname", "class_group", "subunit",
    alt.Tooltip("max_ident:Q", format=".1f", title="Max % identity"),
    alt.Tooltip("max_qcov:Q", format=".2f", title="Max query coverage"),
    alt.Tooltip("max_tcov:Q", format=".2f", title="Max target coverage"),
]
,
        )
        .properties(width=fig_width, height=350) 
    )
    chart = chart.configure_legend(labelFontSize=13, titleFontSize=14, symbolSize=150)
    chart = chart.properties(padding={"bottom": 180})
    return chart

def add_subunit_clean(df, label_col= None):
    df = df.copy()

    if label_col is None:
        if "subunit" in df.columns:
            label_col = "subunit"
        elif "query" in df.columns:
            label_col = "query"
        else:
            raise KeyError("Need a label column: expected 'subunit' or 'query'.")

    def clean(label):
        s = str(label)
        m = re.match(r"^(SDHA|SDHB|SDHC|SDHD|MPODD|SDH\d+|SDHTT\d+)", s)
        return m.group(1) if m else s

    df["subunit_clean"] = df[label_col].apply(clean)
    return df

def make_supplement_master_table(all_hits_table):
    df = all_hits_table.copy()
    core = {"SDHA", "SDHB"}

    # Parse recovered subunits into lists
    def to_list(x):
        if isinstance(x, (list, tuple, set)):
            return [str(s).strip() for s in x if str(s).strip()]
        if isinstance(x, str):
            return [s.strip() for s in x.split(",") if s.strip()]
        return []

    rec_list = df["Recovered subunits"].apply(to_list)

    def classify(xs):
        s = set(xs)
        has_core = core.issubset(s)
        extras = s - core
        if has_core and len(extras) == 0:
            return "Core only"
        if has_core and len(extras) > 0:
            return "Core + extra"
        if len(s) > 0 and not has_core:
            return "Partial / non-core"
        return "No hits"

    # Temporary column for ordering only
    df["_recovery_group"] = rec_list.apply(classify)
    df["n_subunits"] = rec_list.apply(lambda xs: len(set(xs)))
    df["Recovered subunits"] = rec_list.apply(lambda xs: ", ".join(xs))

    order = ["Core only", "Core + extra", "Partial / non-core", "No hits"]
    df["_recovery_group"] = pd.Categorical(df["_recovery_group"], categories=order, ordered=True)

    # Sort by recovery class, then taxonomy, then species
    df = (
        df.sort_values(["_recovery_group", "Taxonomic class", "Species"])
          .reset_index(drop=True)
    )

    # Drop the helper column before returning
    return df.drop(columns=["_recovery_group"])


# In[3]:


# Run MMseqs search
run_search = True # Set False is RESULT_DB already exists 

if run_search:
    if run_wsl(f'test -f "{RESULT_DB}.dbtype"', check=False).returncode == 0:
        print("MMseqs result DB already exists — skipping search")
    else:
        mmseqs_run_search()

# Convert MMseqs results to annotated TSV (taxonomy included)
convertalis_with_taxonomy()

# Load + process hits 
hits = load_hits(TSV_ANNOT)
hits = process_hits(hits)

# Filter 
good_hits = filter_good_hits(hits)
if good_hits.empty:
    raise ValueError("No hits passed filters. Relax thresholds or check the search output.")

# Tree order = taxname order 
leaf_order = get_tree_leaf_order(TREE_FILE)
good_hits = good_hits[good_hits["taxid"].isin(leaf_order)].copy()
taxname_order = build_taxname_order(good_hits, leaf_order)

# Summarise + plot 
plot_data = summarise_for_plot(good_hits)
dotplot = make_dotplot(plot_data, taxname_order, leaf_order)

dotplot


# In[8]:


# Summary Table of Hits 
all_hits_table = (
    plot_data.groupby(["taxid", "taxname"], as_index=False)
    .agg(taxlineage=("taxlineage", "first"), **{"Recovered subunits": ("subunit", lambda s: ", ".join(sorted(set(s.dropna()))))})
    .rename(columns={"taxname": "Species"}))

all_hits_table["Taxonomic class"] = all_hits_table["taxlineage"].apply(class_from_lineage)

supplement_master = make_supplement_master_table(all_hits_table)
supplement_master.to_excel("Supplement_Table_AllSpecies_ComplexII.xlsx", index=False)
display(supplement_master)


# In[7]:


# Species with no Hits Summary
species = pd.read_csv(
    r"C:\Users\jaidy\OneDrive - University of Glasgow\Documents\Honours Project\MMseqs workspace 11.25\data\myzozoaDB\myzozoaDB_2511.manual_annotation.tsv",
    sep="\t",
    dtype={"taxid": "string"}
)

hit_taxids = set(good_hits["taxid"].dropna().unique())
no_hits = species[~species["taxid"].isin(hit_taxids)]

print(f"Species with no hits (rows): {len(no_hits)}")
print("Unique taxids in no_hits:", no_hits["taxid"].nunique(dropna=True))
print("Rows with missing taxid:", no_hits["taxid"].isna().sum())

cols = [c for c in ["organism","taxid","generic_label","strain","genome_assembly_accession"] if c in no_hits.columns]
display(
    no_hits[cols]
      .dropna(subset=["taxid"])
      .drop_duplicates(subset=["taxid"])
      .sort_values("organism")
)


# In[ ]:




