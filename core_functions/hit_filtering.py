import pandas as pd

def filter_good_hits(hits, min_ident_pct=15, min_qcov=0.15, min_tcov=0.15, max_evalue=1e-5):
    """Filter to significant hits (tune thresholds as needed)"""
    return hits[
        (hits["fident_pct"] >= min_ident_pct) &
        (hits["qcov"] >= min_qcov) &
        (hits["tcov"] >= min_tcov) &
        (hits["evalue"] <= max_evalue)
    ].copy()

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
