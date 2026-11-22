#!/usr/bin/env python3
"""Extract per-trait sample N from UKBB/Karj HM3 for GFA/GLEANR/FactorGo."""
from __future__ import annotations
import csv, gzip, os, statistics
from pathlib import Path

PROJ = Path(os.environ.get("NMR_GFA_PROJECT_ROOT", Path.cwd())).resolve()
CROSS = PROJ / "inputs" / "crosswalk" / "matched52_trait_gcst_list.tsv"
UK = PROJ / "inputs" / "sumstats" / "ukbb_eur"
KJ = PROJ / "inputs" / "sumstats" / "karjalainen"
KJ233 = KJ / "233GWAS"
OUT = PROJ / "inputs" / "panels" / "trait_N.tsv"

def sample_N_from_stream(f, header, n_rows: int = 5000) -> float:
    cols = {c: i for i, c in enumerate(header)}
    if "n" in cols:
        ni = cols["n"]
    elif "N" in cols:
        ni = cols["N"]
    else:
        raise RuntimeError(f"no N column: {header[:20]}")
    vals = []
    for _ in range(n_rows):
        line = f.readline()
        if not line:
            break
        parts = line.rstrip("\n").split("\t")
        try:
            vals.append(float(parts[ni]))
        except (ValueError, IndexError):
            continue
    if not vals:
        raise RuntimeError("no N values")
    return float(statistics.median(vals))

def sample_N_hm3(path: Path) -> float:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        return sample_N_from_stream(f, header)

def sample_N_txt(path: Path) -> float:
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        return sample_N_from_stream(f, header)

def main():
    rows = list(csv.DictReader(open(CROSS, encoding="utf-8"), delimiter="\t"))
    out_rows = []
    for r in rows:
        trait = r["trait_code"]
        uk_acc = r["ukbb_eur_accession"]
        kj_acc = r["karjalainen_accession"]
        uk_n = sample_N_hm3(UK / f"{uk_acc}.hm3.tsv.gz")
        kj_path = KJ / f"{kj_acc}.hm3.tsv.gz"
        try:
            kj_n = sample_N_hm3(kj_path)
            kj_src = "hm3"
        except RuntimeError:
            src = KJ233 / f"{kj_acc}.txt"
            if not src.exists():
                raise
            kj_n = sample_N_txt(src)
            kj_src = "233GWAS"
        out_rows.append({
            "trait_code": trait,
            "ukbb_eur_accession": uk_acc,
            "karjalainen_accession": kj_acc,
            "ukbb_N": f"{uk_n:.0f}",
            "karjalainen_N": f"{kj_n:.0f}",
            "karjalainen_N_source": kj_src,
        })
        print(f"{trait}\tuk={uk_n:.0f}\tkj={kj_n:.0f}\tsrc={kj_src}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"WROTE {OUT}", flush=True)

if __name__ == "__main__":
    main()
