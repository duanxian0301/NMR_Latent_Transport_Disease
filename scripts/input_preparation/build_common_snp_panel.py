#!/usr/bin/env python3
"""Build common SNP panel after BOTH UKBB and Karjalainen have 52 hm3.

Outputs:
  results/week1/common_snps.txt
  results/week1/common_snp_panel_stats.md
  results/week1/allele_audit_10k.tsv
  inputs/panels/ukbb_eur_z/panel_{250,500,1000}.z.tsv.gz
  inputs/panels/karjalainen_z/panel_{250,500,1000}.z.tsv.gz
  inputs/panels/ukbb_eur_z/panel_{250,500,1000}.a1.tsv
  inputs/panels/karjalainen_z/panel_{250,500,1000}.a1.tsv
  (+ matching se matrices)

Does NOT run F-GWAS/MR. Does NOT delete sumstats sources.
"""
from __future__ import annotations

import csv
import gzip
import math
import os
import random
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

PROJ = Path(os.environ.get("NMR_GFA_PROJECT_ROOT", Path.cwd())).resolve()
UK = PROJ / "inputs" / "sumstats" / "ukbb_eur"
KJ = PROJ / "inputs" / "sumstats" / "karjalainen"
CROSS = PROJ / "inputs" / "crosswalk" / "matched52_trait_gcst_list.tsv"
OUT_RES = PROJ / "results" / "week1"
OUT_UKZ = PROJ / "inputs" / "panels" / "ukbb_eur_z"
OUT_KJZ = PROJ / "inputs" / "panels" / "karjalainen_z"
PANEL_DEFS = {
    250: PROJ / "inputs" / "crosswalk" / "panel_250_definition.tsv",
    500: PROJ / "inputs" / "crosswalk" / "panel_500_definition.tsv",
    1000: PROJ / "inputs" / "crosswalk" / "panel_1000_definition.tsv",
}
TZ = timezone(timedelta(hours=8))
MIN_COMMON = 100_000
TARGET_COMMON = 200_000


def now() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S %z") + " (Asia/Shanghai)"


def load_crosswalk():
    rows = list(csv.DictReader(open(CROSS, encoding="utf-8"), delimiter="\t"))
    if len(rows) != 52:
        raise SystemExit(f"crosswalk rows={len(rows)} expected 52")
    return rows


def hm3_path(side: str, accession: str) -> Path:
    root = UK if side == "ukbb" else KJ
    p = root / f"{accession}.hm3.tsv.gz"
    if not p.exists():
        raise FileNotFoundError(p)
    return p


def detect_schema(header: list[str]) -> dict:
    cols = {c: i for i, c in enumerate(header)}
    schema = {"kind": None}
    if "rsid" in cols and "beta" in cols and "standard_error" in cols:
        schema.update(
            kind="catalog",
            rsid=cols["rsid"],
            a1=cols.get("effect_allele"),
            a2=cols.get("other_allele"),
            beta=cols["beta"],
            se=cols["standard_error"],
            p=cols.get("p_value", cols.get("neg_log_10_p_value")),
            p_is_neglog="neg_log_10_p_value" in cols,
        )
    elif "SNP" in cols and "BETA" in cols and "SE" in cols:
        schema.update(
            kind="author",
            rsid=cols["SNP"],
            a1=cols.get("A1"),
            a2=cols.get("A2"),
            beta=cols["BETA"],
            se=cols["SE"],
            p=cols.get("P"),
            p_is_neglog=False,
        )
    else:
        raise RuntimeError(f"unknown header: {header[:15]}")
    return schema


def iter_hm3(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        schema = detect_schema(header)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            yield schema, parts


def snp_set(path: Path) -> set[str]:
    s = set()
    for schema, parts in iter_hm3(path):
        i = schema["rsid"]
        if i < len(parts):
            r = parts[i]
            if r and r != "NA":
                s.add(r)
    return s


def snp_alleles(path: Path, wanted: set[str] | None = None) -> dict[str, tuple[str, str]]:
    out = {}
    for schema, parts in iter_hm3(path):
        i = schema["rsid"]
        if i >= len(parts):
            continue
        r = parts[i]
        if not r or r == "NA":
            continue
        if wanted is not None and r not in wanted:
            continue
        a1 = parts[schema["a1"]] if schema["a1"] is not None and schema["a1"] < len(parts) else ""
        a2 = parts[schema["a2"]] if schema["a2"] is not None and schema["a2"] < len(parts) else ""
        out[r] = (a1.upper(), a2.upper())
        if wanted is not None and len(out) >= len(wanted):
            # may still miss some; keep scanning
            pass
    return out


def snp_effects(path: Path, wanted: set[str]) -> dict[str, tuple[float, float, str, str]]:
    """rsid -> (beta, se, a1, a2)"""
    out = {}
    for schema, parts in iter_hm3(path):
        i = schema["rsid"]
        if i >= len(parts):
            continue
        r = parts[i]
        if r not in wanted:
            continue
        try:
            beta = float(parts[schema["beta"]])
            se = float(parts[schema["se"]])
        except (ValueError, IndexError):
            continue
        if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
            continue
        a1 = parts[schema["a1"]].upper() if schema["a1"] is not None else ""
        a2 = parts[schema["a2"]].upper() if schema["a2"] is not None else ""
        out[r] = (beta, se, a1, a2)
    return out


def load_panel_rsids(path: Path) -> list[str]:
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    return [r["rsid"] for r in rows]


def complement(a: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(a, "")


def allele_relation(a1, a2, b1, b2) -> str:
    if not all(len(x) == 1 and x in "ACGT" for x in (a1, a2, b1, b2)):
        return "non_snv_or_missing"
    if a1 == b1 and a2 == b2:
        return "match"
    if a1 == b2 and a2 == b1:
        return "flip"
    ca1, ca2 = complement(a1), complement(a2)
    if ca1 == b1 and ca2 == b2:
        return "complement"
    if ca1 == b2 and ca2 == b1:
        return "complement_flip"
    return "mismatch"


def write_z_panel(side: str, size: int, ordered_rsids: list[str], cross_rows, ref_alleles: dict):
    out_dir = OUT_UKZ if side == "ukbb" else OUT_KJZ
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(ordered_rsids)
    traits = [r["trait_code"] for r in cross_rows]
    # matrix: snp x trait
    zmat = {rs: {} for rs in ordered_rsids}
    semat = {rs: {} for rs in ordered_rsids}
    a1mat = {rs: {} for rs in ordered_rsids}
    n_present = {t: 0 for t in traits}
    n_flipped = 0
    n_missing = 0
    for row in cross_rows:
        trait = row["trait_code"]
        acc = row["ukbb_eur_accession"] if side == "ukbb" else row["karjalainen_accession"]
        eff = snp_effects(hm3_path(side, acc), wanted)
        for rs in ordered_rsids:
            if rs not in eff:
                n_missing += 1
                continue
            beta, se, a1, a2 = eff[rs]
            rel = "match"
            if rs in ref_alleles:
                ra1, ra2 = ref_alleles[rs]
                rel = allele_relation(ra1, ra2, a1, a2)
                if rel in ("flip", "complement_flip"):
                    beta = -beta
                    a1, a2 = a2, a1
                    n_flipped += 1
                elif rel == "complement":
                    a1, a2 = complement(a1), complement(a2)
                elif rel == "mismatch":
                    n_missing += 1
                    continue
            z = beta / se
            zmat[rs][trait] = z
            semat[rs][trait] = se
            a1mat[rs][trait] = a1
            n_present[trait] += 1

    def dump_matrix(path: Path, mat: dict, fmt: str):
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write("rsid\t" + "\t".join(traits) + "\n")
            for rs in ordered_rsids:
                vals = []
                for t in traits:
                    if t in mat[rs]:
                        vals.append(fmt.format(mat[rs][t]) if fmt != "s" else mat[rs][t])
                    else:
                        vals.append("NA")
                f.write(rs + "\t" + "\t".join(vals) + "\n")

    dump_matrix(out_dir / f"panel_{size}.z.tsv.gz", zmat, "{:.6g}")
    dump_matrix(out_dir / f"panel_{size}.se.tsv.gz", semat, "{:.6g}")
    # A1 as plain tsv (small panels)
    with open(out_dir / f"panel_{size}.a1.tsv", "w", encoding="utf-8") as f:
        f.write("rsid\t" + "\t".join(traits) + "\n")
        for rs in ordered_rsids:
            vals = [a1mat[rs].get(t, "NA") for t in traits]
            f.write(rs + "\t" + "\t".join(vals) + "\n")
    return {
        "side": side,
        "panel": size,
        "n_snps": len(ordered_rsids),
        "flips": n_flipped,
        "missing_cells": n_missing,
        "per_trait_present": n_present,
    }


def main() -> int:
    print(f"{now()} BUILD_COMMON_PANEL_START")
    cross = load_crosswalk()
    uk_n = len(list(UK.glob("*.hm3.tsv.gz")))
    kj_n = len(list(KJ.glob("*.hm3.tsv.gz")))
    print(f"ukbb_hm3={uk_n}/52 karj_hm3={kj_n}/52")
    if uk_n < 52 or kj_n < 52:
        raise SystemExit("GATES_NOT_GREEN — refuse to build common panel")
    # verify done sidecars
    for row in cross:
        for side, acc in (("ukbb", row["ukbb_eur_accession"]), ("karjalainen", row["karjalainen_accession"])):
            p = hm3_path("ukbb" if side == "ukbb" else "karjalainen", acc)
            done = Path(str(p) + ".done")
            if not done.exists():
                raise SystemExit(f"missing done sidecar: {done}")

    per_trait_counts = []
    common = None
    uk_sets = []
    kj_sets = []
    for row in cross:
        trait = row["trait_code"]
        ukp = hm3_path("ukbb", row["ukbb_eur_accession"])
        kjp = hm3_path("karjalainen", row["karjalainen_accession"])
        print(f"{now()} SET trait={trait}")
        us = snp_set(ukp)
        ks = snp_set(kjp)
        uk_sets.append(us)
        kj_sets.append(ks)
        inter = us & ks
        per_trait_counts.append(
            {
                "trait_code": trait,
                "ukbb_accession": row["ukbb_eur_accession"],
                "karj_accession": row["karjalainen_accession"],
                "ukbb_hm3_n": len(us),
                "karj_hm3_n": len(ks),
                "pair_intersect_n": len(inter),
            }
        )
        both = us & ks
        common = both if common is None else (common & both)

    assert common is not None
    common_list = sorted(common)
    OUT_RES.mkdir(parents=True, exist_ok=True)
    common_path = OUT_RES / "common_snps.txt"
    with open(common_path, "w", encoding="utf-8") as f:
        for r in common_list:
            f.write(r + "\n")

    # allele audit on random 10k from common (or all if smaller)
    rng = random.Random(42)
    sample_n = min(10_000, len(common_list))
    sample = set(rng.sample(common_list, sample_n)) if sample_n else set()
    # use first trait as UKBB reference alleles vs Karj same trait
    ref_row = cross[0]
    uk_all = snp_alleles(hm3_path("ukbb", ref_row["ukbb_eur_accession"]), sample)
    kj_all = snp_alleles(hm3_path("karjalainen", ref_row["karjalainen_accession"]), sample)
    audit_counts = {"match": 0, "flip": 0, "complement": 0, "complement_flip": 0, "mismatch": 0, "non_snv_or_missing": 0, "absent": 0}
    audit_path = OUT_RES / "allele_audit_10k.tsv"
    with open(audit_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rsid", "uk_a1", "uk_a2", "kj_a1", "kj_a2", "relation", "trait_ref"])
        for r in sorted(sample):
            if r not in uk_all or r not in kj_all:
                audit_counts["absent"] += 1
                w.writerow([r, "", "", "", "", "absent", ref_row["trait_code"]])
                continue
            ua1, ua2 = uk_all[r]
            ka1, ka2 = kj_all[r]
            rel = allele_relation(ua1, ua2, ka1, ka2)
            audit_counts[rel] = audit_counts.get(rel, 0) + 1
            w.writerow([r, ua1, ua2, ka1, ka2, rel, ref_row["trait_code"]])

    # Build GFA-ready small panels restricted to common SNPs, UKBB alleles as reference
    # Reference alleles from UKBB Ala (or first trait) over union of panel rsids
    panel_stats = []
    for size, pdef in PANEL_DEFS.items():
        raw = load_panel_rsids(pdef)
        kept = [r for r in raw if r in common]
        print(f"{now()} PANEL_{size} def={len(raw)} in_common={len(kept)}")
        ref_alleles = snp_alleles(hm3_path("ukbb", cross[0]["ukbb_eur_accession"]), set(kept))
        uk_stat = write_z_panel("ukbb", size, kept, cross, ref_alleles)
        kj_stat = write_z_panel("karjalainen", size, kept, cross, ref_alleles)
        panel_stats.append({"size": size, "def_n": len(raw), "in_common": len(kept), "uk": uk_stat, "kj": kj_stat})

    # stats md
    uk_ns = [r["ukbb_hm3_n"] for r in per_trait_counts]
    kj_ns = [r["karj_hm3_n"] for r in per_trait_counts]
    uk_ns_sorted = sorted(uk_ns)
    kj_ns_sorted = sorted(kj_ns)
    med = lambda xs: xs[len(xs) // 2]
    alert = ""
    if len(common_list) < MIN_COMMON:
        alert = f"**ALERT:** common SNP n={len(common_list)} < {MIN_COMMON}"
    elif len(common_list) < TARGET_COMMON:
        alert = f"**NOTE:** common SNP n={len(common_list)} < target {TARGET_COMMON} (still usable if panels OK)"
    else:
        alert = f"Common SNP n={len(common_list)} meets target ≥{TARGET_COMMON}."

    counts_tsv = OUT_RES / "common_snp_per_trait.tsv"
    with open(counts_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_trait_counts[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(per_trait_counts)

    md = []
    md.append("# Common SNP panel stats")
    md.append("")
    md.append(f"Updated: {now()}")
    md.append("")
    md.append("## Gates")
    md.append(f"- UKBB_EUR HM3: **{uk_n}/52**")
    md.append(f"- Karjalainen HM3: **{kj_n}/52**")
    md.append(f"- Crosswalk matched-52: **{len(cross)}**")
    md.append("")
    md.append("## Intersection (all 52×2)")
    md.append(f"- Common rsIDs: **{len(common_list)}** → `results/week1/common_snps.txt`")
    md.append(f"- UKBB HM3 per-trait: min={min(uk_ns)} median={med(uk_ns_sorted)} max={max(uk_ns)}")
    md.append(f"- Karj HM3 per-trait: min={min(kj_ns)} median={med(kj_ns_sorted)} max={max(kj_ns)}")
    md.append(f"- {alert}")
    md.append("")
    md.append("## Allele audit (random 10k on first trait Ala pair)")
    for k, v in audit_counts.items():
        md.append(f"- {k}: {v}")
    md.append(f"- Detail: `results/week1/allele_audit_10k.tsv`")
    md.append("")
    md.append("## GFA/GLEANR/FactorGo prep panels (restricted to common)")
    for ps in panel_stats:
        md.append(
            f"- panel_{ps['size']}: def={ps['def_n']}, in_common={ps['in_common']}; "
            f"UK flips={ps['uk']['flips']} missing_cells={ps['uk']['missing_cells']}; "
            f"KJ flips={ps['kj']['flips']} missing_cells={ps['kj']['missing_cells']}"
        )
        md.append(f"  - `inputs/panels/ukbb_eur_z/panel_{ps['size']}.z.tsv.gz`")
        md.append(f"  - `inputs/panels/karjalainen_z/panel_{ps['size']}.z.tsv.gz`")
    md.append("")
    md.append("## Explicit non-goals")
    md.append("- No F-GWAS / MR")
    md.append("- No fit-driven trait deletion")
    md.append("- Full-genome Z matrices not written (use documented 250/500/1000 panels for Gate C)")
    md.append("")
    (OUT_RES / "common_snp_panel_stats.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"{now()} BUILD_COMMON_PANEL_DONE common={len(common_list)}")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    sys.exit(main())
