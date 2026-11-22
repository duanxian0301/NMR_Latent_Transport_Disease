#!/usr/bin/env python3
"""S3 preflight audit: official inputs, N, munge, HM3, LD scores, factor-rg table.

Read-only with respect to official P108/S2/exp52_9 objects. Writes only under
SINGLE_TRAIT_CROSS_RESOURCE_BASELINE_20261119/.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import sys
from collections import Counter, OrderedDict
from datetime import datetime

ROOT = os.path.abspath(os.environ.get("NMR_GFA_PROJECT_ROOT", os.getcwd()))
OUT = os.path.join(ROOT, "SINGLE_TRAIT_CROSS_RESOURCE_BASELINE_20261119")
PRE = os.path.join(OUT, "preflight")
LDSC = os.path.abspath(os.environ.get("LDSC_ROOT", os.path.join(ROOT, "external", "ldsc")))
REF = os.path.join(LDSC, "eur_w_ld_chr")
HM3LIST = os.path.join(REF, "w_hm3.snplist")
UKB_MUNGED = os.path.join(ROOT, "results/exp52_9/ldsc_univariate/munged")
TRAIT_N_OFF = os.path.join(ROOT, "results/exp52_9/trait_N.tsv")
TRAIT_N_PANEL = os.path.join(ROOT, "inputs/panels/trait_N.tsv")
S0 = os.path.join(ROOT, "INDEPENDENT_GFA_REPLICATION_20261119/S0_N_AND_TRAIT_IDENTITY.tsv")
NREC = os.path.join(ROOT, "audit_followup_20261117/P2_FINAL/KARJALAINEN_N_RECOVERY.tsv")
ALL30 = os.path.join(ROOT, "audit_followup_20261117/P2_FINAL/ALL30_EXTERNAL_VALIDATION.tsv")
ALL30_P108PKG = os.path.join(
    ROOT, "AGING_EXTENSION_20261119/GPT_READ_PACKAGE_20261119/p108/ALL30_EXTERNAL_VALIDATION.tsv"
)

FAMILIES = {
    "Amino_acids": ["Ala", "Gln", "His", "Ile", "Leu", "Phe", "Tyr", "Val"],
    "Fatty_acids": ["DHA", "LA", "MUFA", "Omega_3", "Omega_6", "PUFA", "SFA", "Total_FA"],
    "Energy_ketone": ["Citrate", "Glucose", "Lactate", "Acetate", "Acetone", "bOHbutyrate", "Pyruvate"],
    "Inflammation_other": ["GlycA", "Albumin", "Creatinine"],
}
TRAIT_TO_FAMILY = {}
for fam, traits in FAMILIES.items():
    for t in traits:
        TRAIT_TO_FAMILY[t] = fam

RETAIN13 = [
    ("F9", "LDL_C"),
    ("F11", "GlycA"),
    ("F13", "Lactate"),
    ("F14", "Gln"),
    ("F15", "LDL_size"),
    ("F17", "Val"),
    ("F18", "Tyr"),
    ("F20", "Gln"),
    ("F22", "Citrate"),
    ("F23", "His"),
    ("F24", "Albumin"),
    ("F25", "Glucose"),
    ("F26", "Phe"),
]

KARJ_MUNGE_SEARCH = [
    os.path.join(ROOT, "results/exp52_9/ldsc_univariate/munged"),
    os.path.join(ROOT, "results/exp52_9"),
    os.path.join(ROOT, "AGING_EXTENSION_20261119/phase1b_dominant_baseline/munged/karjalainen"),
    os.path.join(ROOT, "AGING_EXTENSION_20261119/phase1_13x3/munged"),
    os.path.join(ROOT, "audit_followup_20261117/P108_METHOD_GATE/munged"),
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S %z")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def sha256_file(path, max_bytes=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        if max_bytes is None:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        else:
            h.update(f.read(max_bytes))
    return h.hexdigest()


def inspect_gzip_table(path, n_col_candidates=("N", "n", "sample_size"), max_rows=None, sample_tail=False):
    """Read header, row count (optional), unique values of N if present."""
    out = {
        "exists": os.path.isfile(path),
        "bytes": os.path.getsize(path) if os.path.isfile(path) else 0,
        "header": "",
        "n_rows_scanned": 0,
        "n_rows_total_if_full": "",
        "has_N_col": False,
        "N_unique": "",
        "N_min": "",
        "N_max": "",
        "N_mode": "",
        "N_mode_count": "",
        "gzip_ok": False,
        "error": "",
    }
    if not out["exists"]:
        return out
    try:
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n")
            out["header"] = header
            cols = header.split("\t") if "\t" in header else header.split()
            n_idx = None
            for cand in n_col_candidates:
                if cand in cols:
                    n_idx = cols.index(cand)
                    out["has_N_col"] = True
                    break
            cnt = Counter()
            n = 0
            for line in f:
                n += 1
                if n_idx is not None:
                    parts = line.rstrip("\n").split("\t") if "\t" in line else line.split()
                    if len(parts) > n_idx:
                        try:
                            cnt[int(float(parts[n_idx]))] += 1
                        except Exception:
                            cnt[parts[n_idx]] += 1
                if max_rows is not None and n >= max_rows:
                    break
            out["n_rows_scanned"] = n
            if max_rows is None:
                out["n_rows_total_if_full"] = n
            out["gzip_ok"] = True
            if cnt:
                # numeric if possible
                items = list(cnt.items())
                nums = []
                for k, v in items:
                    try:
                        nums.append((int(k), v, k))
                    except Exception:
                        nums.append((None, v, k))
                out["N_unique"] = len(cnt)
                numeric = [x for x in nums if x[0] is not None]
                if numeric:
                    out["N_min"] = min(x[0] for x in numeric)
                    out["N_max"] = max(x[0] for x in numeric)
                mode = max(items, key=lambda kv: kv[1])
                out["N_mode"] = mode[0]
                out["N_mode_count"] = mode[1]
    except Exception as e:
        out["error"] = repr(e)
    return out


def find_karj_hm3(acc):
    candidates = [
        os.path.join(ROOT, "inputs/sumstats/karjalainen", acc + ".hm3.tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/exp52_9/karjalainen", acc + ".hm3.tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/karjalainen", acc + ".tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/exp52_9/karjalainen", acc + ".tsv.gz"),
    ]
    found = []
    for p in candidates:
        if os.path.isfile(p):
            found.append(p)
    return found, candidates[0]


def find_ukb_hm3(acc):
    candidates = [
        os.path.join(ROOT, "inputs/sumstats/ukbb_eur", acc + ".hm3.tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/exp52_9/ukbb_eur", acc + ".hm3.tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/ukbb_eur", acc + ".tsv.gz"),
        os.path.join(ROOT, "inputs/sumstats/exp52_9/ukbb_eur", acc + ".tsv.gz"),
    ]
    found = [p for p in candidates if os.path.isfile(p)]
    return found


def find_existing_karj_munge(trait, acc):
    hits = []
    names = [
        f"{trait}_{acc}.sumstats.gz",
        f"{acc}.sumstats.gz",
        f"{trait}.sumstats.gz",
    ]
    for folder in KARJ_MUNGE_SEARCH:
        if not os.path.isdir(folder):
            continue
        for nm in names:
            p = os.path.join(folder, nm)
            if os.path.isfile(p):
                hits.append(p)
        # shallow walk one level
        try:
            for ent in os.listdir(folder):
                p2 = os.path.join(folder, ent, f"{trait}_{acc}.sumstats.gz")
                if os.path.isfile(p2):
                    hits.append(p2)
        except Exception:
            pass
    # unique preserve order
    seen = set()
    out = []
    for h in hits:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return out


def assign_family(trait, leftover_lipo):
    if trait in TRAIT_TO_FAMILY:
        return TRAIT_TO_FAMILY[trait]
    leftover_lipo.append(trait)
    return "Lipoproteins_lipids"


def main():
    os.makedirs(PRE, exist_ok=True)
    issues = []
    notes = []

    # --- LDSC / LD scores ---
    ld_ok = True
    ld_rows = []
    ldsc_py = os.path.join(LDSC, "ldsc.py")
    munge_py = os.path.join(LDSC, "munge_sumstats.py")
    version = ""
    if os.path.isfile(ldsc_py):
        with open(ldsc_py, "r", errors="replace") as f:
            for line in f:
                if "__version__" in line and "=" in line:
                    version = line.strip()
                    break
    else:
        ld_ok = False
        issues.append("BLOCKED: official ldsc.py not found at " + ldsc_py)
    if not os.path.isfile(munge_py):
        ld_ok = False
        issues.append("BLOCKED: official munge_sumstats.py not found")
    if not os.path.isfile(HM3LIST):
        ld_ok = False
        issues.append("BLOCKED: missing w_hm3.snplist")
    missing_chr = []
    for chrom in range(1, 23):
        p = os.path.join(REF, f"{chrom}.l2.ldscore.gz")
        ok = os.path.isfile(p)
        ld_rows.append({"chrom": chrom, "path": p, "exists": ok, "bytes": os.path.getsize(p) if ok else 0})
        if not ok:
            missing_chr.append(chrom)
            ld_ok = False
    if missing_chr:
        issues.append("BLOCKED: missing EUR LD score chromosomes: " + ",".join(map(str, missing_chr)))
    notes.append("ldsc_version_string=" + version)
    notes.append("ld_ref=" + REF)
    write_tsv(os.path.join(PRE, "LDSC_ENV.tsv"), ld_rows, ["chrom", "path", "exists", "bytes"])

    # --- official tables ---
    for req in [S0, TRAIT_N_OFF, NREC, ALL30]:
        if not os.path.isfile(req):
            issues.append("BLOCKED: missing required official table " + req)

    s0 = read_tsv(S0)
    if len(s0) != 61:
        issues.append(f"BLOCKED: S0 trait identity has {len(s0)} rows, expected 61")

    # official trait_N
    tn_off = read_tsv(TRAIT_N_OFF)
    tn_panel = read_tsv(TRAIT_N_PANEL) if os.path.isfile(TRAIT_N_PANEL) else []
    rec = read_tsv(NREC)
    rec_by = {r["trait_code"]: r for r in rec}

    # trait_N.tsv schema unknown; inspect
    tn_cols = list(tn_off[0].keys()) if tn_off else []
    notes.append("official_trait_N_cols=" + ",".join(tn_cols))
    # try to index official Karj N
    def lookup_official_n(rows, trait, resource_hints):
        # flexible
        hits = []
        for r in rows:
            keys = {k.lower(): k for k in r.keys()}
            trait_key = None
            for cand in ["trait_code", "trait", "biomarker", "name", "code"]:
                if cand in keys:
                    trait_key = keys[cand]
                    break
            if trait_key is None:
                # first column
                trait_key = list(r.keys())[0]
            if str(r[trait_key]).strip() != trait:
                continue
            # resource filter if present
            res_key = None
            for cand in ["resource", "cohort", "source", "dataset"]:
                if cand in keys:
                    res_key = keys[cand]
                    break
            if res_key and resource_hints:
                val = str(r[res_key]).strip().lower()
                if not any(h in val for h in resource_hints):
                    continue
            hits.append(r)
        return hits

    leftover_lipo = []
    audit_rows = []

    for recs0 in s0:
        trait = recs0["trait_code"].strip()
        ukb_acc = recs0["ukbb_eur_accession_traitN"].strip()
        kj_acc = recs0["karjalainen_accession_traitN"].strip()
        ukb_n_s0 = recs0.get("ukbb_N_traitN", "").strip()
        kj_n_s0 = recs0.get("karjalainen_N_traitN", "").strip()
        kj_n_rec = recs0.get("karjalainen_N_recovered_P3", "").strip()
        family = assign_family(trait, leftover_lipo)

        # official trait_N lookup
        off_kj_hits = lookup_official_n(tn_off, trait, ["karj", "kj", "external"])
        off_ukb_hits = lookup_official_n(tn_off, trait, ["ukb", "ukbb"])
        # If no resource column, rows may contain both N columns
        official_kj_n = ""
        official_ukb_n = ""
        official_kj_raw = ""
        if tn_off:
            # also try exact column names
            keys = {k.lower(): k for k in tn_off[0].keys()}
            # find matching trait row(s)
            all_hits = lookup_official_n(tn_off, trait, None)
            if all_hits:
                r0 = all_hits[0]
                official_kj_raw = str(r0)
                for cand in [
                    "karjalainen_n",
                    "n_karjalainen",
                    "n_karj",
                    "karj_n",
                    "n_kj",
                    "karjalainen",
                ]:
                    if cand in keys and official_kj_n == "":
                        official_kj_n = str(r0.get(keys[cand], "")).strip()
                for cand in ["ukbb_eur_n", "ukb_n", "n_ukb", "n_ukbb", "ukbb_n", "n"]:
                    if cand in keys and official_ukb_n == "":
                        official_ukb_n = str(r0.get(keys[cand], "")).strip()
            # if resource-long format
            if not official_kj_n and off_kj_hits:
                r0 = off_kj_hits[0]
                for cand in ["n", "N", "sample_size"]:
                    if cand in r0 or cand.lower() in {k.lower() for k in r0}:
                        kk = cand if cand in r0 else [k for k in r0 if k.lower() == cand.lower()][0]
                        official_kj_n = str(r0[kk]).strip()
                        break
            if not official_ukb_n and off_ukb_hits:
                r0 = off_ukb_hits[0]
                for cand in ["n", "N", "sample_size"]:
                    if cand in r0 or cand.lower() in {k.lower() for k in r0}:
                        kk = cand if cand in r0 else [k for k in r0 if k.lower() == cand.lower()][0]
                        official_ukb_n = str(r0[kk]).strip()
                        break

        # recovered table
        rec_n = ""
        rec_status = ""
        if trait in rec_by:
            rec_n = str(rec_by[trait].get("recovered_N", "")).strip()
            rec_status = str(rec_by[trait].get("status", "")).strip()

        # Decide analysis N
        # Priority: recovered if official blank or recovered exists as confirmation;
        # official trait_N if present; never invent.
        kj_n_used = ""
        kj_n_source = ""
        kj_n_official_blank = (kj_n_s0 in ("", "NA", "NaN", "nan", ".")) and (official_kj_n in ("", "NA", "NaN", "nan", "."))
        # prefer S0 official column first
        official_recorded = kj_n_s0 if kj_n_s0 not in ("", "NA", "NaN", "nan") else official_kj_n
        if official_recorded not in ("", "NA", "NaN", "nan"):
            try:
                official_recorded_i = str(int(float(official_recorded)))
            except Exception:
                official_recorded_i = official_recorded
        else:
            official_recorded_i = ""

        recovered_i = ""
        if kj_n_rec not in ("", "NA", "NaN", "nan"):
            try:
                recovered_i = str(int(float(kj_n_rec)))
            except Exception:
                recovered_i = kj_n_rec
        elif rec_n not in ("", "NA", "NaN", "nan"):
            try:
                recovered_i = str(int(float(rec_n)))
            except Exception:
                recovered_i = rec_n

        if official_recorded_i and recovered_i:
            kj_n_used = recovered_i  # recovered confirms or is the measurement
            kj_n_source = "official_trait_N_confirmed_by_P3_recovery" if official_recorded_i == recovered_i else "CONFLICT_official_vs_recovered"
            if official_recorded_i != recovered_i:
                issues.append(f"BLOCKED: {trait} official Karj N={official_recorded_i} != recovered {recovered_i}")
        elif recovered_i and not official_recorded_i:
            kj_n_used = recovered_i
            kj_n_source = "P3_recovered_official_trait_N_blank"
        elif official_recorded_i:
            kj_n_used = official_recorded_i
            kj_n_source = "official_trait_N"
        else:
            kj_n_used = ""
            kj_n_source = "MISSING"
            issues.append(f"BLOCKED: cannot determine Karjalainen N for {trait}")

        ukb_n_used = ""
        if ukb_n_s0 not in ("", "NA", "NaN", "nan"):
            try:
                ukb_n_used = str(int(float(ukb_n_s0)))
            except Exception:
                ukb_n_used = ukb_n_s0
        else:
            issues.append(f"BLOCKED: cannot determine UKB N for {trait}")

        # UKB official munge
        ukb_munge = os.path.join(UKB_MUNGED, f"{trait}_{ukb_acc}.sumstats.gz")
        ukb_m = inspect_gzip_table(ukb_munge, max_rows=None)  # full unique N; HM3-size ok
        ukb_n_ok = (
            ukb_m["gzip_ok"]
            and ukb_m["has_N_col"]
            and str(ukb_m["N_unique"]) in ("1", 1)
            and str(ukb_m["N_mode"]) == str(ukb_n_used)
        )
        ukb_reuse = "REUSE_OFFICIAL" if ukb_n_ok else "DO_NOT_REUSE"
        if not ukb_m["exists"]:
            issues.append(f"BLOCKED: missing official UKB munge {ukb_munge}")
            ukb_reuse = "MISSING"
        elif not ukb_n_ok:
            issues.append(
                f"LIMIT: official UKB munge N issue for {trait}: unique={ukb_m['N_unique']} mode={ukb_m['N_mode']} expected={ukb_n_used} err={ukb_m['error']}"
            )

        # UKB HM3 (for completeness only)
        ukb_hm3s = find_ukb_hm3(ukb_acc)

        # Karjalainen HM3
        kj_found, kj_preferred = find_karj_hm3(kj_acc)
        kj_hm3 = kj_found[0] if kj_found else ""
        if not kj_found:
            issues.append(f"BLOCKED: missing complete Karjalainen HM3 for {trait} {kj_acc}")
        kj_h = inspect_gzip_table(kj_hm3, max_rows=5000) if kj_hm3 else {
            "exists": False, "bytes": 0, "header": "", "n_rows_scanned": 0,
            "has_N_col": False, "N_unique": "", "N_min": "", "N_max": "",
            "N_mode": "", "N_mode_count": "", "gzip_ok": False, "error": "missing",
            "n_rows_total_if_full": "",
        }
        # If N column present in first 5k, also note consistency vs used N
        kj_hm3_n_match = ""
        if kj_h.get("has_N_col"):
            if str(kj_h.get("N_unique")) in ("1", 1) and str(kj_h.get("N_mode")) == str(kj_n_used):
                kj_hm3_n_match = "HM3_N_MATCHES_USED"
            else:
                kj_hm3_n_match = f"HM3_N_DIFFERS_unique={kj_h.get('N_unique')}_mode={kj_h.get('N_mode')}"
                notes.append(f"{trait} Karj HM3 has N col {kj_hm3_n_match}")
        else:
            kj_hm3_n_match = "HM3_NO_N_COLUMN"

        # existing Karj munge
        existing_munge = find_existing_karj_munge(trait, kj_acc)
        # if found under official univariate, inspect N
        official_kj_munge = ""
        unofficial_kj_munge = ""
        for p in existing_munge:
            rel = p.replace("\\", "/")
            if "/results/exp52_9/" in rel:
                official_kj_munge = official_kj_munge + (";" if official_kj_munge else "") + p
            else:
                unofficial_kj_munge = unofficial_kj_munge + (";" if unofficial_kj_munge else "") + p

        # header layout for munge flags
        hdr = kj_h.get("header") or ""
        if "effect_allele" in hdr:
            layout = "catalog_harmonized"
        elif "\tA1\t" in ("\t" + hdr.replace(" ", "\t") + "\t") or hdr.startswith("SNP") and "A1" in hdr and "BETA" in hdr:
            layout = "SNP_BETA"
        elif hdr:
            layout = "OTHER:" + hdr[:200]
        else:
            layout = "MISSING"

        row = OrderedDict([
            ("trait", trait),
            ("family", family),
            ("biomarker_name", recs0.get("biomarker_name", "")),
            ("panel_group_s0", recs0.get("group", "")),
            ("ukb_accession", ukb_acc),
            ("karjalainen_accession", kj_acc),
            ("ukb_N_used", ukb_n_used),
            ("karjalainen_N_used", kj_n_used),
            ("karjalainen_N_source", kj_n_source),
            ("karjalainen_N_official_traitN_S0", kj_n_s0),
            ("karjalainen_N_recovered_P3", recovered_i),
            ("official_trait_N_blank", "YES" if kj_n_official_blank or kj_n_s0 in ("", "NA", "NaN", "nan") else "NO"),
            ("ukb_official_munge", ukb_munge if ukb_m["exists"] else ""),
            ("ukb_munge_bytes", ukb_m["bytes"]),
            ("ukb_munge_header", ukb_m["header"]),
            ("ukb_munge_nsnps", ukb_m.get("n_rows_total_if_full", ukb_m["n_rows_scanned"])),
            ("ukb_munge_N_unique", ukb_m["N_unique"]),
            ("ukb_munge_N_mode", ukb_m["N_mode"]),
            ("ukb_munge_gzip_ok", ukb_m["gzip_ok"]),
            ("ukb_munge_reuse", ukb_reuse),
            ("ukb_hm3_found", ";".join(ukb_hm3s)),
            ("karj_hm3_path", kj_hm3),
            ("karj_hm3_n_found", len(kj_found)),
            ("karj_hm3_bytes", kj_h["bytes"]),
            ("karj_hm3_header", kj_h["header"]),
            ("karj_hm3_layout", layout),
            ("karj_hm3_gzip_ok", kj_h["gzip_ok"]),
            ("karj_hm3_rows_scanned_5k", kj_h["n_rows_scanned"]),
            ("karj_hm3_has_N_col", kj_h["has_N_col"]),
            ("karj_hm3_N_note", kj_hm3_n_match),
            ("karj_hm3_error", kj_h["error"]),
            ("official_karj_munge_found", official_kj_munge or "NONE"),
            ("nonofficial_karj_munge_found", unofficial_kj_munge or "NONE"),
        ])
        audit_rows.append(row)

    if len(leftover_lipo) != 35:
        issues.append(f"LIMIT: Lipoproteins_lipids count={len(leftover_lipo)} expected 35: {leftover_lipo}")
    if len(audit_rows) != 61:
        issues.append(f"BLOCKED: audit rows={len(audit_rows)}")

    write_tsv(os.path.join(PRE, "TRAIT_INPUT_AUDIT.tsv"), audit_rows, list(audit_rows[0].keys()))

    # --- ALL30 official factor rg ---
    all30 = read_tsv(ALL30)
    # confirm P108 package copy if present
    all30_same = ""
    if os.path.isfile(ALL30_P108PKG):
        import hashlib as _h

        def _sha(p):
            hh = _h.sha256()
            with open(p, "rb") as f:
                hh.update(f.read())
            return hh.hexdigest()

        all30_same = "IDENTICAL" if _sha(ALL30) == _sha(ALL30_P108PKG) else "DIFFER"
        notes.append("ALL30_P2_vs_P108pkg=" + all30_same)

    rules = set(r.get("rg_usable_rule", "") for r in all30)
    retain = [r for r in all30 if r.get("locked_class") == "Retain" or r.get("factor") in {x[0] for x in RETAIN13}]
    # Use factor id filter to be safe
    retain = [r for r in all30 if r.get("factor") in {x[0] for x in RETAIN13}]
    usable = [r for r in retain if str(r.get("rg_usable", "")).lower() in ("true", "1", "yes")]
    unusable = [r for r in retain if str(r.get("rg_usable", "")).lower() not in ("true", "1", "yes")]
    notes.append("rg_usable_rule=" + " | ".join(sorted(rules)))
    notes.append("retain13_usable=" + ",".join(r["factor"] for r in usable))
    notes.append("retain13_unusable=" + ",".join(r["factor"] for r in unusable))

    factor_rows = []
    for fac, dom in RETAIN13:
        hits = [r for r in all30 if r.get("factor") == fac]
        r = hits[0] if hits else {}
        factor_rows.append({
            "factor": fac,
            "dominant_trait": dom,
            "locked_class": r.get("locked_class", ""),
            "ldsc_rg": r.get("ldsc_rg", ""),
            "ldsc_rg_se": r.get("ldsc_rg_se", ""),
            "rg_usable": r.get("rg_usable", ""),
            "rg_usable_rule": r.get("rg_usable_rule", ""),
            "ldsc_gcov_int": r.get("ldsc_gcov_int", ""),
            "n_snps_ldsc_rg": r.get("n_snps_ldsc_rg", ""),
            "official_source": ALL30,
        })
    write_tsv(os.path.join(PRE, "OFFICIAL_RETAIN13_FACTOR_RG.tsv"), factor_rows,
              ["factor", "dominant_trait", "locked_class", "ldsc_rg", "ldsc_rg_se", "rg_usable",
               "rg_usable_rule", "ldsc_gcov_int", "n_snps_ldsc_rg", "official_source"])

    # Verdict
    blocked = [x for x in issues if x.startswith("BLOCKED")]
    limits = [x for x in issues if not x.startswith("BLOCKED")]
    if blocked:
        verdict = "BLOCKED"
    elif limits:
        verdict = "PASS_WITH_LIMITS"
    else:
        verdict = "PASS"

    n_reuse = sum(1 for r in audit_rows if r["ukb_munge_reuse"] == "REUSE_OFFICIAL")
    n_kj_hm3 = sum(1 for r in audit_rows if r["karj_hm3_gzip_ok"] in (True, "True", "true"))
    n_official_kj_munge = sum(1 for r in audit_rows if r["official_karj_munge_found"] != "NONE")
    n_n_ok = sum(1 for r in audit_rows if r["karjalainen_N_used"] and not str(r["karjalainen_N_source"]).startswith("MISSING") and "CONFLICT" not in str(r["karjalainen_N_source"]))

    md_path = os.path.join(OUT, "PREFLIGHT_NOTE_CN.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# S3 预检说明（61 对同性状跨资源 LDSC 基线）\n\n")
        f.write(f"- 时间（Asia/Shanghai）：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} CST\n")
        f.write(f"- 预检裁决：{verdict}\n")
        f.write(f"- 官方 ldsc：`{ldsc_py}`；版本字符串：`{version}`\n")
        f.write(f"- EUR LD 参考：`{REF}`（1–22 号染色体 l2 + w_hm3.snplist）\n")
        f.write(f"- 官方因子跨资源 rg 来源：`{ALL30}`\n")
        f.write(f"- 可用性规则（原文，冻结）：`{' | '.join(sorted(rules))}`\n")
        f.write(f"- 锁定 13 个 Retain 中，官方可用 rg：{', '.join(r['factor'] for r in usable) or '无'}（n={len(usable)}）\n")
        f.write(f"- 锁定 13 个 Retain 中，官方不可用 rg：{', '.join(r['factor'] for r in unusable) or '无'}（n={len(unusable)}）\n")
        f.write(f"- 官方 UKB munge 可复用（N 与性状 N 一致）：{n_reuse}/61\n")
        f.write(f"- Karjalainen HM3 可读：{n_kj_hm3}/61\n")
        f.write(f"- 官方目录中已有 Karjalainen 性状级 munge：{n_official_kj_munge}/61\n")
        f.write(f"- 已确定 Karjalainen 性状特异 N：{n_n_ok}/61\n")
        f.write(f"- Lipoproteins_lipids 家族计数：{len(leftover_lipo)}（预期 35）\n\n")
        f.write("## 裁决依据\n\n")
        if not issues:
            f.write("未发现阻断项或限制项。\n\n")
        else:
            for x in issues:
                f.write(f"- {x}\n")
            f.write("\n")
        f.write("## 备注\n\n")
        for x in notes:
            f.write(f"- {x}\n")
        f.write("\n")
        f.write("## 下一步（仅在 PASS 或 PASS_WITH_LIMITS 后启动）\n\n")
        f.write("- 复用官方 UKB 单变量 munge，不覆盖官方文件。\n")
        f.write("- 在新目录中，按已核实的性状特异 N 对 Karjalainen HM3 进行 munge；不把因子层常数 N=136016 套到所有代谢物。\n")
        f.write("- 仅运行 61 对同性状双变量 LDSC，不运行 61×61。\n")
        f.write("- 不修改官方 P108、S2、Aging 主结果、手稿或官方 munge / trait_N.tsv。\n")

    # also write issues/notes
    with open(os.path.join(PRE, "PREFLIGHT_ISSUES.txt"), "w", encoding="utf-8") as f:
        f.write("VERDICT\t" + verdict + "\n")
        for x in issues:
            f.write(x + "\n")
        f.write("\nNOTES\n")
        for x in notes:
            f.write(x + "\n")

    print("VERDICT", verdict)
    print("issues", len(issues), "blocked", len(blocked), "limits", len(limits))
    print("wrote", md_path)
    return 0 if verdict != "BLOCKED" else 2


if __name__ == "__main__":
    sys.exit(main())
