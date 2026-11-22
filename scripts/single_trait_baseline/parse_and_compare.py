#!/usr/bin/env python3
"""Parse 61 same-trait UKB x Karjalainen LDSC logs and compare with official 13 factor rg.

Writes only under SINGLE_TRAIT_CROSS_RESOURCE_BASELINE_20261119/.
Does not modify official objects. Does not retune the usable-rg rule.
"""
from __future__ import annotations

import csv
import math
import os
import re
from collections import OrderedDict
from datetime import datetime

ROOT = os.path.abspath(os.environ.get("NMR_GFA_PROJECT_ROOT", os.getcwd()))
OUT = os.path.join(ROOT, "SINGLE_TRAIT_CROSS_RESOURCE_BASELINE_20261119")
AUDIT = os.path.join(OUT, "preflight/TRAIT_INPUT_AUDIT.tsv")
RGDIR = os.path.join(OUT, "ldsc_rg")
KARJ_MUNGE = os.path.join(OUT, "munged/karjalainen")
ALL30 = os.path.join(ROOT, "audit_followup_20261117/P2_FINAL/ALL30_EXTERNAL_VALIDATION.tsv")
S2_MATCH = os.path.join(ROOT, "INDEPENDENT_GFA_REPLICATION_20261119/CROSS_RESOURCE_LOADING_MATCH_RETAIN13.tsv")
S2_CONC = os.path.join(ROOT, "INDEPENDENT_GFA_REPLICATION_20261119/S2_FACTOR_CONCENTRATION.tsv")
S2_STAB = os.path.join(ROOT, "INDEPENDENT_GFA_REPLICATION_20261119/S2_STABILITY_MATCHES.tsv")

RG_USABLE_RULE = (
    "SE<0.30 and 0<=rg<=1.05 and rg/SE>=2 (run_exp52_9_factor_qc.py); "
    "huge-SE rg not treated as usable"
)

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


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


def fmt(x, nd=6):
    if x is None or x == "":
        return ""
    if isinstance(x, bool):
        return "TRUE" if x else "FALSE"
    if isinstance(x, int):
        return str(x)
    try:
        xf = float(x)
        if not math.isfinite(xf):
            return str(x)
        return f"{xf:.10g}"
    except Exception:
        return str(x)


def parse_paren(s):
    """Parse '0.0376 (0.0046)' -> (0.0376, 0.0046)."""
    m = re.search(r"([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*\(\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*\)", s)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def parse_float(s):
    m = re.search(r"([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)", s)
    return float(m.group(1)) if m else None


def section_after(text, title):
    # grab until next blank line followed by a title-like line or end
    pat = re.compile(re.escape(title) + r"\s*-+\s*(.*?)(?=\n[A-Z][^\n]*\n-+|\nSummary of Genetic Correlation Results|\Z)", re.S)
    m = pat.search(text)
    return m.group(1) if m else ""


def parse_h2_block(block):
    out = {
        "h2": None, "h2_se": None, "lambda_gc": None, "mean_chisq": None,
        "intercept": None, "intercept_se": None, "ratio": None, "ratio_se": None,
    }
    for line in block.splitlines():
        if "Total Observed scale h2:" in line:
            out["h2"], out["h2_se"] = parse_paren(line)
        elif line.strip().startswith("Lambda GC:"):
            out["lambda_gc"] = parse_float(line.split(":", 1)[1])
        elif line.strip().startswith("Mean Chi^2:") or line.strip().startswith("Mean Chi2:"):
            out["mean_chisq"] = parse_float(line.split(":", 1)[1])
        elif line.strip().startswith("Intercept:"):
            out["intercept"], out["intercept_se"] = parse_paren(line)
        elif line.strip().startswith("Ratio:"):
            out["ratio"], out["ratio_se"] = parse_paren(line)
    return out


def parse_rg_log(path):
    res = {
        "estimate_available": False,
        "failure_reason": "",
        "ldsc_version": "",
        "n_snps_ukb_sumstats": None,
        "n_snps_karj_sumstats": None,
        "n_snps_valid_alleles": None,
        "n_snps_after_ref_ld": None,
        "n_snps_after_w_ld": None,
        "ukb_h2": None, "ukb_h2_se": None, "ukb_h2_z": None,
        "ukb_lambda_gc": None, "ukb_mean_chisq": None,
        "ukb_ldsc_intercept": None, "ukb_ldsc_intercept_se": None,
        "ukb_ratio": None, "ukb_ratio_se": None,
        "karj_h2": None, "karj_h2_se": None, "karj_h2_z": None,
        "karj_lambda_gc": None, "karj_mean_chisq": None,
        "karj_ldsc_intercept": None, "karj_ldsc_intercept_se": None,
        "karj_ratio": None, "karj_ratio_se": None,
        "rg": None, "rg_se": None, "rg_z": None, "rg_p": None,
        "rg_ci95_low": None, "rg_ci95_high": None,
        "gcov": None, "gcov_se": None, "mean_z1z2": None,
        "gcov_int": None, "gcov_int_se": None,
        "log_path": path,
    }
    if not os.path.isfile(path):
        res["failure_reason"] = "missing_rg_log"
        return res
    text = open(path, encoding="utf-8", errors="replace").read()
    mver = re.search(r"Version\s+([0-9.]+)", text)
    if mver:
        res["ldsc_version"] = mver.group(1)

    # SNP counts in order: first sumstats is UKB (p1), then after LD merge, then p2 Karjalainen
    reads = re.findall(r"Read summary statistics for\s+(\d+)\s+SNPs", text)
    if len(reads) >= 1:
        res["n_snps_ukb_sumstats"] = int(reads[0])
    if len(reads) >= 2:
        res["n_snps_karj_sumstats"] = int(reads[1])
    m = re.search(r"After merging with reference panel LD,\s+(\d+)\s+SNPs remain", text)
    if m:
        res["n_snps_after_ref_ld"] = int(m.group(1))
    m = re.search(r"After merging with regression SNP LD,\s+(\d+)\s+SNPs remain", text)
    if m:
        res["n_snps_after_w_ld"] = int(m.group(1))
    m = re.search(r"(\d+)\s+SNPs with valid alleles", text)
    if m:
        res["n_snps_valid_alleles"] = int(m.group(1))

    # Split phenotype blocks
    p1 = re.search(
        r"Heritability of phenotype 1\s*-+\s*(.*?)(?=Heritability of phenotype 2)",
        text, re.S,
    )
    p2 = re.search(
        r"Heritability of phenotype 2/2\s*-+\s*(.*?)(?=Genetic Covariance)",
        text, re.S,
    )
    gcov = re.search(
        r"Genetic Covariance\s*-+\s*(.*?)(?=Genetic Correlation)",
        text, re.S,
    )
    gcor = re.search(
        r"Genetic Correlation\s*-+\s*(.*?)(?=Summary of Genetic Correlation Results|\Z)",
        text, re.S,
    )
    if p1:
        b = parse_h2_block(p1.group(1))
        res["ukb_h2"], res["ukb_h2_se"] = b["h2"], b["h2_se"]
        res["ukb_lambda_gc"], res["ukb_mean_chisq"] = b["lambda_gc"], b["mean_chisq"]
        res["ukb_ldsc_intercept"], res["ukb_ldsc_intercept_se"] = b["intercept"], b["intercept_se"]
        res["ukb_ratio"], res["ukb_ratio_se"] = b["ratio"], b["ratio_se"]
        if res["ukb_h2"] is not None and res["ukb_h2_se"] not in (None, 0):
            res["ukb_h2_z"] = res["ukb_h2"] / res["ukb_h2_se"]
    if p2:
        b = parse_h2_block(p2.group(1))
        res["karj_h2"], res["karj_h2_se"] = b["h2"], b["h2_se"]
        res["karj_lambda_gc"], res["karj_mean_chisq"] = b["lambda_gc"], b["mean_chisq"]
        res["karj_ldsc_intercept"], res["karj_ldsc_intercept_se"] = b["intercept"], b["intercept_se"]
        res["karj_ratio"], res["karj_ratio_se"] = b["ratio"], b["ratio_se"]
        if res["karj_h2"] is not None and res["karj_h2_se"] not in (None, 0):
            res["karj_h2_z"] = res["karj_h2"] / res["karj_h2_se"]
    if gcov:
        for line in gcov.group(1).splitlines():
            if "Total Observed scale gencov:" in line:
                res["gcov"], res["gcov_se"] = parse_paren(line)
            elif "Mean z1*z2:" in line:
                res["mean_z1z2"] = parse_float(line.split(":", 1)[1])
            elif line.strip().startswith("Intercept:"):
                res["gcov_int"], res["gcov_int_se"] = parse_paren(line)
    if gcor:
        for line in gcor.group(1).splitlines():
            if "Genetic Correlation:" in line:
                res["rg"], res["rg_se"] = parse_paren(line)
            elif line.strip().startswith("Z-score:"):
                res["rg_z"] = parse_float(line.split(":", 1)[1])
            elif re.match(r"\s*P:", line):
                res["rg_p"] = parse_float(line.split(":", 1)[1])

    if res["rg"] is not None and res["rg_se"] is not None:
        res["rg_ci95_low"] = res["rg"] - 1.96 * res["rg_se"]
        res["rg_ci95_high"] = res["rg"] + 1.96 * res["rg_se"]
        if res["rg_z"] is None and res["rg_se"] != 0:
            res["rg_z"] = res["rg"] / res["rg_se"]
        res["estimate_available"] = True
    else:
        if "Genetic Correlation:" not in text:
            res["failure_reason"] = "ldsc_log_missing_Genetic_Correlation"
        else:
            res["failure_reason"] = "could_not_parse_rg_from_log"
    return res


def parse_munge_log(path):
    out = {
        "munge_log": path if os.path.isfile(path) else "",
        "munge_N_used": "",
        "munge_nsnps_written_nonmissing": "",
        "munge_mean_chisq": "",
        "munge_lambda_gc": "",
        "munge_max_chisq": "",
        "munge_gws_snps": "",
    }
    if not os.path.isfile(path):
        return out
    text = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"Using N = ([0-9.]+)", text)
    if m:
        try:
            out["munge_N_used"] = str(int(float(m.group(1))))
        except Exception:
            out["munge_N_used"] = m.group(1)
    m = re.search(r"Writing summary statistics for\s+(\d+)\s+SNPs\s+\((\d+)\s+with nonmissing beta\)", text)
    if m:
        out["munge_nsnps_written_nonmissing"] = int(m.group(2))
    m = re.search(r"Mean chi\^2 = ([0-9.eE+-]+)", text)
    if m:
        out["munge_mean_chisq"] = float(m.group(1))
    m = re.search(r"Lambda GC = ([0-9.eE+-]+)", text)
    if m:
        out["munge_lambda_gc"] = float(m.group(1))
    m = re.search(r"Max chi\^2 = ([0-9.eE+-]+)", text)
    if m:
        out["munge_max_chisq"] = float(m.group(1))
    m = re.search(r"(\d+)\s+Genome-wide significant SNPs", text)
    if m:
        out["munge_gws_snps"] = int(m.group(1))
    return out


def apply_usable_rule(rg, se):
    """Official P2/ALL30 rule, verbatim: SE<0.30 and 0<=rg<=1.05 and rg/SE>=2."""
    if rg is None or se is None:
        return False, "missing_rg_or_se"
    reasons = []
    try:
        rgf = float(rg)
        sef = float(se)
    except Exception:
        return False, "non_numeric_rg_or_se"
    if not math.isfinite(rgf) or not math.isfinite(sef):
        return False, "non_finite_rg_or_se"
    if not (sef < 0.30):
        reasons.append("SE>=0.30")
    if not (0 <= rgf <= 1.05):
        reasons.append("rg_not_in_[0,1.05]")
    if sef == 0 or (not (rgf / sef >= 2)):
        reasons.append("rg/SE<2")
    if reasons:
        return False, ";".join(reasons)
    return True, ""


def bh_qvalues(pvals, m_family):
    """BH q-values with frozen family size m_family. Missing p -> None."""
    indexed = []
    for i, p in enumerate(pvals):
        if p is None:
            continue
        try:
            pf = float(p)
            if math.isfinite(pf) and 0 <= pf <= 1:
                indexed.append((i, pf))
        except Exception:
            continue
    q = [None] * len(pvals)
    if not indexed:
        return q
    indexed.sort(key=lambda x: x[1])
    n = m_family
    running = 1.0
    tmp = [None] * len(indexed)
    for rank_from_end, (k, (i, p)) in enumerate(reversed(list(enumerate(indexed, start=1))), start=0):
        # enumerate(..., start=1) then reversed: k is BH rank i=1..m_valid? We use family m=61
        pass
    # proper: after sorting ascending, q_i = min_{j>=i} (p_j * m / j), capped at 1
    ranks = list(range(1, len(indexed) + 1))
    raw = [indexed[j][1] * n / ranks[j] for j in range(len(indexed))]
    # cumulative min from the tail
    cum = 1.0
    out_sorted = [None] * len(indexed)
    for j in range(len(indexed) - 1, -1, -1):
        cum = min(cum, raw[j])
        out_sorted[j] = min(1.0, cum)
    for (i, _p), qv in zip(indexed, out_sorted):
        q[i] = qv
    return q


def percentile_rank(value, dist):
    """Percentile rank: 100 * (# of dist <= value) / n. None if empty."""
    if value is None or not dist:
        return None
    n = len(dist)
    le = sum(1 for x in dist if x <= value)
    lt = sum(1 for x in dist if x < value)
    # midpoint rank
    return 100.0 * (lt + le) / (2.0 * n)


def ci_overlap(a_lo, a_hi, b_lo, b_hi):
    if None in (a_lo, a_hi, b_lo, b_hi):
        return ""
    return "YES" if (a_lo <= b_hi and b_lo <= a_hi) else "NO"


def lookup_s2(ukb_factor):
    out = {
        "s2_karj_matched_factor": "",
        "s2_abs_cosine": "",
        "s2_congruence_class": "",
        "s2_same_dominant_trait": "",
        "s2_karj_dominant_trait": "",
        "s2_ukb_concentration_share": "",
        "s2_karj_concentration_share": "",
        "s2_ukb_participation_ratio": "",
        "s2_karj_participation_ratio": "",
        "s2_top5_jaccard": "",
        "s2_ukb_odd_even_note": "",
    }
    if os.path.isfile(S2_MATCH):
        rows = read_tsv(S2_MATCH)
        hits = [r for r in rows if r.get("UKB_factor") == ukb_factor and r.get("match_type") == "unique_hungarian"]
        if hits:
            r = hits[0]
            out["s2_karj_matched_factor"] = r.get("Karj_factor", "")
            out["s2_abs_cosine"] = r.get("abs_cosine", "")
            out["s2_congruence_class"] = r.get("congruence_class", "")
            out["s2_same_dominant_trait"] = r.get("same_dominant_trait", "")
            out["s2_karj_dominant_trait"] = r.get("Karj_dominant_trait", "")
            out["s2_ukb_concentration_share"] = r.get("UKB_concentration_share", "")
            out["s2_karj_concentration_share"] = r.get("Karj_concentration_share", "")
            out["s2_ukb_participation_ratio"] = r.get("UKB_participation_ratio", "")
            out["s2_karj_participation_ratio"] = r.get("Karj_participation_ratio", "")
            out["s2_top5_jaccard"] = r.get("top5_jaccard", "")
    return out


def official_oddeven_flag(factor):
    # Locked previous official UKB odd/even result among Retain 13:
    # both odd and even |phi|>=0.85: F11, F18, F22, F23 only.
    # This is recorded, not recomputed.
    both = {"F11", "F18", "F22", "F23"}
    if factor in both:
        return "both_odd_and_even_abs_phi_ge_0.85 (official locked; not rerun)"
    return "not_both_odd_and_even_abs_phi_ge_0.85 (official locked; not rerun)"


def main():
    audit = read_tsv(AUDIT)
    if len(audit) != 61:
        raise SystemExit(f"expected 61 audit rows, got {len(audit)}")

    rows = []
    for a in audit:
        trait = a["trait"]
        acc_k = a["karjalainen_accession"]
        logp = os.path.join(RGDIR, f"{trait}_UKB_vs_Karjalainen.log")
        parsed = parse_rg_log(logp)
        munged = parse_munge_log(os.path.join(KARJ_MUNGE, f"{trait}_{acc_k}.log"))
        usable, why_unusable = apply_usable_rule(parsed["rg"], parsed["rg_se"])
        if not parsed["estimate_available"]:
            usable = False
            if parsed["failure_reason"] and why_unusable:
                why_unusable = parsed["failure_reason"] + ";" + why_unusable
            elif parsed["failure_reason"]:
                why_unusable = parsed["failure_reason"]

        row = OrderedDict()
        row["trait"] = trait
        row["family"] = a["family"]
        row["biomarker_name"] = a.get("biomarker_name", "")
        row["ukb_accession"] = a["ukb_accession"]
        row["karjalainen_accession"] = acc_k
        row["ukb_N"] = a["ukb_N_used"]
        row["karjalainen_N"] = a["karjalainen_N_used"]
        row["karjalainen_N_source"] = a["karjalainen_N_source"]
        row["ukb_munge_path"] = a["ukb_official_munge"]
        row["karj_munge_path"] = os.path.join(KARJ_MUNGE, f"{trait}_{acc_k}.sumstats.gz")
        row["ldsc_log_path"] = logp
        row["ldsc_version"] = parsed["ldsc_version"]
        for k in [
            "ukb_h2", "ukb_h2_se", "ukb_h2_z", "ukb_lambda_gc", "ukb_mean_chisq",
            "ukb_ldsc_intercept", "ukb_ldsc_intercept_se", "ukb_ratio", "ukb_ratio_se",
            "karj_h2", "karj_h2_se", "karj_h2_z", "karj_lambda_gc", "karj_mean_chisq",
            "karj_ldsc_intercept", "karj_ldsc_intercept_se", "karj_ratio", "karj_ratio_se",
            "rg", "rg_se", "rg_z", "rg_p", "rg_ci95_low", "rg_ci95_high",
            "gcov", "gcov_se", "mean_z1z2", "gcov_int", "gcov_int_se",
            "n_snps_ukb_sumstats", "n_snps_karj_sumstats", "n_snps_after_ref_ld",
            "n_snps_after_w_ld", "n_snps_valid_alleles",
        ]:
            row[k] = parsed[k]
        row["karj_munge_N_used"] = munged["munge_N_used"]
        row["karj_munge_nsnps_nonmissing"] = munged["munge_nsnps_written_nonmissing"]
        row["karj_munge_mean_chisq"] = munged["munge_mean_chisq"]
        row["karj_munge_lambda_gc"] = munged["munge_lambda_gc"]
        row["estimate_available"] = "TRUE" if parsed["estimate_available"] else "FALSE"
        row["rg_usable"] = "TRUE" if usable else "FALSE"
        row["rg_usable_rule"] = RG_USABLE_RULE
        row["unusable_or_failure_reason"] = why_unusable
        row["rg_bh_q_family61"] = ""  # filled later
        rows.append((row, parsed["rg_p"]))

    qvals = bh_qvalues([p for _, p in rows], m_family=61)
    for (row, _p), q in zip(rows, qvals):
        row["rg_bh_q_family61"] = q if q is not None else ""

    out_rows = [r for r, _ in rows]
    cols = list(out_rows[0].keys())
    # format numerics
    formatted = []
    for r in out_rows:
        fr = OrderedDict()
        for k in cols:
            fr[k] = fmt(r[k])
        formatted.append(fr)

    tsv61 = os.path.join(OUT, "61_TRAIT_CROSS_RESOURCE_LDSC.tsv")
    write_tsv(tsv61, formatted, cols)

    # --- 13-factor comparison ---
    by_trait = {r["trait"]: r for r in out_rows}
    all30 = read_tsv(ALL30)
    all30_by = {r["factor"]: r for r in all30}

    usable_trait_rgs = []
    all_avail_rgs = []
    for r in out_rows:
        if r["rg"] not in (None, ""):
            all_avail_rgs.append(float(r["rg"]))
        if r["rg_usable"] == "TRUE":
            usable_trait_rgs.append(float(r["rg"]))

    comp = []
    for fac, dom in RETAIN13:
        tr = by_trait[dom]
        of = all30_by.get(fac, {})
        try:
            f_rg = float(of["ldsc_rg"]) if of.get("ldsc_rg", "") not in ("", None) else None
        except Exception:
            f_rg = None
        try:
            f_se = float(of["ldsc_rg_se"]) if of.get("ldsc_rg_se", "") not in ("", None) else None
        except Exception:
            f_se = None
        f_usable_official = str(of.get("rg_usable", "")).strip()
        # keep official flag; also re-evaluate with same rule for transparency but do not override
        f_usable_recheck, f_why = apply_usable_rule(f_rg, f_se)
        f_lo = f_rg - 1.96 * f_se if (f_rg is not None and f_se is not None) else None
        f_hi = f_rg + 1.96 * f_se if (f_rg is not None and f_se is not None) else None
        t_rg = float(tr["rg"]) if tr["rg"] not in (None, "") else None
        t_se = float(tr["rg_se"]) if tr["rg_se"] not in (None, "") else None
        t_lo = float(tr["rg_ci95_low"]) if tr["rg_ci95_low"] not in (None, "") else None
        t_hi = float(tr["rg_ci95_high"]) if tr["rg_ci95_high"] not in (None, "") else None

        diff = (f_rg - t_rg) if (f_rg is not None and t_rg is not None) else None
        if f_rg is not None and t_rg is not None:
            if f_rg > t_rg:
                which = "factor_rg_higher"
            elif f_rg < t_rg:
                which = "dominant_trait_rg_higher"
            else:
                which = "equal"
        else:
            which = "comparison_incomplete"

        s2 = lookup_s2(fac)
        row = OrderedDict()
        row["factor"] = fac
        row["locked_class"] = of.get("locked_class", "")
        row["dominant_trait"] = dom
        row["dominant_family"] = tr["family"]
        row["note_gln_mapped_twice"] = "YES" if fac in ("F14", "F20") else "NO"
        row["official_factor_rg"] = f_rg
        row["official_factor_rg_se"] = f_se
        row["official_factor_rg_ci95_low"] = f_lo
        row["official_factor_rg_ci95_high"] = f_hi
        row["official_factor_rg_usable"] = f_usable_official
        row["official_factor_rg_usable_rule"] = of.get("rg_usable_rule", RG_USABLE_RULE)
        row["official_factor_rg_rule_recheck_pass"] = "TRUE" if f_usable_recheck else "FALSE"
        row["official_factor_rg_recheck_reason"] = f_why
        row["official_factor_gcov_int"] = of.get("ldsc_gcov_int", "")
        row["official_factor_n_snps_ldsc_rg"] = of.get("n_snps_ldsc_rg", "")
        row["official_factor_h2_ukb"] = of.get("h2_ukb", "")
        row["official_factor_h2_ukb_se"] = of.get("h2_ukb_se", "")
        row["official_factor_h2_karj"] = of.get("h2_karj", "")
        row["official_source"] = ALL30
        row["dominant_rg"] = t_rg
        row["dominant_rg_se"] = t_se
        row["dominant_rg_p"] = tr["rg_p"]
        row["dominant_rg_ci95_low"] = t_lo
        row["dominant_rg_ci95_high"] = t_hi
        row["dominant_rg_usable"] = tr["rg_usable"]
        row["dominant_unusable_or_failure_reason"] = tr["unusable_or_failure_reason"]
        row["dominant_gcov_int"] = tr["gcov_int"]
        row["dominant_gcov_int_se"] = tr["gcov_int_se"]
        row["dominant_n_snps_valid_alleles"] = tr["n_snps_valid_alleles"]
        row["dominant_ukb_h2"] = tr["ukb_h2"]
        row["dominant_ukb_h2_se"] = tr["ukb_h2_se"]
        row["dominant_karj_h2"] = tr["karj_h2"]
        row["dominant_karj_h2_se"] = tr["karj_h2_se"]
        row["karjalainen_N_dominant"] = tr["karjalainen_N"]
        row["descriptive_diff_factor_minus_dominant"] = diff
        row["which_rg_numerically_higher"] = which
        row["ci95_overlap_factor_vs_dominant"] = ci_overlap(f_lo, f_hi, t_lo, t_hi)
        row["no_independent_SE_difference_test"] = (
            "YES; shared raw GWAS; different estimands; descriptive contrast only"
        )
        row["dominant_rg_rank_among_61_available"] = ""
        row["official_factor_rg_percentile_among_61_available"] = percentile_rank(f_rg, all_avail_rgs)
        row["official_factor_rg_percentile_among_usable_single_traits"] = percentile_rank(f_rg, usable_trait_rgs)
        row["dominant_usable_but_factor_unusable"] = (
            "YES" if (tr["rg_usable"] == "TRUE" and str(f_usable_official).lower() in ("false", "0", "no")) else "NO"
        )
        row["factor_usable_but_dominant_unusable"] = (
            "YES" if (str(f_usable_official).lower() in ("true", "1", "yes") and tr["rg_usable"] != "TRUE") else "NO"
        )
        for k, v in s2.items():
            row[k] = v
        row["official_ukb_oddeven_stability"] = official_oddeven_flag(fac)
        row["selection_conditioned"] = (
            "YES; the 13 Retain factors were already filtered on Karjalainen sign and related official QC"
        )
        comp.append(row)

    # ranks of dominant rg among 61 available
    avail = [(r["trait"], float(r["rg"])) for r in out_rows if r["rg"] not in (None, "")]
    avail_sorted = sorted(avail, key=lambda x: -x[1])
    rank_map = {t: i + 1 for i, (t, _) in enumerate(avail_sorted)}
    n_avail = len(avail)
    for row in comp:
        dom = row["dominant_trait"]
        if dom in rank_map:
            row["dominant_rg_rank_among_61_available"] = f"{rank_map[dom]}/{n_avail}"

    cols13 = list(comp[0].keys())
    formatted13 = []
    for r in comp:
        formatted13.append(OrderedDict((k, fmt(r[k])) for k in cols13))
    tsv13 = os.path.join(OUT, "13_FACTOR_VS_DOMINANT_CROSS_RESOURCE.tsv")
    write_tsv(tsv13, formatted13, cols13)

    # summary numbers for the interpretation draft
    def stats(xs):
        xs = sorted(xs)
        if not xs:
            return {}
        n = len(xs)
        mid = n // 2
        med = xs[mid] if n % 2 else 0.5 * (xs[mid - 1] + xs[mid])
        return {"n": n, "min": xs[0], "max": xs[-1], "median": med, "mean": sum(xs) / n}

    n_avail_est = sum(r["estimate_available"] == "TRUE" for r in formatted)
    n_usable = sum(r["rg_usable"] == "TRUE" for r in formatted)
    usable_rgs = [float(r["rg"]) for r in formatted if r["rg_usable"] == "TRUE"]
    all_rgs = [float(r["rg"]) for r in formatted if r["rg"] != ""]
    path_sum = os.path.join(OUT, "logs/04_parse_summary.txt")
    os.makedirs(os.path.dirname(path_sum), exist_ok=True)
    with open(path_sum, "w", encoding="utf-8") as f:
        f.write(f"parsed_at\t{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} CST (Asia/Shanghai)\n")
        f.write(f"n_traits\t{len(formatted)}\n")
        f.write(f"n_estimate_available\t{n_avail_est}\n")
        f.write(f"n_rg_usable\t{n_usable}\n")
        st = stats(usable_rgs)
        f.write(f"usable_rg_stats\t{st}\n")
        st2 = stats(all_rgs)
        f.write(f"all_available_rg_stats\t{st2}\n")
        f.write(f"rg_usable_rule\t{RG_USABLE_RULE}\n")
        f.write("BH_family\t61 same-trait pairs; descriptive only; do not subset by q\n")
        f.write("\n--- 13 factor vs dominant ---\n")
        for r in formatted13:
            f.write(
                f"{r['factor']}\tf_rg={r['official_factor_rg']}\tf_se={r['official_factor_rg_se']}\t"
                f"f_usable={r['official_factor_rg_usable']}\tdom={r['dominant_trait']}\t"
                f"t_rg={r['dominant_rg']}\tt_se={r['dominant_rg_se']}\tt_usable={r['dominant_rg_usable']}\t"
                f"diff={r['descriptive_diff_factor_minus_dominant']}\twhich={r['which_rg_numerically_higher']}\t"
                f"ci_overlap={r['ci95_overlap_factor_vs_dominant']}\t"
                f"pct61={r['official_factor_rg_percentile_among_61_available']}\t"
                f"s2={r['s2_congruence_class']}|cos={r['s2_abs_cosine']}\n"
            )
        # failures
        f.write("\n--- unusable or failed single traits ---\n")
        for r in formatted:
            if r["rg_usable"] != "TRUE":
                f.write(f"{r['trait']}\trg={r['rg']}\tse={r['rg_se']}\treason={r['unusable_or_failure_reason']}\n")

    print("wrote", tsv61)
    print("wrote", tsv13)
    print("wrote", path_sum)
    print("n_available", n_avail_est, "n_usable", n_usable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
