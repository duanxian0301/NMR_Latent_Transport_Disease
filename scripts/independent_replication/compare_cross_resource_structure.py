#!/usr/bin/env python3
"""Phase S2: cross-resource loading match and subspace validation.
Write only under INDEPENDENT_GFA_REPLICATION_20261119. Official paths are read-only.
Protocol: S2_INDEPENDENT_STRUCTURE_PROTOCOL.md (frozen before this run).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr, spearmanr

PROJECT_ROOT = Path(os.environ.get("NMR_GFA_PROJECT_ROOT", Path.cwd())).resolve()
ROOT = PROJECT_ROOT / "INDEPENDENT_GFA_REPLICATION_20261119"
OUT = ROOT
INP = ROOT / "s2_inputs"
S1 = ROOT / "s1_fit"
ALL30_PATH = PROJECT_ROOT / "audit_followup_20261117" / "P2_FINAL" / "ALL30_EXTERNAL_VALIDATION.tsv"

SEED = 20261119
N_NULL = 1000
CUMVAR_M = 0.99
ST_THRESH = 0.95

RETAIN13 = ["F9", "F11", "F13", "F14", "F15", "F17", "F18", "F20", "F22", "F23", "F24", "F25", "F26"]
INTERCEPT_ONLY = ["F1", "F4", "F12", "F19", "F28"]
WATCH = ["F5", "F7", "F21", "F27"]
DROP = ["F2", "F3", "F6", "F8", "F10", "F16", "F29", "F30"]
CLASS_MAP = {}
for f in RETAIN13:
    CLASS_MAP[f] = "Retain"
for f in INTERCEPT_ONLY:
    CLASS_MAP[f] = "Intercept-only"
for f in WATCH:
    CLASS_MAP[f] = "Watch"
for f in DROP:
    CLASS_MAP[f] = "Drop"

FAMILIES = {
    "Amino_acids": ["Ala", "Gln", "His", "Ile", "Leu", "Phe", "Tyr", "Val"],
    "Fatty_acids": ["DHA", "LA", "MUFA", "Omega_3", "Omega_6", "PUFA", "SFA", "Total_FA"],
    "Energy_ketone": ["Citrate", "Glucose", "Lactate", "Acetate", "Acetone", "bOHbutyrate", "Pyruvate"],
    "Inflammation_other": ["GlycA", "Albumin", "Creatinine"],
}


def locked_class(fac: str) -> str:
    return CLASS_MAP.get(fac, "unknown")


def phi_class(abs_phi: float) -> str:
    if abs_phi >= 0.95:
        return "almost_identical"
    if abs_phi >= 0.85:
        return "fair_similarity"
    if abs_phi >= 0.70:
        return "poor_to_borderline"
    return "no_adequate_correspondence"


def load_loading(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    first = df.columns[0]
    if first.lower() in {"trait", "trait_code", "unname:d: 0".lower()} or not first.startswith(("F", "Karj_")):
        df = df.set_index(first)
    else:
        # header missing index name: first column already values? or unnamed
        if df.columns[0].startswith(("F", "Karj_")):
            df = pd.read_csv(path, sep="\t", index_col=0)
        else:
            df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    # drop accidental unnamed
    df = df.loc[:, [c for c in df.columns if not c.lower().startswith("unnamed")]]
    return df.apply(pd.to_numeric, errors="coerce")


def load_R(path: Path, traits: list[str]) -> pd.DataFrame:
    R = pd.read_csv(path, sep="\t", index_col=0)
    R.index = R.index.astype(str)
    R.columns = R.columns.astype(str)
    return R.loc[traits, traits].apply(pd.to_numeric, errors="coerce")


def col_cosine(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    An = A / np.clip(np.linalg.norm(A, axis=0, keepdims=True), 1e-15, None)
    Bn = B / np.clip(np.linalg.norm(B, axis=0, keepdims=True), 1e-15, None)
    return An.T @ Bn


def pearson_mat(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    Ac = A - A.mean(axis=0, keepdims=True)
    Bc = B - B.mean(axis=0, keepdims=True)
    return col_cosine(Ac, Bc)


def spearman_mat(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    def rank_cols(X):
        # average ranks
        order = np.argsort(np.argsort(X, axis=0), axis=0)
        # handle ties poorly but loadings are effectively unique
        return order.astype(float)

    return pearson_mat(rank_cols(A), rank_cols(B))


def find_single_trait(F: np.ndarray, thresh: float = ST_THRESH) -> np.ndarray:
    Fn = F / np.clip(np.linalg.norm(F, axis=0, keepdims=True), 1e-15, None)
    return np.max(np.abs(Fn), axis=0) > thresh


def concentration(F: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ss = np.sum(F ** 2, axis=0)
    share = np.max(F ** 2, axis=0) / np.clip(ss, 1e-15, None)
    p = (F ** 2) / np.clip(ss, 1e-15, None)
    pr = 1.0 / np.clip(np.sum(p ** 2, axis=0), 1e-15, None)
    return share, pr


def dominant_traits(F: np.ndarray, names: list[str], k: int = 10) -> list[list[str]]:
    idx = np.argsort(-np.abs(F), axis=0)
    out = []
    for j in range(F.shape[1]):
        out.append([names[int(i)] for i in idx[:k, j]])
    return out


def hungarian_abs(S: np.ndarray) -> list[tuple[int, int, float]]:
    # S is n_left x n_right (Karj x UKB). maximize |S|
    cost = -np.abs(S)
    r, c = linear_sum_assignment(cost)
    return [(int(i), int(j), float(S[i, j])) for i, j in zip(r, c)]


def orthonormalize(F: np.ndarray) -> np.ndarray:
    Q, _ = np.linalg.qr(F, mode="reduced")
    # ensure full column rank numerically
    return Q[:, : F.shape[1]]


def leading_left_sv(F: np.ndarray, k: int) -> np.ndarray:
    k = min(k, F.shape[0], F.shape[1])
    U, _, _ = np.linalg.svd(F, full_matrices=False)
    return U[:, :k]


def canonical(Fa: np.ndarray, Fb: np.ndarray) -> np.ndarray:
    Qa = orthonormalize(Fa)
    Qb = orthonormalize(Fb)
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return np.clip(s, 0.0, 1.0)


def subspace_metrics(Fa: np.ndarray, Fb: np.ndarray) -> dict:
    ka, kb = Fa.shape[1], Fb.shape[1]
    rho = canonical(Fa, Fb)
    n = len(rho)
    theta = np.arccos(rho)
    sins = np.sin(theta)
    extra = abs(ka - kb)
    chordal = float(np.sqrt(np.sum(sins ** 2) + extra))
    proj_min = float(np.sum(rho ** 2) / min(ka, kb))
    proj_geom = float(np.sum(rho ** 2) / np.sqrt(ka * kb))
    # unequal-dimension Procrustes on raw loadings: C = Fa.T @ Fb
    C = Fa.T @ Fb
    s_proc = np.linalg.svd(C, compute_uv=False)
    proc_sum = float(np.sum(s_proc))
    nrm = float(np.linalg.norm(Fa) * np.linalg.norm(Fb))
    proc_cong = proc_sum / nrm if nrm > 0 else np.nan
    out = {
        "K_a": ka,
        "K_b": kb,
        "n_angles": n,
        "mean_rho": float(np.mean(rho)),
        "median_rho": float(np.median(rho)),
        "min_rho": float(np.min(rho)),
        "mean_rho_top3": float(np.mean(rho[: min(3, n)])),
        "mean_rho_top5": float(np.mean(rho[: min(5, n)])),
        "sum_rho_sq": float(np.sum(rho ** 2)),
        "projection_overlap_over_minK": proj_min,
        "projection_overlap_over_geomK": proj_geom,
        "grassmann_chordal": chordal,
        "procrustes_sum_sv": proc_sum,
        "procrustes_congruence": proc_cong,
    }
    for i, v in enumerate(rho, 1):
        out[f"rho_{i}"] = float(v)
    return out


def rand_orthonormal(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    X = rng.normal(size=(n, k))
    Q, _ = np.linalg.qr(X, mode="reduced")
    return Q[:, :k]


def eig_m(R: np.ndarray, cumvar: float = CUMVAR_M) -> tuple[int, np.ndarray, np.ndarray]:
    # R may not be exactly PD; use eigh of symmetric part
    Rs = 0.5 * (R + R.T)
    w, V = np.linalg.eigh(Rs)
    # descending
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    w_pos = np.clip(w, 0, None)
    tot = w_pos.sum()
    c = np.cumsum(w_pos) / tot
    m = int(np.searchsorted(c, cumvar) + 1)
    m = max(m, 1)
    m = min(m, R.shape[0])
    return m, w, V[:, :m]


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if (sa | sb) else np.nan


def main():
    ukb = load_loading(INP / "FULL_GFA_LOADINGS_61x30_Fhat.tsv")
    ukb_sc = load_loading(INP / "FULL_GFA_LOADINGS_61x30.tsv")
    karj = load_loading(S1 / "karjalainen_independent_F_hat.tsv")
    karj_sc = load_loading(S1 / "karjalainen_independent_F_hat_times_scale.tsv")
    odd = load_loading(INP / "F_hat_odd_seed101.tsv")
    even = load_loading(INP / "F_hat_even_seed101.tsv")

    traits = list(ukb.index)
    assert len(traits) == 61, len(traits)
    for df, name in [(karj, "karj"), (karj_sc, "karj_sc"), (ukb_sc, "ukb_sc"), (odd, "odd"), (even, "even")]:
        missing = [t for t in traits if t not in df.index]
        extra = [t for t in df.index if t not in traits]
        if missing or extra:
            raise SystemExit(f"{name} trait mismatch missing={missing} extra={extra}")
        # reorder
        # done below

    ukb = ukb.loc[traits]
    ukb_sc = ukb_sc.loc[traits]
    karj = karj.loc[traits]
    karj_sc = karj_sc.loc[traits]
    odd = odd.loc[traits]
    even = even.loc[traits]

    # families
    lipo = [t for t in traits if all(t not in v for v in FAMILIES.values())]
    families = dict(FAMILIES)
    families["Lipoproteins_lipids"] = lipo
    covered = [t for mems in families.values() for t in mems]
    assert sorted(covered) == sorted(traits), (set(traits) - set(covered), set(covered) - set(traits))
    assert len(covered) == 61

    Ruk = load_R(INP / "R_used_for_gfa_fit_seed101.tsv", traits).values
    Rkj = load_R(S1 / "R_used_for_karj_independent_gfa.tsv", traits).values
    m_ukb, ev_ukb, Vred_ukb = eig_m(Ruk)
    m_karj, ev_karj, Vred_karj = eig_m(Rkj)
    cum_ukb = float(np.clip(ev_ukb, 0, None)[:m_ukb].sum() / np.clip(ev_ukb, 0, None).sum())
    cum_karj = float(np.clip(ev_karj, 0, None)[:m_karj].sum() / np.clip(ev_karj, 0, None).sum())

    Fu = ukb.values
    Fk = karj.values
    names_u = list(ukb.columns)
    names_k = list(karj.columns)
    assert names_k == [f"Karj_F{i}" for i in range(1, 16)]
    assert names_u == [f"F{i}" for i in range(1, 31)]

    # similarities on raw F_hat (primary)
    cos = col_cosine(Fk, Fu)
    pea = pearson_mat(Fk, Fu)
    spe = spearman_mat(Fk, Fu)
    cos_sc = col_cosine(karj_sc.values, ukb_sc.values)

    share_k, pr_k = concentration(Fk)
    share_u, pr_u = concentration(Fu)
    single_k = find_single_trait(Fk)
    single_u = find_single_trait(Fu)
    top10_k = dominant_traits(Fk, traits, 10)
    top10_u = dominant_traits(Fu, traits, 10)
    top5_k = [x[:5] for x in top10_k]
    top5_u = [x[:5] for x in top10_u]
    dom_k = [x[0] for x in top10_k]
    dom_u = [x[0] for x in top10_u]
    normmax_k = np.max(np.abs(Fk / np.clip(np.linalg.norm(Fk, axis=0, keepdims=True), 1e-15, None)), axis=0)
    normmax_u = np.max(np.abs(Fu / np.clip(np.linalg.norm(Fu, axis=0, keepdims=True), 1e-15, None)), axis=0)

    # long similarity matrix
    rows = []
    for i, kn in enumerate(names_k):
        for j, un in enumerate(names_u):
            rows.append({
                "Karj_factor": kn,
                "UKB_factor": un,
                "UKB_locked_class": locked_class(un),
                "signed_cosine": float(cos[i, j]),
                "abs_cosine": float(abs(cos[i, j])),
                "signed_pearson": float(pea[i, j]),
                "abs_pearson": float(abs(pea[i, j])),
                "signed_spearman": float(spe[i, j]),
                "abs_spearman": float(abs(spe[i, j])),
                "signed_cosine_Fhat_times_scale": float(cos_sc[i, j]),
                "abs_cosine_Fhat_times_scale": float(abs(cos_sc[i, j])),
            })
    sim_df = pd.DataFrame(rows)
    sim_df.to_csv(OUT / "CROSS_RESOURCE_LOADING_SIMILARITY_MATRIX.tsv", sep="\t", index=False)

    def match_table(S: np.ndarray, right_names: list[str], tag: str) -> pd.DataFrame:
        pairs = hungarian_abs(S)
        recs = []
        used_right = set()
        for i, j, signed in pairs:
            used_right.add(j)
            abs_phi = abs(signed)
            sign_flip = -1 if signed < 0 else 1
            # signed metrics after flip
            recs.append(_pair_record(i, j, S, right_names, sign_flip, "unique_hungarian", tag))
        recs_df = pd.DataFrame(recs)
        return recs_df, used_right

    def _pair_record(i, j, S_unused, right_names, sign_flip, match_type, tag):
        un = right_names[j]
        kn = names_k[i]
        signed_cos = float(cos[i, j] * (sign_flip if right_names is names_u else np.sign(S_unused[i, j] or 1)))
        # For retain-only, use the corresponding columns in the global matrices when possible
        if un in names_u:
            jj = names_u.index(un)
            sc = float(cos[i, jj])
            sp = float(pea[i, jj])
            ss = float(spe[i, jj])
            scs = float(cos_sc[i, jj])
            t5u = top5_u[jj]
            t10u = top10_u[jj]
            du = dom_u[jj]
            su = "single" if single_u[jj] else "multi"
            shu = float(share_u[jj])
            pru = float(pr_u[jj])
            nmu = float(normmax_u[jj])
            ukb_cls = locked_class(un)
            # spearman of |loadings|
            rho_abs, _ = spearmanr(np.abs(Fk[:, i]), np.abs(Fu[:, jj]))
        else:
            raise RuntimeError(un)
        flip = -1 if sc < 0 else 1
        return {
            "match_set": tag,
            "match_type": match_type,
            "Karj_factor": kn,
            "UKB_factor": un,
            "UKB_locked_class": ukb_cls,
            "sign_flip_applied_to_Karj": flip,
            "signed_cosine": sc * flip,
            "abs_cosine": abs(sc),
            "congruence_class": phi_class(abs(sc)),
            "signed_pearson": sp * flip,
            "abs_pearson": abs(sp),
            "signed_spearman": ss * flip,
            "abs_spearman": abs(ss),
            "abs_cosine_Fhat_times_scale": abs(scs),
            "Karj_single_trait_0p95": "single" if single_k[i] else "multi",
            "UKB_single_trait_0p95": su,
            "Karj_norm_col_max": float(normmax_k[i]),
            "UKB_norm_col_max": nmu,
            "Karj_concentration_share": float(share_k[i]),
            "UKB_concentration_share": shu,
            "Karj_participation_ratio": float(pr_k[i]),
            "UKB_participation_ratio": pru,
            "Karj_dominant_trait": dom_k[i],
            "UKB_dominant_trait": du,
            "same_dominant_trait": bool(dom_k[i] == du),
            "top5_jaccard": jaccard(top5_k[i], t5u),
            "top10_jaccard": jaccard(top10_k[i], t10u),
            "spearman_abs_loadings_all61": float(rho_abs),
            "Karj_top5": ",".join(top5_k[i]),
            "UKB_top5": ",".join(t5u),
        }

    # unique hungarian on 15 x 30
    match_u, used_u = match_table(cos, names_u, "all30")
    # best (non-unique)
    best_rows = []
    for i in range(len(names_k)):
        j = int(np.argmax(np.abs(cos[i, :])))
        best_rows.append(_pair_record(i, j, cos, names_u, 1, "best_nonunique", "all30"))
    best_df = pd.DataFrame(best_rows)

    # retain 13 only
    ret_idx = [names_u.index(f) for f in RETAIN13]
    cos_ret = cos[:, ret_idx]
    pairs_ret = hungarian_abs(cos_ret)
    ret_rows = []
    for i, jloc, _ in pairs_ret:
        j = ret_idx[jloc]
        ret_rows.append(_pair_record(i, j, cos, names_u, 1, "unique_hungarian", "retain13_selection_conditioned"))
    # two Karj unmatched in unique retain match
    matched_k_ret = {i for i, _, _ in pairs_ret}
    for i in range(len(names_k)):
        if i not in matched_k_ret:
            j = ret_idx[int(np.argmax(np.abs(cos_ret[i, :])))]
            rec = _pair_record(i, j, cos, names_u, 1, "unmatched_best_in_retain13", "retain13_selection_conditioned")
            ret_rows.append(rec)
    ret_df = pd.DataFrame(ret_rows)

    match_all = pd.concat([match_u, best_df], ignore_index=True)
    match_all.to_csv(OUT / "CROSS_RESOURCE_LOADING_MATCH.tsv", sep="\t", index=False)
    ret_df.to_csv(OUT / "CROSS_RESOURCE_LOADING_MATCH_RETAIN13.tsv", sep="\t", index=False)

    # concentration table
    conc_rows = []
    for i, kn in enumerate(names_k):
        conc_rows.append({
            "side": "Karjalainen_independent",
            "factor": kn,
            "locked_class": "NA",
            "single_trait_0p95": "single" if single_k[i] else "multi",
            "norm_col_max": float(normmax_k[i]),
            "concentration_share": float(share_k[i]),
            "participation_ratio": float(pr_k[i]),
            "dominant_trait": dom_k[i],
            "top5": ",".join(top5_k[i]),
        })
    for j, un in enumerate(names_u):
        conc_rows.append({
            "side": "UKB_official_seed101",
            "factor": un,
            "locked_class": locked_class(un),
            "single_trait_0p95": "single" if single_u[j] else "multi",
            "norm_col_max": float(normmax_u[j]),
            "concentration_share": float(share_u[j]),
            "participation_ratio": float(pr_u[j]),
            "dominant_trait": dom_u[j],
            "top5": ",".join(top5_u[j]),
        })
    pd.DataFrame(conc_rows).to_csv(OUT / "S2_FACTOR_CONCENTRATION.tsv", sep="\t", index=False)

    # ---- subspace ----
    spaces = []
    spaces.append(("UKB_K30_vs_Karj_K15", "primary", Fu, Fk))
    for k in (3, 5, 10, 13):
        spaces.append((
            f"UKB_leadingsv{k}_vs_Karj_leadingsv{min(k, 15)}",
            "pre_frozen_ukb_only_sv",
            leading_left_sv(Fu, k),
            leading_left_sv(Fk, min(k, 15)),
        ))
    spaces.append(("UKB_Retain13_vs_Karj_K15", "selection_conditioned", Fu[:, ret_idx], Fk))
    spaces.append(("UKB_Retain13_vs_Karj_leadingsv13", "selection_conditioned", Fu[:, ret_idx], leading_left_sv(Fk, 13)))

    obs_metrics = []
    for label, kind, A, B in spaces:
        met = subspace_metrics(A, B)
        met["comparison"] = label
        met["comparison_kind"] = kind
        obs_metrics.append(met)
    # primary A,B for nulls
    A_pri, B_pri = Fu, Fk
    obs_pri = [m for m in obs_metrics if m["comparison"] == "UKB_K30_vs_Karj_K15"][0]

    rng = np.random.default_rng(SEED)
    n_trait = 61
    k_b = Fk.shape[1]

    def summarize_null(draws: list[dict], prefix: str, observed: dict) -> dict:
        keys = ["mean_rho", "mean_rho_top3", "mean_rho_top5", "min_rho", "projection_overlap_over_minK", "grassmann_chordal", "procrustes_congruence"]
        out = {}
        for key in keys:
            arr = np.array([d[key] for d in draws], dtype=float)
            obs = observed[key]
            # for distances (chordal), excess is SMALLER than null; for similarities, LARGER
            if key == "grassmann_chordal":
                p = (1 + int(np.sum(arr <= obs))) / (1 + len(arr))
                tail = "left"
            else:
                p = (1 + int(np.sum(arr >= obs))) / (1 + len(arr))
                tail = "right"
            out[f"{prefix}_{key}_obs"] = obs
            out[f"{prefix}_{key}_null_mean"] = float(np.mean(arr))
            out[f"{prefix}_{key}_null_sd"] = float(np.std(arr, ddof=1))
            out[f"{prefix}_{key}_null_p05"] = float(np.quantile(arr, 0.05 if tail == "left" else 0.95))
            out[f"{prefix}_{key}_p"] = float(p)
            out[f"{prefix}_{key}_tail"] = tail
        return out

    null_weak = []
    null_red_ukb = []
    null_red_karj = []
    for b in range(N_NULL):
        Qw = rand_orthonormal(n_trait, k_b, rng)
        null_weak.append(subspace_metrics(A_pri, Qw))
        # random frame in redundancy subspace: Vred @ rand k
        Z = rng.normal(size=(m_ukb, k_b))
        Qr, _ = np.linalg.qr(Z, mode="reduced")
        Fr = Vred_ukb @ Qr[:, :k_b]
        null_red_ukb.append(subspace_metrics(A_pri, Fr))
        Z2 = rng.normal(size=(m_karj, k_b))
        Qr2, _ = np.linalg.qr(Z2, mode="reduced")
        Fr2 = Vred_karj @ Qr2[:, :k_b]
        null_red_karj.append(subspace_metrics(A_pri, Fr2))

    # permutation nulls of Karj loadings
    # family-internal
    family_index = {fam: [traits.index(t) for t in mems] for fam, mems in families.items()}

    def permute_within_families(F, rng):
        Fp = F.copy()
        for idxs in family_index.values():
            perm = rng.permutation(len(idxs))
            Fp[idxs, :] = F[np.array(idxs)[perm], :]
        return Fp

    def permute_global(F, rng):
        perm = rng.permutation(F.shape[0])
        return F[perm, :]

    null_fam = []
    null_glob = []
    for b in range(N_NULL):
        null_fam.append(subspace_metrics(A_pri, permute_within_families(B_pri, rng)))
        null_glob.append(subspace_metrics(A_pri, permute_global(B_pri, rng)))

    cal = {
        "seed": SEED,
        "B": N_NULL,
        "m_rule": "smallest m with cumulative eigenvalue proportion >= 0.99",
        "m_ukb_R": int(m_ukb),
        "m_ukb_cumvar": cum_ukb,
        "m_karj_R": int(m_karj),
        "m_karj_cumvar": cum_karj,
        "family_sizes": {k: len(v) for k, v in families.items()},
    }
    cal.update(summarize_null(null_weak, "weak", obs_pri))
    cal.update(summarize_null(null_red_ukb, "redundancy_ukbR", obs_pri))
    cal.update(summarize_null(null_red_karj, "redundancy_karjR", obs_pri))
    cal.update(summarize_null(null_fam, "family_perm", obs_pri))
    cal.update(summarize_null(null_glob, "global_perm", obs_pri))

    # write subspace table: observed comparisons + primary null summary columns
    sub_df = pd.DataFrame(obs_metrics)
    # attach primary null p-values as extra rows? better a second file for calibration numbers
    # put p-values on the primary row by merge
    for key, val in cal.items():
        if key.startswith(("weak_", "redundancy_", "family_", "global_")):
            sub_df.loc[sub_df["comparison"] == "UKB_K30_vs_Karj_K15", key] = val
    sub_df.to_csv(OUT / "CROSS_RESOURCE_SUBSPACE_VALIDATION.tsv", sep="\t", index=False)
    # compact calibration json/tsv
    cal_rows = [{"key": k, "value": (json.dumps(v) if isinstance(v, (dict, list)) else v)} for k, v in cal.items()]
    pd.DataFrame(cal_rows).to_csv(OUT / "CROSS_RESOURCE_NULL_SUMMARY.tsv", sep="\t", index=False)

    # ---- SNP-split stability (UKB reuse) ----
    def align_to_ref(Fnew: np.ndarray, names_new: list[str], Fref: np.ndarray, names_ref: list[str], tag: str) -> pd.DataFrame:
        S = col_cosine(Fnew, Fref)
        pairs = hungarian_abs(S)
        recs = []
        for i, j, signed in pairs:
            abs_phi = abs(float(signed))
            recs.append({
                "comparison": tag,
                "match_type": "unique_hungarian",
                "left_factor": names_new[i],
                "right_factor": names_ref[j],
                "right_locked_class": locked_class(names_ref[j]) if names_ref[j] in CLASS_MAP else "NA",
                "abs_cosine": abs_phi,
                "signed_cosine_after_flip": abs_phi,  # after flip, signed = abs
                "congruence_class": phi_class(abs_phi),
                "left_single_trait_0p95": "single" if find_single_trait(Fnew)[i] else "multi",
                "right_single_trait_0p95": "single" if find_single_trait(Fref)[j] else "multi",
                "left_dominant": dominant_traits(Fnew, traits, 1)[i][0],
                "right_dominant": dominant_traits(Fref, traits, 1)[j][0],
                "same_dominant_trait": dominant_traits(Fnew, traits, 1)[i][0] == dominant_traits(Fref, traits, 1)[j][0],
                "left_concentration_share": float(concentration(Fnew)[0][i]),
                "right_concentration_share": float(concentration(Fref)[0][j]),
            })
        return pd.DataFrame(recs)

    stab = pd.concat([
        align_to_ref(odd.values, list(odd.columns), Fu, names_u, "UKB_odd_vs_official101"),
        align_to_ref(even.values, list(even.columns), Fu, names_u, "UKB_even_vs_official101"),
        align_to_ref(Fk, names_k, odd.values, list(odd.columns), "Karj15_vs_UKB_odd"),
        align_to_ref(Fk, names_k, even.values, list(even.columns), "Karj15_vs_UKB_even"),
    ], ignore_index=True)
    stab.to_csv(OUT / "S2_STABILITY_MATCHES.tsv", sep="\t", index=False)

    # UKB factors with both odd and even unique match |phi|>=0.85
    odd_m = stab[stab["comparison"] == "UKB_odd_vs_official101"]
    even_m = stab[stab["comparison"] == "UKB_even_vs_official101"]
    odd_ok = set(odd_m.loc[odd_m["abs_cosine"] >= 0.85, "right_factor"])
    even_ok = set(even_m.loc[even_m["abs_cosine"] >= 0.85, "right_factor"])
    both_ok = sorted(odd_ok & even_ok, key=lambda x: int(x[1:]))

    # Karj unique matches vs official
    karj_unique = match_u.copy()
    fair_multi = karj_unique[(karj_unique["abs_cosine"] >= 0.85) & (karj_unique["Karj_single_trait_0p95"] == "multi")]
    n_fair_multi = int(len(fair_multi))
    n_fair_any = int((karj_unique["abs_cosine"] >= 0.85).sum())
    n_ident_any = int((karj_unique["abs_cosine"] >= 0.95).sum())

    p_mean_rho_red = float(cal["redundancy_ukbR_mean_rho_p"])
    gate_A = n_fair_multi >= 3
    gate_B = p_mean_rho_red < 0.05
    gate_pass = bool(gate_A or gate_B)

    gate_txt = []
    gate_txt.append("S2 Karjalainen SNP-split gate")
    gate_txt.append(f"evaluated {pd.Timestamp.now().tz_localize('Asia/Shanghai')}")
    gate_txt.append(f"n_unique_matches_abs_cosine_ge_0.95 = {n_ident_any}")
    gate_txt.append(f"n_unique_matches_abs_cosine_ge_0.85 = {n_fair_any}")
    gate_txt.append(f"n_unique_matches_abs_cosine_ge_0.85_and_Karj_multi = {n_fair_multi}")
    gate_txt.append(f"gate_A (>=3 multi-trait unique matches at 0.85) = {gate_A}")
    gate_txt.append(f"primary_null = redundancy_ukbR mean_rho p = {p_mean_rho_red}")
    gate_txt.append(f"gate_B (p < 0.05) = {gate_B}")
    gate_txt.append(f"gate_pass = {gate_pass}")
    gate_txt.append("Karjalainen_SNP_split = " + ("AUTHORIZED_BY_GATE" if gate_pass else "NOT_RUN"))
    gate_txt.append("UKB_odd_even_both_abs_cosine_ge_0.85 = " + ",".join(both_ok) if both_ok else "UKB_odd_even_both_abs_cosine_ge_0.85 = none")
    (OUT / "S2_KARJ_SNPSPLIT_GATE.txt").write_text("\n".join(gate_txt) + "\n", encoding="utf-8")

    # juxtaposition extract from official ALL30 (read-only)
    if ALL30_PATH.exists():
        all30 = pd.read_csv(ALL30_PATH, sep="\t")
        keep = all30.copy()
        keep.to_csv(OUT / "S2_OFFICIAL_ALL30_JUXTAPOSITION_COPY.tsv", sep="\t", index=False)
        retain = all30[all30["factor"].isin(RETAIN13)].copy()
        n_usable = int(retain["rg_usable"].astype(str).isin(["True", "TRUE", "true", "1"]).sum()) if "rg_usable" in retain.columns else np.nan
        usable_factors = retain.loc[retain["rg_usable"].astype(str).isin(["True", "TRUE", "true", "1"]), "factor"].tolist() if "rg_usable" in retain.columns else []
        (OUT / "S2_OFFICIAL_RETAIN13_USABLE_RG.txt").write_text(
            f"usable_rg {n_usable}/13\nfactors {','.join(usable_factors)}\n",
            encoding="utf-8",
        )

    # machine-readable run log
    log = {
        "n_trait": 61,
        "UKB_K": 30,
        "Karj_K": 15,
        "odd_K": int(odd.shape[1]),
        "even_K": int(even.shape[1]),
        "m_ukb_R": int(m_ukb),
        "m_karj_R": int(m_karj),
        "n_fair_multi": n_fair_multi,
        "n_fair_any": n_fair_any,
        "n_ident_any": n_ident_any,
        "gate_A": gate_A,
        "gate_B": gate_B,
        "gate_pass": gate_pass,
        "p_mean_rho_redundancy_ukbR": p_mean_rho_red,
        "obs_mean_rho_primary": obs_pri["mean_rho"],
        "UKB_split_stable_both_0.85": both_ok,
        "families": {k: v for k, v in families.items()},
    }
    (OUT / "S2_RUN_LOG.json").write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(log, indent=2, ensure_ascii=False))
    print("WROTE S2 tables under", OUT)


if __name__ == "__main__":
    main()
