#!/usr/bin/env python3
"""Offline integrity and public-safety checks for the staged archive."""
from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "RELEASE_CHECK_RESULTS.tsv"
TEXT_EXTENSIONS = {".md", ".txt", ".tsv", ".csv", ".json", ".yaml", ".yml", ".cff", ".py", ".r", ".R", ".sh", ".ps1", ".gitignore", ".gitattributes"}
FORBIDDEN_PATHS = re.compile(
    r"(?:(?<![A-Za-z0-9])[A-Za-z]:[\\/]|/mnt/|/workspace/|/home/|C:\\Users\\|Administrator)",
    re.IGNORECASE,
)
FORBIDDEN_TERMS = re.compile(r"(?:GROK_BOT|BEGIN (?:RSA|OPENSSH|EC) PRIVATE KEY|api[_-]?key\s*[:=]|password\s*[:=])", re.IGNORECASE)
FORBIDDEN_SUFFIXES = {".docx", ".xlsx", ".pdf", ".rds", ".qs", ".bed", ".bim", ".fam", ".bcor"}
EXPECTED_FIGURE_HASHES = {
    "figures/Figure2_source_data.tsv": "7752b7bf849c1bf7460586c57db7976dea49146c0501b073bf94e7f2edab9ca1",
    "figures/Figure3_source_data.tsv": "17976f9370e7f3e2e2ffde55b803df113719fe64b991ff32b0aed7dce2eb0b0d",
    "figures/Figure4_source_data.tsv": "1e5b8947781d89181d07a8e4d5f8db2dd989b3fdd00d6ebaa0c9ecbf5ad579a7",
}


checks: list[dict[str, str]] = []


def add(check: str, passed: bool, observed: object, expected: object) -> None:
    checks.append(
        {
            "check": check,
            "status": "PASS" if passed else "FAIL",
            "observed": str(observed),
            "expected": str(expected),
        }
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


required = [
    "README.md",
    "LICENSE",
    "LICENSE-DATA",
    "CITATION.cff",
    ".zenodo.json",
    "config/config.example.yaml",
    "metadata/input_manifest.tsv",
    "metadata/trait_crosswalk_61.tsv",
    "derived/ukb_gfa_loadings_61x30.tsv",
    "derived/karjalainen_independent_gfa_loadings_61x15.tsv",
    "derived/factor_qc_decision.tsv",
    "figures/Figure2_source_data.tsv",
    "figures/Figure3_source_data.tsv",
    "figures/Figure4_source_data.tsv",
]
missing = [rel for rel in required if not (ROOT / rel).is_file()]
add("required files", not missing, ", ".join(missing) if missing else "all present", "all present")

files = [path for path in ROOT.rglob("*") if path.is_file() and ".git" not in path.parts]
large = [path.relative_to(ROOT).as_posix() for path in files if path.stat().st_size > 10 * 1024 * 1024]
add("files larger than 10 MiB", not large, ", ".join(large) if large else "none", "none")

forbidden_files = [path.relative_to(ROOT).as_posix() for path in files if path.suffix.lower() in {x.lower() for x in FORBIDDEN_SUFFIXES}]
add("forbidden binary/restricted file types", not forbidden_files, ", ".join(forbidden_files) if forbidden_files else "none", "none")

path_hits = []
term_hits = []
lfs_hits = []
utf8_failures = []
for path in files:
    if path == RESULTS or path.resolve() == Path(__file__).resolve():
        continue
    if path.suffix not in TEXT_EXTENSIONS and path.name not in {".gitignore", ".gitattributes", "LICENSE", "LICENSE-DATA"}:
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        utf8_failures.append(path.relative_to(ROOT).as_posix())
        continue
    rel = path.relative_to(ROOT).as_posix()
    if FORBIDDEN_PATHS.search(text):
        path_hits.append(rel)
    if FORBIDDEN_TERMS.search(text):
        term_hits.append(rel)
    if text.startswith("version https://git-lfs.github.com/spec/v1"):
        lfs_hits.append(rel)
add("UTF-8 text", not utf8_failures, ", ".join(utf8_failures) if utf8_failures else "all readable", "all readable")
add("absolute/local paths", not path_hits, ", ".join(path_hits) if path_hits else "none", "none")
add("secrets/internal bot terms", not term_hits, ", ".join(term_hits) if term_hits else "none", "none")
add("Git LFS pointer files", not lfs_hits, ", ".join(lfs_hits) if lfs_hits else "none", "none")

crosswalk = read_tsv(ROOT / "metadata" / "trait_crosswalk_61.tsv")
add("crosswalk row count", len(crosswalk) == 61, len(crosswalk), 61)
add("crosswalk unique traits", len({row["trait_code"] for row in crosswalk}) == 61, len({row["trait_code"] for row in crosswalk}), 61)

ukb_load = read_tsv(ROOT / "derived" / "ukb_gfa_loadings_61x30.tsv")
ukb_fields = list(ukb_load[0]) if ukb_load else []
factor_cols = [field for field in ukb_fields if re.fullmatch(r"F\d+", field)]
add("UKB loading dimensions", len(ukb_load) == 61 and len(factor_cols) == 30, f"{len(ukb_load)}x{len(factor_cols)}", "61x30")

karj_load = read_tsv(ROOT / "derived" / "karjalainen_independent_gfa_loadings_61x15.tsv")
karj_fields = list(karj_load[0]) if karj_load else []
karj_cols = [field for field in karj_fields if re.fullmatch(r"Karj_F\d+", field)]
add("Karjalainen loading dimensions", len(karj_load) == 61 and len(karj_cols) == 15, f"{len(karj_load)}x{len(karj_cols)}", "61x15")

qc = read_tsv(ROOT / "derived" / "factor_qc_decision.tsv")
retain = [row for row in qc if row.get("decision", "").strip().lower() == "retain"]
add("retained factor count", len(retain) == 13, len(retain), 13)

single = read_tsv(ROOT / "derived" / "single_trait_cross_resource_ldsc_61.tsv")
add("single-trait baseline rows", len(single) == 61, len(single), 61)

for rel, expected in EXPECTED_FIGURE_HASHES.items():
    observed = sha256(ROOT / rel)
    add(f"source hash {rel}", observed == expected, observed, expected)

with RESULTS.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(checks)

failed = [row for row in checks if row["status"] != "PASS"]
for row in checks:
    print(f"{row['status']}\t{row['check']}\t{row['observed']}")
print(f"SUMMARY\t{len(checks) - len(failed)} passed\t{len(failed)} failed")
sys.exit(1 if failed else 0)

