#!/usr/bin/env Rscript
# UKB inter-factor LAVA for the 13 locked 52+9 primary factors.
# Question: do globally distinct frozen factors remain locally distinct?
# This is a confirmatory analysis of frozen factor GWAS, not a new local
# factor model, not a retain/drop rule, and not a disease analysis.
#
# Writes only under results/exp52_9/local_genetic/lava_interfactor/.
# Does not modify results/exp52_9/local_genetic/lava/ (same-factor UKB vs Karj).
# Does not write into any other project. Does not use 1000G as the LD reference.
#
# Pre-specified protocol (locked before seeing results):
# - Phenotypes: UKB projected factor GWAS only (N_UKB = 413897), 13 main factors.
# - Karjalainen inter-factor LAVA is not part of this run (optional later replication).
# - LD reference: official UKB EUR v1.1, prefix
#   supplied through LAVA_UKB_REFERENCE_PREFIX, chr1-22.
# - Loci: official blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile (2495 autosomal).
# - MHC (chr6 25-34 Mb) is flagged and excluded from primary FDR, not dropped
#   from catalogs.
# - Sample overlap: 13 x 13 matrix of 1s. All 13 factor GWAS are GLS projections
#   from the same UKB sample; off-diagonal overlap must not be set to 0.
# - Univariate filter: a pair is tested only if both factors have local h2 P < 0.05
#   in that locus. Default param.lim (+/- 1.25) is not retuned; out-of-bound rho
#   are left as LAVA sets them (NA, excluded from BH).
# - Primary FDR: BH within each factor pair, non-MHC bivariate tests.
# - Secondary FDR: BH on all non-MHC bivariate tests pooled.
# - Descriptive magnitude bins (not a decision rule for the locked set):
#     modest |rho| < 0.5; moderate 0.5 <= |rho| < 0.8; high sharing |rho| >= 0.8.
# - Interpretation of "local collapse" is relative to genome-wide inter-factor
#   LDSC rg for the same pair when that comparator is available. High local |rho|
#   for a pair that already has high genome-wide |rg| is not extra local collapse.
# - This run does not compute the LDSC comparator itself.
# - No inter-factor results will be used to drop or add factors. No factor-disease
#   LAVA and no per-locus GFA / Genomic SEM are started here.
#
# 2026-11-17 repair: the first run called run.univ() with all 13 phenotypes
#   and therefore only kept loci where every factor survived process.locus.
#   Those outputs were archived under lava_interfactor_joint13_20261117/.
#   This script now analyses whatever subset of the 13 remains in the locus.

options(stringsAsFactors = FALSE)
.libPaths(c(Sys.getenv("LAVA_R_LIBRARY", unset = file.path("external", "lava", "R_library_0.1.5")), .libPaths()))
suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
OUT  <- file.path(PROJ, "results/exp52_9/local_genetic")
LAVA_DIR <- file.path(OUT, "lava_interfactor")
SUM_DIR <- file.path(OUT, "lava_sumstats")
LOC_FILE <- Sys.getenv("LAVA_LOCUS_FILE", unset = file.path(PROJ, "external/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
REF_PREFIX <- Sys.getenv("LAVA_UKB_REFERENCE_PREFIX", unset = file.path(PROJ, "external/lava/reference_ukb_v1.1/selected/lava-ukb-v1.1"))
INFO_IN <- file.path(OUT, "lava_input.info.txt")
UNIV_THRESH <- 0.05
MHC_START <- 25000000
MHC_STOP <- 34000000
MAIN_FACTORS <- c("F9","F11","F13","F14","F15","F17","F18","F20","F22","F23","F24","F25","F26")
PHENS <- paste0("UKB_", MAIN_FACTORS)
TZ_LABEL <- "Asia/Shanghai"

dir.create(LAVA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(LAVA_DIR, "work"), recursive = TRUE, showWarnings = FALSE)

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
logmsg <- function(...) {
  cat(sprintf("[%s] %s\n", now_str(), paste0(..., collapse = "")))
  flush(stdout())
}

needed_info <- paste0(REF_PREFIX, "_chr", 1:22, ".info")
needed_bcor <- paste0(REF_PREFIX, "_chr", 1:22, ".bcor")
miss <- c(needed_info, needed_bcor)[!file.exists(c(needed_info, needed_bcor))]
if (length(miss) > 0) {
  stop("UKB v1.1 LD reference incomplete; missing: ", paste(basename(miss), collapse = ", "),
       ". Do not launch until chr1-22 are present. Do not restart existing curl/aria2 downloads from this script.")
}
if (!file.exists(LOC_FILE)) stop("Missing locus file: ", LOC_FILE)
if (!file.exists(INFO_IN)) stop("Missing lava_input.info.txt; run run_exp52_9_local_genetic.R first")

loci_all <- read.loci(LOC_FILE)
loci_all <- loci_all[loci_all$CHR %in% 1:22, ]
loci_all$mhc <- loci_all$CHR == 6 & loci_all$START < MHC_STOP & loci_all$STOP > MHC_START
logmsg("loci=", nrow(loci_all), " MHC loci=", sum(loci_all$mhc), " phenos=", paste(PHENS, collapse = ","))

info_all <- fread(INFO_IN)
info <- info_all[phenotype %in% PHENS]
if (nrow(info) != length(PHENS)) {
  stop("input.info missing UKB factors: ",
       paste(setdiff(PHENS, info$phenotype), collapse = ", "))
}
info <- info[match(PHENS, phenotype)]
if (any(!file.exists(info$filename))) {
  stop("Missing sumstats: ", paste(info$filename[!file.exists(info$filename)], collapse = ", "))
}

safe_df <- function(x) {
  if (is.null(x)) return(NULL)
  as.data.frame(x, stringsAsFactors = FALSE)
}
append_table <- function(x, file) {
  if (is.null(x) || nrow(as.data.table(x)) == 0) return(invisible(NULL))
  fwrite(as.data.table(x), file, sep = "\t", append = file.exists(file), col.names = !file.exists(file))
}
pair_id <- function(a, b) {
  aa <- sub("^UKB_", "", a)
  bb <- sub("^UKB_", "", b)
  or <- order(c(aa, bb))
  paste(c(aa, bb)[or], collapse = "|")
}

univ_file <- file.path(LAVA_DIR, "lava_univ_ukb_interfactor.tsv")
bivar_file <- file.path(LAVA_DIR, "lava_bivar_ukb_interfactor.tsv")
status_file <- file.path(LAVA_DIR, "lava_locus_status.tsv")
for (f in c(univ_file, bivar_file, status_file)) if (file.exists(f)) file.remove(f)

work <- file.path(LAVA_DIR, "work", "ukb_13")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(work, "input.info.txt")
ov_path <- file.path(work, "sample.overlap.txt")
info[, cases := NA_real_]
info[, controls := NA_real_]
write.table(as.data.frame(info[, .(phenotype, cases, controls, filename)]),
            file = info_path, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
ov <- matrix(1, nrow = length(PHENS), ncol = length(PHENS),
             dimnames = list(PHENS, PHENS))
write.table(ov, ov_path, quote = FALSE, sep = "\t", col.names = NA)

logmsg("process.input UKB 13 factors, complete sample overlap")
input <- process.input(
  input.info.file = info_path,
  sample.overlap.file = ov_path,
  ref.prefix = REF_PREFIX,
  phenos = PHENS,
  input.dir = NULL
)
logmsg("process.input done; starting locus loop")

n_locus <- nrow(loci_all)
for (i in seq_len(n_locus)) {
  loc_meta <- loci_all[i, , drop = FALSE]
  loc_id <- loc_meta$LOC
  chr <- loc_meta$CHR
  if (i == 1L || i %% 100L == 0L || i == n_locus) {
    logmsg("locus ", i, "/", n_locus, " LOC=", loc_id, " chr=", chr)
  }
  locus <- tryCatch(
    process.locus(loc_meta, input, phenos = PHENS),
    error = function(e) {
      append_table(data.table(chr = chr, locus = loc_id,
                              start = loc_meta$START, stop = loc_meta$STOP, mhc = loc_meta$mhc,
                              status = "process_locus_error", detail = conditionMessage(e)), status_file)
      NULL
    }
  )
  if (is.null(locus)) next
  # Use only phenotypes retained in this locus. Passing the full 13-factor
  # list makes run.univ() return NULL as soon as any factor was dropped
  # (negative local variance), which would restrict results to the rare
  # loci where all 13 are jointly analyzable.
  avail <- intersect(as.character(locus$phenos), PHENS)
  if (length(avail) == 0L) next
  if (i == 1L || i %% 100L == 0L || i == n_locus) {
    logmsg("locus ", i, "/", n_locus, " LOC=", loc_id, " chr=", chr,
           " n_avail=", length(avail))
  }
  univ <- tryCatch(safe_df(run.univ(locus, phenos = avail)), error = function(e) {
    append_table(data.table(chr = chr, locus = loc_id,
                            start = loc_meta$START, stop = loc_meta$STOP, mhc = loc_meta$mhc,
                            status = "run_univ_error", detail = conditionMessage(e)), status_file)
    NULL
  })
  if (!is.null(univ) && nrow(univ) > 0) {
    univ$chr <- chr
    univ$locus <- loc_id
    univ$start <- loc_meta$START
    univ$stop <- loc_meta$STOP
    univ$mhc <- loc_meta$mhc
    univ$n_available_phenos <- length(avail)
    append_table(univ, univ_file)
  }
  if (is.null(univ) || nrow(univ) == 0) next
  phen_col <- if ("phen" %in% names(univ)) "phen" else if ("phenotype" %in% names(univ)) "phenotype" else NA
  if (is.na(phen_col)) next
  keep <- univ[[phen_col]][is.finite(univ[["p"]]) & univ[["p"]] < UNIV_THRESH]
  keep <- unique(as.character(keep))
  keep <- keep[keep %in% PHENS]
  if (length(keep) < 2L) next
  b <- tryCatch(safe_df(run.bivar(locus, phenos = keep)), error = function(e) {
    append_table(data.table(chr = chr, locus = loc_id,
                            start = loc_meta$START, stop = loc_meta$STOP, mhc = loc_meta$mhc,
                            status = "run_bivar_error", detail = conditionMessage(e)), status_file)
    NULL
  })
  if (is.null(b) || nrow(b) == 0) next
  # Attach univ P for both members of each pair.
  p1_col <- if ("phen1" %in% names(b)) "phen1" else if ("phenotype1" %in% names(b)) "phenotype1" else NA
  p2_col <- if ("phen2" %in% names(b)) "phen2" else if ("phenotype2" %in% names(b)) "phenotype2" else NA
  if (!is.na(p1_col) && !is.na(p2_col)) {
    pmap <- setNames(univ[["p"]], as.character(univ[[phen_col]]))
    b$univ_p_phen1 <- as.numeric(pmap[as.character(b[[p1_col]])])
    b$univ_p_phen2 <- as.numeric(pmap[as.character(b[[p2_col]])])
    b$pair <- mapply(pair_id, as.character(b[[p1_col]]), as.character(b[[p2_col]]))
  } else {
    b$univ_p_phen1 <- NA_real_
    b$univ_p_phen2 <- NA_real_
    b$pair <- NA_character_
  }
  b$chr <- chr
  b$locus <- loc_id
  b$start <- loc_meta$START
  b$stop <- loc_meta$STOP
  b$mhc <- loc_meta$mhc
  b$n_univ_pass <- length(keep)
  b$n_available_phenos <- length(avail)
  append_table(b, bivar_file)
}

logmsg("locus loop done; computing FDR and summaries")

rho_vec <- function(dt) {
  if ("rho" %in% names(dt)) dt$rho else if ("rg" %in% names(dt)) dt$rg else rep(NA_real_, nrow(dt))
}

if (file.exists(bivar_file)) {
  bv <- fread(bivar_file)
  pcol <- if ("p" %in% names(bv)) "p" else if ("pval" %in% names(bv)) "pval" else NA
  if (!("pair" %in% names(bv)) || all(is.na(bv$pair))) {
    p1_col <- if ("phen1" %in% names(bv)) "phen1" else if ("phenotype1" %in% names(bv)) "phenotype1" else NA
    p2_col <- if ("phen2" %in% names(bv)) "phen2" else if ("phenotype2" %in% names(bv)) "phenotype2" else NA
    if (!is.na(p1_col) && !is.na(p2_col)) {
      bv[, pair := mapply(pair_id, as.character(get(p1_col)), as.character(get(p2_col)))]
    }
  }
  bv[, rho_use := rho_vec(bv)]
  if (!is.na(pcol) && "pair" %in% names(bv)) {
    bv[, fdr_within_pair_nonMHC := NA_real_]
    bv[mhc == FALSE, fdr_within_pair_nonMHC := p.adjust(get(pcol), method = "BH"), by = pair]
    nonmhc <- bv$mhc == FALSE
    bv[, fdr_pooled_nonMHC := NA_real_]
    if (any(nonmhc)) bv[nonmhc, fdr_pooled_nonMHC := p.adjust(get(pcol), method = "BH")]
    fwrite(bv, file.path(LAVA_DIR, "lava_bivar_ukb_interfactor_fdr.tsv"), sep = "\t")

    pair_summ <- bv[, {
      nm <- mhc == FALSE
      fdr <- is.finite(fdr_within_pair_nonMHC) & fdr_within_pair_nonMHC < 0.05 & nm
      r <- rho_use
      .(
        n_bivar = .N,
        n_bivar_nonMHC = sum(nm),
        n_fdr05_within_nonMHC = sum(fdr),
        n_fdr05_pooled_nonMHC = sum(nm & is.finite(fdr_pooled_nonMHC) & fdr_pooled_nonMHC < 0.05),
        median_rho_nonMHC = median(r[nm], na.rm = TRUE),
        median_abs_rho_nonMHC = median(abs(r[nm]), na.rm = TRUE),
        n_abs_rho_ge_0_5_nonMHC = sum(nm & is.finite(r) & abs(r) >= 0.5),
        n_abs_rho_ge_0_8_nonMHC = sum(nm & is.finite(r) & abs(r) >= 0.8),
        n_rho_pos_nonMHC = sum(nm & is.finite(r) & r > 0),
        n_rho_neg_nonMHC = sum(nm & is.finite(r) & r < 0),
        n_rho_na_nonMHC = sum(nm & !is.finite(r)),
        median_rho_fdr = median(r[fdr], na.rm = TRUE),
        median_abs_rho_fdr = median(abs(r[fdr]), na.rm = TRUE),
        n_abs_rho_ge_0_8_fdr = sum(fdr & is.finite(r) & abs(r) >= 0.8),
        n_abs_rho_ge_0_5_fdr = sum(fdr & is.finite(r) & abs(r) >= 0.5)
      )
    }, by = pair]
    setorder(pair_summ, pair)
    fwrite(pair_summ, file.path(LAVA_DIR, "lava_bivar_summary_by_pair.tsv"), sep = "\t")

    # Per-factor view: all non-MHC tests that include the factor.
    fac_rows <- rbindlist(lapply(MAIN_FACTORS, function(fac) {
      hit <- grepl(paste0("(^|\\|)", fac, "(\\||$)"), bv$pair)
      dt <- bv[hit]
      nm <- dt$mhc == FALSE
      r <- dt$rho_use
      fdr <- nm & is.finite(dt$fdr_within_pair_nonMHC) & dt$fdr_within_pair_nonMHC < 0.05
      data.table(
        factor = fac,
        n_bivar_nonMHC = sum(nm),
        n_fdr05_within_nonMHC = sum(fdr),
        median_abs_rho_nonMHC = median(abs(r[nm]), na.rm = TRUE),
        median_abs_rho_fdr = median(abs(r[fdr]), na.rm = TRUE),
        n_abs_rho_ge_0_8_nonMHC = sum(nm & is.finite(r) & abs(r) >= 0.8),
        n_abs_rho_ge_0_8_fdr = sum(fdr & is.finite(r) & abs(r) >= 0.8)
      )
    }))
    fwrite(fac_rows, file.path(LAVA_DIR, "lava_bivar_summary_by_factor.tsv"), sep = "\t")
  }
}

writeLines(
  c(paste0("done=", now_str(), " ", TZ_LABEL),
    paste0("univ.thresh=", UNIV_THRESH),
    "analysis=ukb_interfactor_13",
    "sample.overlap=complete"),
  file.path(LAVA_DIR, "LAVA_INTERFACTOR_DONE.txt")
)
logmsg("ALL UKB INTERFACTOR LAVA DONE")
