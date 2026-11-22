#!/usr/bin/env Rscript
# 52+9 main-retained factor locus-level cross-cohort replication.
# Writes ONLY under results/exp52_9/local_genetic/.
# Does not overwrite gls_Rreg, factor_qc, GFA, week1, or inputs.
# Does not run GLS, GFA, GLEANR, disease F-GWAS/MR, or cell-level GFA.

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
EXP  <- file.path(PROJ, "results/exp52_9")
GLS  <- file.path(EXP, "gls_Rreg")
OUT  <- file.path(EXP, "local_genetic")
LAVA_OUT <- file.path(OUT, "lava_sumstats")

UK_FG <- file.path(GLS, "fgwas_ukbb_eur_gls_Rreg_genomewide.tsv.gz")
KJ_FG <- file.path(GLS, "fgwas_karjalainen_gls_Rreg_genomewide.tsv.gz")
MAP_PATH <- file.path(EXP, "cache/snp_map.qs")
OLD_LEADS <- file.path(GLS, "cross_cohort_leads_gls_Rreg.tsv")
QC_DEC <- file.path(EXP, "factor_qc/factor_qc_decision.tsv")

N_UKB <- 413897
N_KARJ <- 136016
MHC_CHR <- 6L
MHC_START <- 25e6
MHC_END <- 34e6
CLUMP_KB_MAIN <- 1000
WIN_BP <- 1000L * 1000L
Z_UKB_MAIN <- 8
Z_KARJ_MAIN <- 5
TZ_LABEL <- "Asia/Shanghai"

MAIN_FACTORS <- c("F9","F11","F13","F14","F15","F17","F18","F20","F22","F23","F24","F25","F26")
SUPP_FACTORS <- c("F1","F4","F12","F19","F28")
ALL_FACTORS <- c(MAIN_FACTORS, SUPP_FACTORS)

STATUS <- file.path(OUT, "status.md")
REPORT <- file.path(OUT, "LOCAL_GENETIC_REPORT.md")
DONE_FLAG <- file.path(OUT, "LOCAL_GENETIC_DONE.txt")
PIDFILE <- file.path(OUT, "local_genetic.pid")

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")

logmsg <- function(...) {
  msg <- sprintf("[%s] %s\n", now_str(), paste0(..., collapse = ""))
  cat(msg)
  flush(stdout())
  try(flush.console(), silent = TRUE)
}

file_ok <- function(path, min_bytes = 50) {
  file.exists(path) && isTRUE(file.info(path)$size >= min_bytes)
}

assert_write_path <- function(path) {
  allowed <- normalizePath(OUT, winslash = "/", mustWork = FALSE)
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
  path_n <- normalizePath(path, winslash = "/", mustWork = FALSE)
  ok <- startsWith(tolower(parent), tolower(allowed)) ||
    startsWith(tolower(path_n), tolower(allowed))
  if (!isTRUE(ok)) stop("Refusing to write outside results/exp52_9/local_genetic/: ", path)
  invisible(TRUE)
}

write_status <- function(step, detail = "") {
  assert_write_path(STATUS)
  lines <- c(
    "# 52+9 locus-level cross-cohort replication — live status",
    "",
    paste0("Updated: ", now_str(), " (", TZ_LABEL, ")"),
    "",
    paste0("**Step:** ", step),
    "",
    detail,
    "",
    "## Policy",
    "- Writes only under `results/exp52_9/local_genetic/`.",
    "- Main analysis: 13 intercept-passing factors. Supplement (not main): F1, F4, F12, F19, F28.",
    "- Watch set F5, F7, F21, F27 is not tabulated.",
    "- Existing gls_Rreg / factor_qc / GFA / week1 / inputs are read-only.",
    "- MHC chr6 25–34 Mb is flagged, not dropped from lead catalogs (old lead protocol also did not drop MHC).",
    "- No GLS/GFA rerun, no GLEANR, no disease F-GWAS/MR, no cell-level GFA."
  )
  writeLines(lines, STATUS, useBytes = TRUE)
}

# Identical to results/exp52_9 gls_Rreg and week1 Step 2 distance clumping.
clump_distance <- function(chr, bp, z, kb = 1000) {
  n <- length(z)
  keep <- rep(FALSE, n)
  ok <- is.finite(chr) & is.finite(bp) & is.finite(z) & is.finite(kb)
  if (!any(ok)) return(keep)
  win <- as.numeric(kb) * 1000
  ord <- which(ok)[order(abs(z[ok]), decreasing = TRUE)]
  occupied <- new.env(hash = TRUE)
  for (i in ord) {
    ckey <- as.character(chr[i])
    b <- bp[i]
    b0 <- floor(b / win)
    blocked <- FALSE
    for (db in (b0 - 1):(b0 + 1)) {
      k <- paste(ckey, db, sep = ":")
      if (exists(k, envir = occupied, inherits = FALSE)) {
        if (abs(occupied[[k]] - b) <= win) { blocked <- TRUE; break }
      }
    }
    if (blocked) next
    keep[i] <- TRUE
    assign(paste(ckey, b0, sep = ":"), b, envir = occupied)
  }
  keep
}

in_mhc <- function(chr, bp) {
  !is.na(chr) & !is.na(bp) & chr == MHC_CHR & bp >= MHC_START & bp <= MHC_END
}

p_from_z <- function(z) 2 * pnorm(-abs(z))

fmt_pct <- function(x, d = 3) {
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  if (any(ok)) out[ok] <- sprintf(paste0("%.", d, "f"), x[ok])
  out
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(LAVA_OUT, recursive = TRUE, showWarnings = FALSE)
assert_write_path(file.path(OUT, "dummy"))
cat(as.integer(Sys.getpid()), file = PIDFILE)

write_status("START", "Reading inputs and preparing lead catalogs.")
logmsg("pid=", Sys.getpid(), " start local-genetic analysis")
logmsg("OUT=", OUT)

stopifnot(file_ok(UK_FG, 1e6), file_ok(KJ_FG, 1e6), file_ok(MAP_PATH, 1e3), file_ok(OLD_LEADS))

map <- qread(MAP_PATH)
setDT(map)
need_map <- c("SNP", "CHR", "BP", "A1", "A2")
if (!all(need_map %in% names(map))) stop("snp_map missing columns: ", paste(setdiff(need_map, names(map)), collapse = ","))
map <- unique(map[, ..need_map], by = "SNP")
logmsg("snp_map rows=", nrow(map))

zcols <- paste0("Z_", ALL_FACTORS)
logmsg("Reading UKB factor GWAS (selected Z columns)...")
uk <- fread(UK_FG, select = c("SNP", zcols))
logmsg("UKB rows=", nrow(uk), " cols=", ncol(uk))
logmsg("Reading Karjalainen factor GWAS (selected Z columns)...")
kj <- fread(KJ_FG, select = c("SNP", zcols))
logmsg("Karj rows=", nrow(kj), " cols=", ncol(kj))

setnames(uk, zcols, paste0("ZU_", ALL_FACTORS))
setnames(kj, zcols, paste0("ZK_", ALL_FACTORS))

dt <- merge(uk, kj, by = "SNP")
rm(uk, kj)
gc(verbose = FALSE)
dt <- merge(dt, map, by = "SNP", all.x = TRUE)
logmsg("merged rows=", nrow(dt), " missing CHR=", sum(is.na(dt$CHR)))
if (nrow(dt) != 1202568L) logmsg("WARNING: merged SNP count is ", nrow(dt), " expected 1202568")

old_leads <- fread(OLD_LEADS)
old8 <- old_leads[lead_set == "ukb_clumped_1mb_zge8"]

ukb_rows <- list()
karj5_rows <- list()
karj_p5_rows <- list()
pair_rows <- list()
sum_rows <- list()
val_rows <- list()

write_status("CLUMP", "Distance-clumping UKB |Z|>=8 and Karj |Z|>=5 (1 Mb), then overlapping loci.")

for (fac in ALL_FACTORS) {
  aset <- if (fac %in% MAIN_FACTORS) "main" else "supplement_intercept_insufficient"
  zu <- dt[[paste0("ZU_", fac)]]
  zk <- dt[[paste0("ZK_", fac)]]
  ok <- is.finite(zu) & is.finite(zk) & is.finite(dt$CHR) & is.finite(dt$BP)
  cl_u <- clump_distance(dt$CHR, dt$BP, zu, kb = CLUMP_KB_MAIN)
  cl_k <- clump_distance(dt$CHR, dt$BP, zk, kb = CLUMP_KB_MAIN)
  uk_mask <- cl_u & ok & abs(zu) >= Z_UKB_MAIN
  kj_mask <- cl_k & ok & abs(zk) >= Z_KARJ_MAIN
  kj_p5_mask <- cl_k & ok & p_from_z(zk) < 5e-8

  uk_leads <- dt[uk_mask, .(
    factor = fac, analysis_set = aset, SNP, CHR, BP, A1, A2,
    Z_UKB = zu[uk_mask], Z_KARJ = zk[uk_mask]
  )]
  if (nrow(uk_leads) > 0) {
    uk_leads[, `:=`(
      sign_concord = as.integer(sign(Z_UKB) == sign(Z_KARJ)),
      p_UKB = p_from_z(Z_UKB),
      p_KARJ = p_from_z(Z_KARJ),
      chi2_UKB = Z_UKB^2,
      in_mhc = in_mhc(CHR, BP),
      clump_kb = CLUMP_KB_MAIN,
      z_threshold = Z_UKB_MAIN
    )]
    setorder(uk_leads, CHR, BP)
  }

  kj_leads <- dt[kj_mask, .(
    factor = fac, analysis_set = aset, SNP, CHR, BP, A1, A2,
    Z_KARJ = zk[kj_mask], Z_UKB = zu[kj_mask]
  )]
  if (nrow(kj_leads) > 0) {
    kj_leads[, `:=`(
      sign_concord = as.integer(sign(Z_UKB) == sign(Z_KARJ)),
      p_KARJ = p_from_z(Z_KARJ),
      p_UKB = p_from_z(Z_UKB),
      chi2_KARJ = Z_KARJ^2,
      in_mhc = in_mhc(CHR, BP),
      clump_kb = CLUMP_KB_MAIN,
      z_threshold = Z_KARJ_MAIN,
      pass_p5e8 = p_from_z(Z_KARJ) < 5e-8
    )]
    setorder(kj_leads, CHR, BP)
  }

  kj_p5 <- dt[kj_p5_mask, .(
    factor = fac, analysis_set = aset, SNP, CHR, BP, A1, A2,
    Z_KARJ = zk[kj_p5_mask], Z_UKB = zu[kj_p5_mask]
  )]
  if (nrow(kj_p5) > 0) {
    kj_p5[, `:=`(
      sign_concord = as.integer(sign(Z_UKB) == sign(Z_KARJ)),
      p_KARJ = p_from_z(Z_KARJ),
      p_UKB = p_from_z(Z_UKB),
      chi2_KARJ = Z_KARJ^2,
      in_mhc = in_mhc(CHR, BP),
      clump_kb = CLUMP_KB_MAIN,
      p_threshold = 5e-8
    )]
    setorder(kj_p5, CHR, BP)
  }

  pairs <- data.table()
  if (nrow(uk_leads) > 0 && nrow(kj_leads) > 0) {
    u <- uk_leads[, .(ukb_SNP = SNP, ukb_chr = CHR, ukb_bp = BP, ukb_z = Z_UKB, ukb_in_mhc = in_mhc)]
    k <- kj_leads[, .(karj_SNP = SNP, karj_chr = CHR, karj_bp = BP, karj_z = Z_KARJ, karj_in_mhc = in_mhc)]
    u[, ukb_id := .I]
    k[, karj_id := .I]
    pr <- list()
    for (chr in intersect(unique(u$ukb_chr), unique(k$karj_chr))) {
      uu <- u[ukb_chr == chr]
      kk <- k[karj_chr == chr]
      if (nrow(uu) == 0 || nrow(kk) == 0) next
      for (i in seq_len(nrow(uu))) {
        d <- abs(kk$karj_bp - uu$ukb_bp[i])
        hit <- which(d <= WIN_BP)
        if (length(hit) == 0) next
        pr[[length(pr) + 1L]] <- data.table(
          factor = fac, analysis_set = aset,
          ukb_id = uu$ukb_id[i], karj_id = kk$karj_id[hit],
          ukb_SNP = uu$ukb_SNP[i], ukb_chr = chr, ukb_bp = uu$ukb_bp[i], ukb_z = uu$ukb_z[i],
          karj_SNP = kk$karj_SNP[hit], karj_chr = chr, karj_bp = kk$karj_bp[hit], karj_z = kk$karj_z[hit],
          distance_bp = d[hit],
          ukb_in_mhc = uu$ukb_in_mhc[i], karj_in_mhc = kk$karj_in_mhc[hit]
        )
      }
    }
    if (length(pr) > 0) {
      pairs <- rbindlist(pr)
      pairs[, same_sign := as.integer(sign(ukb_z) == sign(karj_z))]
      pairs[, absd := distance_bp]
      pairs[, is_closest_for_ukb := absd == min(absd), by = ukb_id]
      pairs[, is_closest_for_karj := absd == min(absd), by = karj_id]
      # tie-break: if two equally close, keep the first after sorting by distance then SNP
      setorder(pairs, ukb_id, distance_bp, karj_SNP)
      pairs[, is_closest_for_ukb := seq_len(.N) == 1L, by = ukb_id]
      setorder(pairs, karj_id, distance_bp, ukb_SNP)
      pairs[, is_closest_for_karj := seq_len(.N) == 1L, by = karj_id]
      setorder(pairs, ukb_chr, ukb_bp, distance_bp)
      pairs[, absd := NULL]
    }
  }

  n_ukb <- nrow(uk_leads)
  n_karj <- nrow(kj_leads)
  n_karj_p5 <- nrow(kj_p5)
  ukb_overlap_ids <- if (nrow(pairs) > 0) unique(pairs$ukb_id) else integer(0)
  karj_overlap_ids <- if (nrow(pairs) > 0) unique(pairs$karj_id) else integer(0)
  n_ukb_ov <- length(ukb_overlap_ids)
  n_karj_ov <- length(karj_overlap_ids)
  closest <- if (nrow(pairs) > 0) pairs[is_closest_for_ukb == TRUE] else pairs
  n_closest <- nrow(closest)
  n_closest_same <- if (n_closest > 0) sum(closest$same_sign == 1L) else 0L
  sign_closest <- if (n_closest > 0) mean(closest$same_sign) else NA_real_

  chi_u <- if (n_ukb > 0) uk_leads$chi2_UKB else numeric(0)
  chi_k <- if (n_karj > 0) kj_leads$chi2_KARJ else numeric(0)
  sum_chi_u <- sum(chi_u)
  sum_chi_k <- sum(chi_k)
  max_share_u <- if (sum_chi_u > 0) max(chi_u) / sum_chi_u else NA_real_
  top3_share_u <- if (sum_chi_u > 0) sum(sort(chi_u, decreasing = TRUE)[seq_len(min(3L, length(chi_u)))]) / sum_chi_u else NA_real_
  max_share_k <- if (sum_chi_k > 0) max(chi_k) / sum_chi_k else NA_real_

  sign_u <- if (n_ukb > 0) mean(uk_leads$sign_concord) else NA_real_
  n_same_u <- if (n_ukb > 0) sum(uk_leads$sign_concord == 1L) else 0L

  old <- old8[factor == fac]
  val_n <- if (nrow(old) == 1L) as.integer(old$n_leads[1]) == n_ukb else NA
  val_s <- if (nrow(old) == 1L && n_ukb > 0 && is.finite(old$sign_concord[1])) {
    isTRUE(all.equal(as.numeric(old$sign_concord[1]), sign_u, tolerance = 1e-8, scale = 1))
  } else if (nrow(old) == 1L && n_ukb == 0L && is.na(old$sign_concord[1])) {
    TRUE
  } else {
    NA
  }

  sum_rows[[fac]] <- data.table(
    factor = fac,
    analysis_set = aset,
    n_loci_ukb_zge8 = n_ukb,
    n_ukb_mhc = if (n_ukb > 0) sum(uk_leads$in_mhc) else 0L,
    n_same_sign_ukb_leads = n_same_u,
    sign_concord_ukb_leads = sign_u,
    n_ukb_leads_karj_p05 = if (n_ukb > 0) sum(uk_leads$p_KARJ < 0.05) else 0L,
    n_ukb_leads_karj_p1e3 = if (n_ukb > 0) sum(uk_leads$p_KARJ < 1e-3) else 0L,
    n_ukb_leads_karj_p5e8 = if (n_ukb > 0) sum(uk_leads$p_KARJ < 5e-8) else 0L,
    n_loci_karj_zge5 = n_karj,
    n_karj_mhc = if (n_karj > 0) sum(kj_leads$in_mhc) else 0L,
    n_loci_karj_p5e8 = n_karj_p5,
    n_ukb_overlap = n_ukb_ov,
    n_ukb_only = n_ukb - n_ukb_ov,
    n_karj_overlap = n_karj_ov,
    n_karj_only = n_karj - n_karj_ov,
    n_overlap_pairs = nrow(pairs),
    n_closest_pairs = n_closest,
    n_closest_same_sign = n_closest_same,
    sign_concord_closest_pairs = sign_closest,
    max_chi2_share_ukb = max_share_u,
    top3_chi2_share_ukb = top3_share_u,
    max_chi2_share_karj = max_share_k,
    validation_n_ukb_match = val_n,
    validation_sign_match = val_s,
    old_n_leads_zge8 = if (nrow(old) == 1L) as.integer(old$n_leads[1]) else NA_integer_,
    old_sign_concord_zge8 = if (nrow(old) == 1L) as.numeric(old$sign_concord[1]) else NA_real_
  )

  val_rows[[fac]] <- data.table(
    factor = fac, analysis_set = aset,
    new_n_ukb_zge8 = n_ukb, old_n_ukb_zge8 = if (nrow(old) == 1L) as.integer(old$n_leads[1]) else NA_integer_,
    n_match = val_n,
    new_sign = sign_u, old_sign = if (nrow(old) == 1L) as.numeric(old$sign_concord[1]) else NA_real_,
    sign_match = val_s
  )

  ukb_rows[[fac]] <- uk_leads
  karj5_rows[[fac]] <- kj_leads
  karj_p5_rows[[fac]] <- kj_p5
  pair_rows[[fac]] <- pairs

  logmsg(sprintf(
    "%s [%s]: UKB |Z|>=8 n=%d (sign=%.3f, mhc=%d, overlap=%d, only=%d); Karj |Z|>=5 n=%d (p5e8=%d, overlap=%d, only=%d); closest same-sign=%d/%d; val_n=%s val_s=%s",
    fac, aset, n_ukb, ifelse(is.na(sign_u), NA, sign_u),
    if (n_ukb > 0) sum(uk_leads$in_mhc) else 0L, n_ukb_ov, n_ukb - n_ukb_ov,
    n_karj, n_karj_p5, n_karj_ov, n_karj - n_karj_ov,
    n_closest_same, n_closest,
    as.character(val_n), as.character(val_s)
  ))
}

bind_or_empty <- function(lst, proto) {
  keep <- lst[vapply(lst, function(x) is.data.table(x) && nrow(x) > 0, logical(1))]
  if (length(keep) == 0) return(proto[0])
  rbindlist(keep, fill = TRUE)
}

proto_ukb <- data.table(
  factor = character(), analysis_set = character(), SNP = character(),
  CHR = integer(), BP = numeric(), A1 = character(), A2 = character(),
  Z_UKB = numeric(), Z_KARJ = numeric(), sign_concord = integer(),
  p_UKB = numeric(), p_KARJ = numeric(), chi2_UKB = numeric(),
  in_mhc = logical(), clump_kb = integer(), z_threshold = numeric()
)
proto_karj <- data.table(
  factor = character(), analysis_set = character(), SNP = character(),
  CHR = integer(), BP = numeric(), A1 = character(), A2 = character(),
  Z_KARJ = numeric(), Z_UKB = numeric(), sign_concord = integer(),
  p_KARJ = numeric(), p_UKB = numeric(), chi2_KARJ = numeric(),
  in_mhc = logical(), clump_kb = integer(), z_threshold = numeric(),
  pass_p5e8 = logical()
)
proto_karj_p5 <- copy(proto_karj)
proto_karj_p5[, `:=`(z_threshold = NULL, pass_p5e8 = NULL, p_threshold = numeric())]
proto_pairs <- data.table(
  factor = character(), analysis_set = character(),
  ukb_id = integer(), karj_id = integer(),
  ukb_SNP = character(), ukb_chr = integer(), ukb_bp = numeric(), ukb_z = numeric(),
  karj_SNP = character(), karj_chr = integer(), karj_bp = numeric(), karj_z = numeric(),
  distance_bp = numeric(), ukb_in_mhc = logical(), karj_in_mhc = logical(),
  same_sign = integer(), is_closest_for_ukb = logical(), is_closest_for_karj = logical()
)

ukb_all <- bind_or_empty(ukb_rows, proto_ukb)
karj5_all <- bind_or_empty(karj5_rows, proto_karj)
karj_p5_all <- bind_or_empty(karj_p5_rows, proto_karj_p5)
pairs_all <- bind_or_empty(pair_rows, proto_pairs)
sum_all <- rbindlist(sum_rows)
val_all <- rbindlist(val_rows)

write_tab <- function(x, path) {
  assert_write_path(path)
  fwrite(x, path, sep = "\t")
  logmsg("wrote ", path, " rows=", nrow(x))
}

write_tab(ukb_all[analysis_set == "main"], file.path(OUT, "ukb_leads_zge8_main.tsv"))
write_tab(ukb_all[analysis_set != "main"], file.path(OUT, "ukb_leads_zge8_supplement.tsv"))
write_tab(karj5_all[analysis_set == "main"], file.path(OUT, "karj_leads_zge5_main.tsv"))
write_tab(karj5_all[analysis_set != "main"], file.path(OUT, "karj_leads_zge5_supplement.tsv"))
write_tab(karj_p5_all[analysis_set == "main"], file.path(OUT, "karj_leads_p5e8_main.tsv"))
write_tab(karj_p5_all[analysis_set != "main"], file.path(OUT, "karj_leads_p5e8_supplement.tsv"))
write_tab(pairs_all[analysis_set == "main"], file.path(OUT, "locus_overlap_pairs_main.tsv"))
write_tab(pairs_all[analysis_set != "main"], file.path(OUT, "locus_overlap_pairs_supplement.tsv"))
write_tab(sum_all[analysis_set == "main"], file.path(OUT, "factor_locus_summary_main.tsv"))
write_tab(sum_all[analysis_set != "main"], file.path(OUT, "factor_locus_summary_supplement.tsv"))
write_tab(val_all, file.path(OUT, "lead_count_validation_vs_gls.tsv"))

# LAVA quantitative-trait inputs for the 13 main factors only (not run here if UKB v1.1 is incomplete).
write_status("LAVA_PREP", "Writing LAVA sumstats for 13 main factors; not launching LAVA if the UKB v1.1 LD reference is incomplete.")
logmsg("Preparing LAVA sumstats for 13 main factors")
lava_info <- list()
for (fac in MAIN_FACTORS) {
  zu <- dt[[paste0("ZU_", fac)]]
  zk <- dt[[paste0("ZK_", fac)]]
  ok_u <- is.finite(zu) & !is.na(dt$SNP) & nchar(dt$SNP) > 0 &
    toupper(dt$A1) %in% c("A","C","G","T") & toupper(dt$A2) %in% c("A","C","G","T") &
    toupper(dt$A1) != toupper(dt$A2)
  ok_k <- is.finite(zk) & !is.na(dt$SNP) & nchar(dt$SNP) > 0 &
    toupper(dt$A1) %in% c("A","C","G","T") & toupper(dt$A2) %in% c("A","C","G","T") &
    toupper(dt$A1) != toupper(dt$A2)
  uk_lava <- dt[ok_u, .(SNP = tolower(SNP), A1 = toupper(A1), A2 = toupper(A2), Z = zu[ok_u], N = N_UKB)]
  kj_lava <- dt[ok_k, .(SNP = tolower(SNP), A1 = toupper(A1), A2 = toupper(A2), Z = zk[ok_k], N = N_KARJ)]
  uk_lava <- unique(uk_lava, by = "SNP")
  kj_lava <- unique(kj_lava, by = "SNP")
  uk_path <- file.path(LAVA_OUT, paste0("UKB_", fac, ".lava.tsv.gz"))
  kj_path <- file.path(LAVA_OUT, paste0("KARJ_", fac, ".lava.tsv.gz"))
  assert_write_path(uk_path); assert_write_path(kj_path)
  fwrite(uk_lava, uk_path, sep = "\t")
  fwrite(kj_lava, kj_path, sep = "\t")
  lava_info[[length(lava_info) + 1L]] <- data.table(
    phenotype = paste0("UKB_", fac), cases = NA_real_, controls = NA_real_,
    filename = uk_path, n_snps = nrow(uk_lava), N = N_UKB, cohort = "UKB", factor = fac
  )
  lava_info[[length(lava_info) + 1L]] <- data.table(
    phenotype = paste0("KARJ_", fac), cases = NA_real_, controls = NA_real_,
    filename = kj_path, n_snps = nrow(kj_lava), N = N_KARJ, cohort = "KARJ", factor = fac
  )
  logmsg("LAVA sumstats ", fac, " UKB=", nrow(uk_lava), " KARJ=", nrow(kj_lava))
  rm(uk_lava, kj_lava)
}
lava_info_dt <- rbindlist(lava_info)
write_tab(lava_info_dt, file.path(OUT, "lava_input_info_full.tsv"))
# Official LAVA info has only phenotype, cases, controls, filename
write_tab(lava_info_dt[, .(phenotype, cases, controls, filename)], file.path(OUT, "lava_input.info.txt"))

# Check UKB v1.1 LD reference completeness (chr1-22).
ref_dir <- Sys.getenv("LAVA_UKB_REFERENCE_DIR", unset = file.path(PROJ, "external/lava/reference_ukb_v1.1/selected"))
wanted <- 1:22
have_bcor <- integer(0)
if (dir.exists(ref_dir)) {
  bcor <- list.files(ref_dir, pattern = "^lava-ukb-v1\\.1_chr([0-9]+)\\.bcor$")
  have_bcor <- as.integer(sub("^lava-ukb-v1\\.1_chr([0-9]+)\\.bcor$", "\\1", bcor))
  have_bcor <- sort(unique(have_bcor[have_bcor %in% 1:22]))
}
missing_chr <- setdiff(wanted, have_bcor)
lava_ref_complete <- length(missing_chr) == 0L
lava_ref_status <- data.table(
  reference = "UKB_EUR_v1.1",
  ref_dir = ref_dir,
  n_chr_present = length(have_bcor),
  chr_present = paste(have_bcor, collapse = ","),
  chr_missing = if (length(missing_chr) == 0) "" else paste(missing_chr, collapse = ","),
  complete_chr1_22 = lava_ref_complete,
  blocks_file = Sys.getenv("LAVA_LOCUS_FILE", unset = file.path(PROJ, "external/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")),
  lava_library = Sys.getenv("LAVA_R_LIBRARY", unset = file.path(PROJ, "external/lava/R_library_0.1.5")),
  note = if (lava_ref_complete) {
    "Reference complete; LAVA can be launched with run_exp52_9_lava.R"
  } else {
    "chr3-9 bcor/info are missing from the official UKB v1.1 LD reference. LAVA was not launched. Do not substitute 1000G as the main LD reference without a pre-specified sensitivity analysis (prior project audit: 1000G can bias local h2)."
  }
)
write_tab(lava_ref_status, file.path(OUT, "lava_reference_status.tsv"))
logmsg("LAVA ref complete=", lava_ref_complete, " missing=", paste(missing_chr, collapse = ","))

# ----- Report -----
write_status("REPORT", "Writing LOCAL_GENETIC_REPORT.md")
main_sum <- sum_all[analysis_set == "main"]
supp_sum <- sum_all[analysis_set != "main"]
n_val_fail <- sum(!val_all$n_match, na.rm = TRUE) + sum(!val_all$sign_match, na.rm = TRUE)

md_table <- function(d, cols) {
  x <- d[, ..cols]
  hdr <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  body <- apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

fmt_main <- copy(main_sum)
fmt_main[, sign_concord_ukb_leads := fmt_pct(sign_concord_ukb_leads)]
fmt_main[, sign_concord_closest_pairs := fmt_pct(sign_concord_closest_pairs)]
fmt_main[, max_chi2_share_ukb := fmt_pct(max_chi2_share_ukb)]
fmt_main[, top3_chi2_share_ukb := fmt_pct(top3_chi2_share_ukb)]

lines <- c(
  "# 52+9 主保留因子的位点水平跨队列复制",
  "",
  paste0("完成时间：", now_str(), "（", TZ_LABEL, "）"),
  "",
  "## 范围与预指定规则",
  "",
  "- 主面板为 52+9（61 性状）。本次只分析已锁定的 13 个主保留因子：F9、F11、F13、F14、F15、F17、F18、F20、F22、F23、F24、F25、F26。",
  "- F1、F4、F12、F19、F28 为拦截不足但有复制证据的因子，仅作补充表，不进入主结论。Watch 集 F5、F7、F21、F27 本次未制主表。",
  "- UKB lead 目录：与 Step 2 / QC 相同的 1 Mb 距离 clump，阈值 |Z_UKB|≥8。",
  "- Karjalainen 侧独立做 1 Mb 距离 clump，主重叠规则为 |Z_KARJ|≥5；同时报告 p<5×10⁻⁸ 的敏感计数。",
  "- 重叠定义：同一染色体上 UKB lead 与 Karj lead 的距离 ≤1 Mb。一个 UKB lead 若对应多个 Karj lead，符号一致性以距离最近的一对计。",
  "- MHC（chr6:25–34 Mb，GRCh37）在 lead 目录中标记、不删除。旧的 gls_Rreg lead 汇总同样未从 clump 中排除 MHC（MHC 只在 LDSC 中被去掉）。",
  "- 是否被少数大效应位点主导，只用描述性指标：独立 lead 数，以及最强 lead 的 χ² 占全部 UKB lead χ² 之和的比例（另附前三个 lead 的合计比例）。这不是遗传力分解。",
  "- 阈值在看结果前已预指定，事后未改。",
  "- 局部位点重叠失败不能解读为全局因子不复制。",
  "",
  "## 已有计算与本次新跑",
  "",
  "已有、只读复用：",
  "",
  "- `results/exp52_9/gls_Rreg/fgwas_ukbb_eur_gls_Rreg_genomewide.tsv.gz` 与 Karjalainen 对应文件（因子 GWAS，未重跑 GLS）。",
  "- `cross_cohort_leads_gls_Rreg.tsv`：仅有 UKB 侧 1 Mb clump 的**汇总计数**（全部 SNP / p<5e-8 / |Z|≥8 三套），含符号一致率，但没有 SNP 水平目录，也没有 Karjalainen 独立 clump，也没有双向区间重叠。",
  "- `snp_map.qs`（SNP、CHR、BP、A1、A2）。",
  "",
  "本次新跑：",
  "",
  "- 13 个主因子（外加 5 个补充因子，因读取基因组文件的开销已经发生，未再单独开昂贵作业）的 UKB |Z|≥8 lead 目录，含 chr、pos、rsID、两侧 Z、符号、Karj p、MHC 标记。",
  "- Karjalainen 独立 1 Mb clump（主规则 |Z|≥5，并输出 p<5e-8 目录与计数）。",
  "- 双向 1 Mb 重叠、最近配对的符号一致性，以及按因子汇总。",
  "- 与旧汇总表核对 UKB |Z|≥8 的 lead 数和符号一致率。",
  "- 为 13 个主因子准备了 LAVA 输入（SNP、A1、A2、Z、N），但因官方 UKB v1.1 LD 参考缺少 chr3–9，本次**没有启动** LAVA。",
  "",
  "## 与旧 lead 汇总的核对",
  ""
)

if (n_val_fail == 0 && all(val_all[factor %in% MAIN_FACTORS]$n_match) && all(val_all[factor %in% MAIN_FACTORS]$sign_match)) {
  lines <- c(lines, "13 个主因子以及 5 个补充因子的 UKB |Z|≥8 lead 数、符号一致率，均与 `cross_cohort_leads_gls_Rreg.tsv` 中 `ukb_clumped_1mb_zge8` 一致。说明本次目录使用了同一套距离 clump 规则和同一批因子 GWAS。")
} else {
  lines <- c(lines, "**核对未全部通过。** 详见 `lead_count_validation_vs_gls.tsv`，主表数字以本次目录为准并应先排查差异。")
}

lines <- c(
  lines,
  "",
  md_table(val_all[factor %in% MAIN_FACTORS], c("factor","new_n_ukb_zge8","old_n_ukb_zge8","n_match","new_sign","old_sign","sign_match")),
  "",
  "## 主分析：13 个保留因子",
  "",
  "下表为预指定主结果。`n_ukb_overlap` / `n_ukb_only` 以 UKB |Z|≥8 lead 为起点；`n_karj_overlap` / `n_karj_only` 以 Karj |Z|≥5 lead 为起点。`sign_concord_ukb_leads` 是同一 SNP 上 UKB lead 与其 Karjalainen Z 的符号一致率（与旧 QC 定义相同）。`sign_concord_closest_pairs` 是重叠区间内最近 UKB–Karj lead 对的符号一致率。",
  "",
  md_table(
    fmt_main,
    c("factor","n_loci_ukb_zge8","sign_concord_ukb_leads","n_ukb_overlap","n_ukb_only",
      "n_loci_karj_zge5","n_loci_karj_p5e8","n_karj_overlap","n_karj_only",
      "sign_concord_closest_pairs","n_ukb_mhc","max_chi2_share_ukb","top3_chi2_share_ukb")
  ),
  ""
)

# Short per-factor notes from the numbers actually computed.
for (i in seq_len(nrow(main_sum))) {
  r <- main_sum[i]
  dominance <- if (!is.na(r$max_chi2_share_ukb) && r$max_chi2_share_ukb >= 0.5) {
    "从 lead χ² 份额看，较接近由少数大效应位点主导"
  } else if (!is.na(r$max_chi2_share_ukb) && r$max_chi2_share_ukb >= 0.25) {
    "最强 lead 占有可观的 χ² 份额，但并非单一位点完全主导"
  } else {
    "lead χ² 较分散，不太像由极少数大效应位点单独驱动"
  }
  lines <- c(lines, paste0(
    "- **", r$factor, "**：UKB |Z|≥8 lead ", r$n_loci_ukb_zge8,
    " 个（MHC ", r$n_ukb_mhc, "；同 SNP 符号一致 ", r$n_same_sign_ukb_leads, "/", r$n_loci_ukb_zge8,
    "，", fmt_pct(r$sign_concord_ukb_leads), "）。",
    "与 Karj |Z|≥5 lead 重叠 ", r$n_ukb_overlap, " 个，仅 UKB ", r$n_ukb_only,
    " 个。Karj |Z|≥5 lead ", r$n_loci_karj_zge5, " 个（其中 p<5e-8：", r$n_loci_karj_p5e8,
    "），重叠 ", r$n_karj_overlap, "、仅 Karj ", r$n_karj_only, "。",
    "最近配对符号一致 ", r$n_closest_same_sign, "/", r$n_closest_pairs,
    if (is.na(r$sign_concord_closest_pairs)) "" else paste0("（", fmt_pct(r$sign_concord_closest_pairs), "）"),
    "。最强 UKB lead 的 χ² 占比 ", fmt_pct(r$max_chi2_share_ukb),
    "，前三个合计 ", fmt_pct(r$top3_chi2_share_ukb), "。", dominance, "。"
  ))
}

if (nrow(supp_sum) > 0) {
  fmt_s <- copy(supp_sum)
  fmt_s[, sign_concord_ukb_leads := fmt_pct(sign_concord_ukb_leads)]
  fmt_s[, sign_concord_closest_pairs := fmt_pct(sign_concord_closest_pairs)]
  fmt_s[, max_chi2_share_ukb := fmt_pct(max_chi2_share_ukb)]
  lines <- c(
    lines,
    "",
    "## 补充表：拦截不足但有复制证据的五因子（非主分析）",
    "",
    "这些因子未通过预指定的 intercept 门槛，不进入主分析，也未单独重跑昂贵计算；只是在同一批已读入的因子 GWAS 上顺带抽出 lead 目录与重叠，避免以后还要再读一遍约 1 GB 的压缩文件。",
    "",
    md_table(
      fmt_s,
      c("factor","n_loci_ukb_zge8","sign_concord_ukb_leads","n_ukb_overlap","n_ukb_only",
        "n_loci_karj_zge5","n_loci_karj_p5e8","n_karj_overlap","n_karj_only",
        "sign_concord_closest_pairs","max_chi2_share_ukb")
    )
  )
}

lines <- c(
  lines,
  "",
  "## LAVA 状态",
  "",
  paste0("- LAVA 0.1.5 库由环境变量 `LAVA_R_LIBRARY` 指定。"),
  paste0("- 标准区块文件由环境变量 `LAVA_LOCUS_FILE` 指定。"),
  paste0("- 主 LD 参考由 `LAVA_UKB_REFERENCE_DIR` 指向官方 UK Biobank 欧洲 v1.1 二进制 LD，而不是 1000 Genomes。先前项目审计认为 1000G 可能偏倚局部遗传力。"),
  paste0("- 当前已有染色体：", if (length(have_bcor) == 0) "无" else paste(have_bcor, collapse = ", "), "。"),
  paste0("- 缺少染色体：", if (length(missing_chr) == 0) "无" else paste(missing_chr, collapse = ", "), "。"),
  if (lava_ref_complete) {
    "- 参考已齐，可用 `run_exp52_9_lava.R` 启动。本次脚本在写出输入后未自动开始 LAVA，以便先确认 lead 分析完整。"
  } else {
    "- **因缺少 chr3–9 的 .bcor/.info，本次没有启动 LAVA。** 不要为了补参考而去重启已经存在的 curl/aria2/Catalog 下载。补齐 chr3–9 后，用已写好的 `lava_sumstats/` 与 `run_exp52_9_lava.R` 即可。预指定流程是：单变量局部 h2 p<0.05 过滤后，再做同一因子 UKB 对 Karjalainen 的双变量局部 rg；不做因子之间、也不做因子–疾病 LAVA。样本重叠按 0 处理。"
  },
  paste0("- 13 个主因子的 LAVA 输入已写出，N_UKB=", N_UKB, "，N_KARJ=", N_KARJ, "。"),
  "",
  "## 输出文件",
  "",
  "- UKB lead 目录（主）：`ukb_leads_zge8_main.tsv`",
  "- Karj |Z|≥5 目录（主）：`karj_leads_zge5_main.tsv`",
  "- Karj p<5e-8 目录（主，敏感）：`karj_leads_p5e8_main.tsv`",
  "- 重叠对（主）：`locus_overlap_pairs_main.tsv`",
  "- 按因子汇总（主）：`factor_locus_summary_main.tsv`",
  "- 补充表：对应的 `*_supplement.tsv`",
  "- 与旧汇总核对：`lead_count_validation_vs_gls.tsv`",
  "- LAVA 输入：`lava_sumstats/`，`lava_input.info.txt`，`lava_reference_status.tsv`",
  "- 本报告：`LOCAL_GENETIC_REPORT.md`",
  "",
  "## 脚本",
  "",
  "- `run_exp52_9_local_genetic.R`",
  "- `run_exp52_9_lava.R`（LAVA 启动脚本；参考不齐时不要跑）",
  ""
)

assert_write_path(REPORT)
writeLines(lines, REPORT, useBytes = TRUE)
logmsg("wrote report ", REPORT)

assert_write_path(DONE_FLAG)
writeLines(c(paste0("done=", now_str()), paste0("pid=", Sys.getpid()), paste0("lava_ref_complete=", lava_ref_complete)), DONE_FLAG)
write_status("DONE", paste0("Lead catalogs and summaries written. LAVA launched: no. Reference complete: ", lava_ref_complete, "."))
logmsg("ALL DONE")
