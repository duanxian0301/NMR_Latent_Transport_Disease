#!/usr/bin/env Rscript
# Independent 52+9 (61-trait) genome-wide GLS projection (gls_Rreg).
# Freeze: GFA seed 101, k=30. Project F genome-wide onto UKBB_EUR and
# Karjalainen with official-style GLS (Sigma = cohort residual R,
# cond_num=1e3, plus_resid=FALSE). Then the same calibration and
# one-to-one cross-cohort metrics as matched-52 Step 2.
#
# Outputs ONLY under results/exp52_9/ (GLS/QC under gls_Rreg/; new
# Karjalainen Z cache and Karjalainen 61x61 R in cache/ and R/).
# Does NOT overwrite frozen matched-52 / main GFA / week1 cache /
# gls_Rreg / step2_ldclump objects. Does NOT refit GFA, lock kmax,
# start GLEANR GRID, disease F-GWAS/MR, or make a retain/drop
# factor-QC decision.
#
# Pipeline (resumable / idempotent):
#   1. Build 61-trait Karjalainen Z on the 52+9 common SNP panel
#      (reuse frozen 52-trait Karj Z read-only; add 9 expansion
#      traits from HM3 with the same allele alignment as the UKB Z).
#   2. Estimate a new 61x61 official-style R_ldsc for Karjalainen
#      and condition at 1e3 in results/exp52_9/R/.
#   3. Genome-wide GLS, both cohorts, frozen seed-101 fit RDS
#      (F_hat * scale; not F_hat alone).
#   4. lambda_GC, LDSC intercept/h2, Pearson/sign/leads/LDSC rg.

suppressPackageStartupMessages({
  library(GFA)
  library(data.table)
  library(qs)
  library(Matrix)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
EXP  <- file.path(PROJ, "results/exp52_9")
OUT  <- file.path(EXP, "gls_Rreg")
CACHE <- file.path(EXP, "cache")
RDIR <- file.path(EXP, "R")
GFA_DIR <- file.path(EXP, "gfa")

# Frozen matched-52 / week1 objects: READ-ONLY, never write.
FROZEN_Z_KJ <- file.path(PROJ, "results/week1/fgwas/cache/Z_karjalainen.qs")
FROZEN_Z_UK <- file.path(PROJ, "results/week1/fgwas/cache/Z_ukbb_eur.qs")
FROZEN_MAP <- file.path(PROJ, "results/week1/fgwas/cache/snp_map.qs")
FROZEN_R_UK <- file.path(PROJ, "results/week1/fgwas/gls_Rreg/R_ldsc_cond1e3_ukbb_eur.tsv")
FROZEN_R_KJ <- file.path(PROJ, "results/week1/fgwas/gls_Rreg/R_ldsc_cond1e3_karjalainen.tsv")
FROZEN_STEP2 <- file.path(PROJ, "results/week1/fgwas/step2_ldclump")

COMMON <- file.path(EXP, "common_snps.txt")
TRAIT_N <- file.path(EXP, "trait_N.tsv")
CROSS61 <- file.path(PROJ, "inputs/crosswalk/exp52_9_trait_gcst_list.tsv")
KARJ9 <- file.path(PROJ, "inputs/sumstats/exp52_9/karjalainen")
KARJ52 <- file.path(PROJ, "inputs/sumstats/karjalainen")

Z_UK_PATH <- file.path(CACHE, "Z_ukbb_eur.qs")
Z_KJ_PATH <- file.path(CACHE, "Z_karjalainen.qs")
MAP_PATH <- file.path(CACHE, "snp_map.qs")
ZSTATS_KJ <- file.path(CACHE, "Z_build_stats_karjalainen.tsv")

R_UK_PATH <- file.path(RDIR, "R_ldsc_cond1e3_ukbb_eur.tsv")
R_KJ_PATH <- file.path(RDIR, "R_ldsc_cond1e3_karjalainen.tsv")
R_KJ_PAIR <- file.path(RDIR, "R_ldsc_pairs_karjalainen.tsv")
R_KJ_RAW <- file.path(RDIR, "R_ldsc_raw_karjalainen.tsv")
R_META <- file.path(RDIR, "R_meta.tsv")

FIT_RDS <- file.path(GFA_DIR, "ukbb_ldclump_seed101.rds")
F_HAT_TSV <- file.path(GFA_DIR, "ukbb_ldclump_seed101_F_hat.tsv")

LD_DIR <- Sys.getenv("LDSC_EUR_LD_DIR", unset = file.path(PROJ, "external/ldsc/eur_w_ld_chr"))

LOW_H2Z <- c("Omega_3", "DHA", "XL_HDL_TG", "XS_VLDL_L", "XS_VLDL_FC")
N_UKB <- 413897
N_KARJ <- 136016
N_KARJ_FALLBACK <- 136016
MHC_CHR <- 6L
MHC_START <- 25e6
MHC_END <- 34e6
CLUMP_KB_MAIN <- 1000
CHISQ_TWO_STEP <- 30
COND_NUM <- 1000
METHOD <- "gls_Rreg"
SEED_PRIMARY <- 101
EXPECTED_K <- 30L
TZ_LABEL <- "Asia/Shanghai"

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "partial"), recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

UK_FG <- file.path(OUT, "fgwas_ukbb_eur_gls_Rreg_genomewide.tsv.gz")
KJ_FG <- file.path(OUT, "fgwas_karjalainen_gls_Rreg_genomewide.tsv.gz")
QC_TAB <- file.path(OUT, "qc_lambda_ldsc_gls_Rreg.tsv")
SUM_OUT <- file.path(OUT, "cross_cohort_gls_Rreg.tsv")
LEAD_OUT <- file.path(OUT, "cross_cohort_leads_gls_Rreg.tsv")
STATUS <- file.path(OUT, "status.md")
REPORT <- file.path(OUT, "EXP52_9_GLS_RREG_REPORT.md")
DONE_FLAG <- file.path(OUT, "GLS_RREG_DONE.txt")
PIDFILE <- file.path(OUT, "exp52_9_gls_Rreg.pid")
CK_Z_KJ <- file.path(OUT, "checkpoint_Z_karjalainen.txt")
CK_R_KJ <- file.path(OUT, "checkpoint_R_karjalainen.txt")
CK_GLS_UK <- file.path(OUT, "checkpoint_gls_ukbb_eur.txt")
CK_GLS_KJ <- file.path(OUT, "checkpoint_gls_karjalainen.txt")
CK_QC <- file.path(OUT, "checkpoint_qc_cross.txt")

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
  allowed <- normalizePath(EXP, winslash = "/", mustWork = FALSE)
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
  path_n <- normalizePath(path, winslash = "/", mustWork = FALSE)
  ok <- startsWith(tolower(parent), tolower(allowed)) ||
    startsWith(tolower(path_n), tolower(allowed))
  if (!isTRUE(ok)) {
    stop("Refusing to write outside results/exp52_9/: ", path)
  }
  banned <- c(
    file.path(PROJ, "results/week1"),
    file.path(PROJ, "inputs/sumstats")
  )
  for (b in banned) {
    bn <- normalizePath(b, winslash = "/", mustWork = FALSE)
    if (startsWith(tolower(parent), tolower(bn)) ||
        startsWith(tolower(path_n), tolower(bn))) {
      stop("Refusing to write under banned path ", b, " : ", path)
    }
  }
  invisible(TRUE)
}

write_status <- function(step, detail = "") {
  assert_write_path(STATUS)
  lines <- c(
    "# 52+9 genome-wide gls_Rreg projection — live status",
    "",
    paste0("Updated: ", now_str(), " (", TZ_LABEL, ")"),
    "",
    paste0("**Step:** ", step),
    "",
    detail,
    "",
    "## Policy",
    "- Independent 52+9 analysis. GLS/QC written only under `results/exp52_9/gls_Rreg/`.",
    "- Karjalainen Z cache: `results/exp52_9/cache/Z_karjalainen.qs` (new).",
    "- Karjalainen 61x61 R: `results/exp52_9/R/R_ldsc_cond1e3_karjalainen.tsv` (new).",
    "- Frozen matched-52 / week1 cache / gls_Rreg / step2_ldclump / main GFA are not overwritten.",
    "- Frozen GFA seed 101 (k=30) is read-only; GFA is not refit; kmax is not locked.",
    "- Official-style gls_Rreg: Sigma = cohort residual R, cond_num=1e3, plus_resid=FALSE; F from fit RDS as F_hat * scale.",
    "- All 61 traits kept, including five UKB univariate LDSC h2 Z that did not pass Z>4 (flagged, not dropped).",
    "- Glycine is not included.",
    "- No retain/drop factor-QC decision, no GLEANR GRID, no disease F-GWAS / MR.",
    "",
    "- Script: `results/exp52_9/run_exp52_9_gls_Rreg.R`",
    "- Log: `results/exp52_9/gls_Rreg/exp52_9_gls_Rreg_rscript.log`"
  )
  con <- file(STATUS, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}

lambda_gc <- function(z) {
  z <- z[is.finite(z)]
  median(z^2, na.rm = TRUE) / qchisq(0.5, 1)
}

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(as.numeric(x), digits = d, format = "f"))

mat_cond <- function(M) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  max(ev) / max(min(ev), .Machine$double.eps)
}

write_R_matrix <- function(M, path) {
  assert_write_path(path)
  dt <- as.data.table(M, keep.rownames = "trait")
  fwrite(dt, path, sep = "\t")
}

read_R_matrix <- function(path) {
  dt <- fread(path)
  rn <- dt$trait
  M <- as.matrix(dt[, -1, with = FALSE])
  rownames(M) <- rn
  storage.mode(M) <- "double"
  M
}

condition_pd <- function(S, cond_num, corr = FALSE) {
  ev0 <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  raw_cond <- max(ev0) / max(min(ev0), .Machine$double.eps)
  pd <- Matrix::nearPD(S, corr = corr, keepDiag = !corr, posd.tol = 1 / cond_num)
  M <- as.matrix(pd$mat)
  colnames(M) <- colnames(S)
  rownames(M) <- rownames(S)
  ev1 <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(
    mat = M,
    raw_min_ev = min(ev0), raw_max_ev = max(ev0), raw_cond = raw_cond,
    pd_min_ev = min(ev1), pd_max_ev = max(ev1),
    pd_cond = max(ev1) / max(min(ev1), .Machine$double.eps),
    diag_min = min(diag(M)), diag_max = max(diag(M)),
    offdiag_max_abs = max(abs(M[lower.tri(M)]))
  )
}

complement <- function(x) chartr("ACGTacgt", "TGCAtgca", x)
allele_relation <- function(a1, a2, b1, b2) {
  a1 <- toupper(a1); a2 <- toupper(a2); b1 <- toupper(b1); b2 <- toupper(b2)
  ok <- nchar(a1) == 1 & nchar(a2) == 1 & nchar(b1) == 1 & nchar(b2) == 1 &
    a1 %chin% c("A", "C", "G", "T") & a2 %chin% c("A", "C", "G", "T") &
    b1 %chin% c("A", "C", "G", "T") & b2 %chin% c("A", "C", "G", "T")
  rel <- rep("mismatch", length(a1))
  rel[!ok] <- "non_snv"
  rel[ok & a1 == b1 & a2 == b2] <- "match"
  rel[ok & a1 == b2 & a2 == b1] <- "flip"
  ca1 <- complement(a1); ca2 <- complement(a2)
  rel[ok & ca1 == b1 & ca2 == b2] <- "complement"
  rel[ok & ca1 == b2 & ca2 == b1] <- "complement_flip"
  rel
}

detect_kind <- function(path) {
  hdr <- names(fread(path, nrows = 0))
  if ("rsid" %in% hdr && "beta" %in% hdr && "standard_error" %in% hdr) return("catalog")
  if ("SNP" %in% hdr && "BETA" %in% hdr && "SE" %in% hdr) return("author")
  stop("Unknown schema: ", path, " hdr=", paste(hdr, collapse = ","))
}

read_hm3_z <- function(path) {
  kind <- detect_kind(path)
  if (kind == "catalog") {
    dt <- fread(path, select = c("rsid", "effect_allele", "other_allele", "beta", "standard_error"),
                colClasses = list(character = c("rsid", "effect_allele", "other_allele")))
    setnames(dt, c("rsid", "a1", "a2", "beta", "se"))
  } else {
    dt <- fread(path, select = c("SNP", "A1", "A2", "BETA", "SE"),
                colClasses = list(character = c("SNP", "A1", "A2")))
    setnames(dt, c("rsid", "a1", "a2", "beta", "se"))
  }
  dt[, `:=`(a1 = toupper(a1), a2 = toupper(a2))]
  dt <- dt[is.finite(beta) & is.finite(se) & se > 0]
  dt <- dt[!duplicated(rsid)]
  dt[, z := beta / se]
  dt[, `:=`(beta = NULL, se = NULL)]
  list(dt = dt, kind = kind)
}

read_ldscores <- function(ld_dir) {
  files <- list.files(ld_dir, pattern = "^[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
  if (length(files) == 0) stop("No LD score files in ", ld_dir)
  ld <- rbindlist(lapply(files, function(f) {
    fread(f, select = c("CHR", "SNP", "BP", "L2"))
  }), use.names = TRUE)
  ld <- ld[!duplicated(SNP)]
  mfiles <- list.files(ld_dir, pattern = "^[0-9]+\\.l2\\.M_5_50$", full.names = TRUE)
  M <- sum(vapply(mfiles, function(f) as.numeric(readLines(f, n = 1, warn = FALSE)), numeric(1)))
  list(ld = ld, M = M, n_ld_snps = nrow(ld))
}

estimate_R_ldsc <- function(Z, ldscores, ld_size, Nvec, side, out_dir) {
  pair_path <- file.path(out_dir, sprintf("R_ldsc_pairs_%s.tsv", side))
  raw_path <- file.path(out_dir, sprintf("R_ldsc_raw_%s.tsv", side))
  assert_write_path(pair_path)
  assert_write_path(raw_path)
  traits <- colnames(Z)
  M <- length(traits)
  stopifnot(identical(names(Nvec), traits))
  pairs <- data.table(trait1_idx = sequence(seq_len(M)),
                      trait2_idx = rep(seq_len(M), times = seq_len(M)))
  pairs[, `:=`(trait1 = traits[trait1_idx], trait2 = traits[trait2_idx])]
  done <- data.table(trait1 = character(), trait2 = character(), intercept = numeric())
  if (file_ok(pair_path)) {
    done <- tryCatch(fread(pair_path), error = function(e) done)
    if (!all(c("trait1", "trait2", "intercept") %in% names(done))) {
      done <- data.table(trait1 = character(), trait2 = character(), intercept = numeric())
    }
  }
  key_done <- paste(done$trait1, done$trait2, sep = "\t")
  todo <- pairs[!paste(trait1, trait2, sep = "\t") %chin% key_done]
  logmsg("R_ldsc ", side, ": ", nrow(done), " pairs already saved, ", nrow(todo),
         " remaining (total ", nrow(pairs), ")")
  if (nrow(todo) > 0) {
    J <- nrow(Z)
    Nmat <- matrix(rep(as.numeric(Nvec), each = J), nrow = J)
    colnames(Nmat) <- traits
    h2_cache <- vector("list", length(traits))
    names(h2_cache) <- traits
    logmsg("R_ldsc ", side, ": precomputing ", length(traits), " univariate snp_ldsc fits")
    t_univ <- proc.time()[3]
    for (tr in traits) {
      ixu <- which(is.finite(Z[, tr]) & is.finite(ldscores))
      h2_cache[[tr]] <- tryCatch(
        snp_ldsc(ld_score = ldscores[ixu], ld_size = ld_size,
                 chi2 = Z[ixu, tr]^2, sample_size = Nmat[ixu, tr],
                 blocks = NULL, chi2_thr1 = 30, chi2_thr2 = Inf),
        error = function(e) {
          logmsg("univar ldsc fail ", side, " ", tr, ": ", e$message)
          NULL
        }
      )
    }
    logmsg("R_ldsc ", side, ": univariates done in ", round(proc.time()[3] - t_univ), "s")
    t0 <- proc.time()[3]
    buf <- list()
    for (r in seq_len(nrow(todo))) {
      i <- todo$trait1[r]
      j <- todo$trait2[r]
      ix <- which(is.finite(Z[, i]) & is.finite(Z[, j]) & is.finite(ldscores))
      val <- tryCatch({
        rg <- ldsc_rg(
          ld_score = ldscores[ix], ld_size = ld_size,
          z1 = Z[ix, i], z2 = Z[ix, j],
          h2_1 = h2_cache[[i]], h2_2 = h2_cache[[j]],
          sample_size_1 = Nmat[ix, i], sample_size_2 = Nmat[ix, j],
          blocks = NULL, ncores = 1
        )
        as.numeric(rg[[1]])
      }, error = function(e) {
        logmsg("R_ldsc pair fail ", side, " ", i, "-", j, ": ", e$message)
        NA_real_
      })
      buf[[length(buf) + 1]] <- data.table(trait1 = i, trait2 = j, intercept = val)
      if (r %% 20 == 0 || r == nrow(todo)) {
        chunk <- rbindlist(buf)
        if (!file.exists(pair_path) || file.info(pair_path)$size < 10) {
          fwrite(chunk, pair_path, sep = "\t")
        } else {
          fwrite(chunk, pair_path, sep = "\t", append = TRUE)
        }
        buf <- list()
        elapsed <- proc.time()[3] - t0
        rate <- elapsed / r
        eta <- rate * (nrow(todo) - r)
        logmsg("R_ldsc ", side, " ", r, "/", nrow(todo),
               " elapsed=", round(elapsed), "s eta=", round(eta), "s")
      }
    }
    if (length(buf) > 0) {
      chunk <- rbindlist(buf)
      if (!file.exists(pair_path) || file.info(pair_path)$size < 10) {
        fwrite(chunk, pair_path, sep = "\t")
      } else {
        fwrite(chunk, pair_path, sep = "\t", append = TRUE)
      }
    }
    done <- fread(pair_path)
  }
  Se <- matrix(NA_real_, M, M, dimnames = list(traits, traits))
  for (r in seq_len(nrow(done))) {
    a <- done$trait1[r]; b <- done$trait2[r]; v <- done$intercept[r]
    if (!a %in% traits || !b %in% traits) next
    Se[a, b] <- v
    Se[b, a] <- v
  }
  if (any(!is.finite(Se))) {
    nbad <- sum(!is.finite(Se))
    stop("R_ldsc ", side, " has ", nbad, " non-finite entries after assembly")
  }
  write_R_matrix(Se, raw_path)
  Se
}

wls_coef <- function(y, x, wt, intercept = NULL) {
  if (is.null(intercept)) {
    as.numeric(lm.wfit(cbind(1, x), y, wt)$coefficients)
  } else {
    c(intercept, as.numeric(lm.wfit(cbind(x), y - intercept, wt)$coefficients[1]))
  }
}

ldsc_fit_once <- function(y, x, wld, intercept = NULL, niter = 3) {
  ic0 <- if (is.null(intercept)) max(mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE) * 0, 1) else intercept
  sl0 <- max((mean(y, na.rm = TRUE) - ic0) / max(mean(x, na.rm = TRUE), 1e-8), 0)
  pred <- ic0 + sl0 * x
  coef <- c(ic0, sl0)
  for (i in seq_len(niter)) {
    pred <- pmax(pred, 1e-6)
    wt <- 1 / (2 * pred^2 * pmax(wld, 1))
    wt[!is.finite(wt)] <- 0
    coef <- wls_coef(y, x, wt, intercept = intercept)
    pred <- coef[1] + coef[2] * x
  }
  list(intercept = coef[1], slope = coef[2])
}

block_indices <- function(n, n.blocks = 200) {
  if (n < n.blocks * 2) n.blocks <- max(2, floor(n / 50))
  cuts <- floor(seq(0, n, length.out = n.blocks + 1))
  blocks <- vector("list", n.blocks)
  for (b in seq_len(n.blocks)) blocks[[b]] <- seq.int(cuts[b] + 1, cuts[b + 1])
  blocks
}

chr_blocks <- function(chr) split(seq_along(chr), as.character(chr))

jackknife_se <- function(full_coef, block_coefs) {
  m <- length(block_coefs)
  arr <- do.call(rbind, block_coefs)
  bar <- colMeans(arr)
  dev <- sweep(arr, 2, bar, "-")
  sqrt(((m - 1) / m) * colSums(dev^2))
}

ldsc_univar <- function(snp, z, chr, bp, l2, wld, N, M, drop_mhc = TRUE,
                        two_step = FALSE, n.blocks = 200) {
  dt <- data.table(SNP = snp, z = z, CHR = chr, BP = bp, L2 = l2, W = wld)
  dt <- dt[is.finite(z) & is.finite(L2) & is.finite(W) & L2 > 0]
  if (drop_mhc) dt <- dt[!(CHR == MHC_CHR & BP >= MHC_START & BP <= MHC_END)]
  setorder(dt, CHR, BP)
  y <- dt$z^2
  x <- dt$L2
  w <- dt$W
  intercept_arg <- if (two_step) {
    sub <- y <= CHISQ_TWO_STEP
    if (sum(sub) < 1000) sub <- rep(TRUE, length(y))
    ldsc_fit_once(y[sub], x[sub], w[sub], intercept = NULL)$intercept
  } else {
    NULL
  }
  full <- ldsc_fit_once(y, x, w, intercept = intercept_arg)
  blocks <- chr_blocks(dt$CHR)
  if (length(blocks) < 5) blocks <- block_indices(nrow(dt), min(n.blocks, 22))
  bcoefs <- lapply(blocks, function(idx) {
    keep <- setdiff(seq_len(nrow(dt)), idx)
    if (length(keep) < 100) return(c(full$intercept, full$slope))
    r <- ldsc_fit_once(y[keep], x[keep], w[keep], intercept = intercept_arg)
    c(r$intercept, r$slope)
  })
  se <- jackknife_se(c(full$intercept, full$slope), bcoefs)
  h2 <- full$slope * M / N
  list(
    n_snp = nrow(dt),
    intercept = full$intercept,
    intercept_se = se[1],
    slope = full$slope,
    slope_se = se[2],
    h2 = h2,
    h2_se = se[2] * M / N,
    mean_chisq = mean(y),
    lambda_gc = lambda_gc(dt$z),
    M = M,
    N = N,
    two_step = two_step
  )
}

ldsc_bivar <- function(snp, z1, z2, chr, bp, l2, wld, N1, N2, h21, h22, M,
                       drop_mhc = TRUE, n.blocks = 200) {
  dt <- data.table(SNP = snp, z1 = z1, z2 = z2, CHR = chr, BP = bp, L2 = l2, W = wld)
  dt <- dt[is.finite(z1) & is.finite(z2) & is.finite(L2) & is.finite(W) & L2 > 0]
  if (drop_mhc) dt <- dt[!(CHR == MHC_CHR & BP >= MHC_START & BP <= MHC_END)]
  setorder(dt, CHR, BP)
  y <- dt$z1 * dt$z2
  x <- dt$L2
  w <- dt$W
  full <- ldsc_fit_once(y, x, w, intercept = NULL)
  pred <- full$intercept + full$slope * x
  wt <- 1 / (2 * pmax(abs(pred), 1e-3)^2 * pmax(w, 1))
  coef <- wls_coef(y, x, wt, intercept = NULL)
  full$intercept <- coef[1]
  full$slope <- coef[2]
  blocks <- chr_blocks(dt$CHR)
  if (length(blocks) < 5) blocks <- block_indices(nrow(dt), min(n.blocks, 22))
  bcoefs <- lapply(blocks, function(idx) {
    keep <- setdiff(seq_len(nrow(dt)), idx)
    if (length(keep) < 100) return(c(full$intercept, full$slope))
    r <- ldsc_fit_once(y[keep], x[keep], w[keep], intercept = NULL)
    pr <- r$intercept + r$slope * x[keep]
    wt2 <- 1 / (2 * pmax(abs(pr), 1e-3)^2 * pmax(w[keep], 1))
    cf <- wls_coef(y[keep], x[keep], wt2, intercept = NULL)
    cf
  })
  se <- jackknife_se(c(full$intercept, full$slope), bcoefs)
  gcov <- full$slope * M / sqrt(N1 * N2)
  gcov_se <- se[2] * M / sqrt(N1 * N2)
  rg <- gcov / sqrt(max(h21, 1e-8) * max(h22, 1e-8))
  rg_se <- gcov_se / sqrt(max(h21, 1e-8) * max(h22, 1e-8))
  list(
    n_snp = nrow(dt),
    gcov_intercept = full$intercept,
    gcov_intercept_se = se[1],
    gcov = gcov,
    gcov_se = gcov_se,
    rg = rg,
    rg_se = rg_se,
    mean_z1z2 = mean(y)
  )
}

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

run_gls <- function(Z, fit, R, plus_resid = FALSE) {
  if (any(!is.finite(Z))) {
    Z[!is.finite(Z)] <- 0
  }
  if (is.null(fit$F_hat) || is.null(fit$scale)) {
    stop("Fit object missing F_hat or scale; refusing to project from F_hat alone")
  }
  myF <- fit$F_hat * fit$scale
  if (plus_resid) {
    myR <- R + diag(as.numeric(fit$fit$residuals_sd)^2)
  } else {
    myR <- R
  }
  ev <- eigen(myR, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) stop("Sigma is not PD; min ev=", min(ev))
  S1 <- matrix(1, nrow = nrow(Z), ncol = ncol(Z))
  ret <- GFA:::loadings_gls(X = Z, S = S1, R = myR, F_hat = myF)
  ret$P <- 2 * pnorm(-abs(ret$L / ret$S))
  ret$Sigma_min_ev <- min(ev)
  ret$Sigma_cond <- max(ev) / min(ev)
  ret
}

write_fgwas <- function(ret, snps, path) {
  assert_write_path(path)
  L <- as.matrix(ret$L)
  S <- as.matrix(ret$S)
  P <- as.matrix(ret$P)
  k <- ncol(L)
  dt <- data.table(SNP = snps)
  for (j in seq_len(k)) {
    dt[, paste0("BETA_F", j) := L[, j]]
    dt[, paste0("SE_F", j) := S[, j]]
    dt[, paste0("P_F", j) := P[, j]]
    dt[, paste0("Z_F", j) := L[, j] / S[, j]]
  }
  fwrite(dt, path, sep = "\t")
  invisible(path)
}

cross_metrics <- function(snps, z_u, z_k, map, ld, M, h2u, h2k, factor, method) {
  dt <- data.table(SNP = snps, zu = z_u, zk = z_k)
  dt <- merge(dt, map[, .(SNP, CHR, BP)], by = "SNP", all.x = TRUE)
  dt <- merge(dt, ld[, .(SNP, L2)], by = "SNP", all.x = TRUE)
  ok <- is.finite(dt$zu) & is.finite(dt$zk)
  pear <- cor(dt$zu[ok], dt$zk[ok], method = "pearson")
  spear <- cor(dt$zu[ok], dt$zk[ok], method = "spearman")
  sign_all <- mean(sign(dt$zu[ok]) == sign(dt$zk[ok]))
  thresh <- c(2, 3, 5, 8)
  sign_thr <- sapply(thresh, function(t) {
    w <- ok & abs(dt$zu) >= t
    if (sum(w) < 10) return(NA_real_)
    mean(sign(dt$zu[w]) == sign(dt$zk[w]))
  })
  cor_thr <- sapply(thresh, function(t) {
    w <- ok & abs(dt$zu) >= t
    if (sum(w) < 10) return(NA_real_)
    cor(dt$zu[w], dt$zk[w])
  })
  cl <- clump_distance(dt$CHR, dt$BP, dt$zu, kb = CLUMP_KB_MAIN)
  lead_all <- cl & ok
  lead_p5 <- cl & ok & (2 * pnorm(-abs(dt$zu)) < 5e-8)
  rep_row <- function(mask, label) {
    n <- sum(mask)
    if (n == 0) {
      return(data.table(lead_set = label, n_leads = 0, sign_concord = NA_real_,
                        n_same_sign = 0L, n_p05 = 0L, n_p1e3 = 0L, n_p5e8 = 0L,
                        median_z_ukb = NA_real_, median_z_karj = NA_real_))
    }
    data.table(
      lead_set = label, n_leads = n,
      sign_concord = mean(sign(dt$zu[mask]) == sign(dt$zk[mask])),
      n_same_sign = sum(sign(dt$zu[mask]) == sign(dt$zk[mask])),
      n_p05 = sum(2 * pnorm(-abs(dt$zk[mask])) < 0.05),
      n_p1e3 = sum(2 * pnorm(-abs(dt$zk[mask])) < 1e-3),
      n_p5e8 = sum(2 * pnorm(-abs(dt$zk[mask])) < 5e-8),
      median_z_ukb = median(abs(dt$zu[mask])),
      median_z_karj = median(abs(dt$zk[mask]))
    )
  }
  leads <- rbind(
    rep_row(lead_all, "ukb_clumped_1mb_all"),
    rep_row(lead_p5, "ukb_clumped_1mb_p5e8"),
    rep_row(cl & ok & abs(dt$zu) >= 8, "ukb_clumped_1mb_zge8")
  )
  biv <- tryCatch(
    ldsc_bivar(dt$SNP, dt$zu, dt$zk, dt$CHR, dt$BP, dt$L2, dt$L2,
               N_UKB, N_KARJ, h2u, h2k, M, drop_mhc = TRUE),
    error = function(e) { logmsg("bivar ldsc fail ", method, " ", factor, ": ", e$message); NULL }
  )
  list(
    summary = data.table(
      method = method, factor = factor,
      n_common = sum(ok),
      pearson_z = pear,
      spearman_z = spear,
      sign_concord_all = sign_all,
      sign_concord_zukb_ge2 = sign_thr[1],
      sign_concord_zukb_ge3 = sign_thr[2],
      sign_concord_zukb_ge5 = sign_thr[3],
      sign_concord_zukb_ge8 = sign_thr[4],
      pearson_zukb_ge2 = cor_thr[1],
      pearson_zukb_ge3 = cor_thr[2],
      pearson_zukb_ge5 = cor_thr[3],
      pearson_zukb_ge8 = cor_thr[4],
      ldsc_rg = if (is.null(biv)) NA else biv$rg,
      ldsc_rg_se = if (is.null(biv)) NA else biv$rg_se,
      ldsc_gcov_intercept = if (is.null(biv)) NA else biv$gcov_intercept,
      ldsc_gcov_intercept_se = if (is.null(biv)) NA else biv$gcov_intercept_se,
      n_ldsc_rg = if (is.null(biv)) NA else biv$n_snp
    ),
    leads = cbind(data.table(method = method, factor = factor), leads)
  )
}

valid_kj_z <- function(traits, common) {
  if (!(file_ok(CK_Z_KJ, 10) && file_ok(Z_KJ_PATH, 1e6))) return(FALSE)
  obj <- tryCatch(qs::qread(Z_KJ_PATH), error = function(e) NULL)
  if (is.null(obj) || is.null(obj$Z) || !is.matrix(obj$Z)) return(FALSE)
  Z <- obj$Z
  isTRUE(nrow(Z) == length(common) && ncol(Z) == 61 &&
           setequal(colnames(Z), traits) &&
           setequal(rownames(Z), common))
}

valid_kj_r <- function(traits) {
  if (!(file_ok(CK_R_KJ, 10) && file_ok(R_KJ_PATH, 200) && file_ok(R_KJ_RAW, 200))) return(FALSE)
  Rchk <- tryCatch(read_R_matrix(R_KJ_PATH), error = function(e) NULL)
  isTRUE(!is.null(Rchk) && nrow(Rchk) == 61 && ncol(Rchk) == 61 &&
           setequal(rownames(Rchk), traits) && setequal(colnames(Rchk), traits))
}

outputs_complete <- function() {
  file_ok(DONE_FLAG, 10) && file_ok(REPORT, 200) && file_ok(QC_TAB, 200) &&
    file_ok(SUM_OUT, 200) && file_ok(LEAD_OUT, 200) &&
    file_ok(UK_FG, 1e6) && file_ok(KJ_FG, 1e6)
}

snp_map_ok <- function(common) {
  if (!file_ok(MAP_PATH, 1000)) return(FALSE)
  sm <- tryCatch(qs::qread(MAP_PATH), error = function(e) NULL)
  if (is.null(sm)) return(FALSE)
  if (!is.data.table(sm)) sm <- as.data.table(sm)
  isTRUE(all(c("SNP", "CHR", "BP", "A1", "A2") %in% names(sm)) &&
           setequal(sm$SNP, common) && nrow(sm) == length(common))
}

# =============================================================================
assert_write_path(PIDFILE)
cat(sprintf("pid=%s\nstarted=%s\nscript=results/exp52_9/run_exp52_9_gls_Rreg.R\n",
            Sys.getpid(), now_str()), file = PIDFILE)
logmsg("START 52+9 genome-wide gls_Rreg; pid=", Sys.getpid(),
       " freeze=GFA seed 101 k=30 plus_resid=FALSE cond_num=1e3")
logmsg("Checkpoints present: Z_kj=", file_ok(CK_Z_KJ, 10),
       " R_kj=", file_ok(CK_R_KJ, 10),
       " gls_uk=", file_ok(CK_GLS_UK, 10),
       " gls_kj=", file_ok(CK_GLS_KJ, 10),
       " qc=", file_ok(CK_QC, 10))
write_status("START", "Preflight and policy checks. First step: Karjalainen 61-trait Z (build or reuse), then Karjalainen 61x61 R, then GLS both sides, then QC/cross-cohort.")

frozen_block <- c(FROZEN_Z_KJ, FROZEN_Z_UK, FROZEN_MAP, FROZEN_R_UK, FROZEN_R_KJ)
logmsg("Frozen paths will not be written: ", paste(basename(frozen_block), collapse = ", "))
if (normalizePath(RDIR, winslash = "/", mustWork = FALSE) ==
    normalizePath(file.path(PROJ, "results/week1/fgwas/gls_Rreg"), winslash = "/", mustWork = FALSE)) {
  stop("Refusing to write R into the frozen week1 gls_Rreg directory")
}

if (outputs_complete()) {
  logmsg("GLS_RREG_DONE.txt and key outputs already present; idempotent exit")
  write_status("DONE", "Idempotent skip: 52+9 gls_Rreg already finished.")
  quit(save = "no", status = 0)
}

need <- c(COMMON, TRAIT_N, Z_UK_PATH, MAP_PATH, R_UK_PATH, FIT_RDS, F_HAT_TSV,
          FROZEN_Z_KJ, LD_DIR)
missing <- need[!file.exists(need)]
if (length(missing)) stop("Missing required inputs:\n", paste(missing, collapse = "\n"))

traitN <- fread(TRAIT_N)
stopifnot("trait_code" %in% names(traitN))
if ("Glycine" %in% traitN$trait_code || "Gly" %in% traitN$trait_code) {
  stop("Glycine found in trait_N.tsv; this 52+9 run must exclude it")
}
traits <- traitN$trait_code
if (length(traits) != 61) stop("Expected 61 traits, got ", length(traits))
if (!all(LOW_H2Z %in% traits)) stop("Low-h2Z flag traits missing from panel: ",
                                    paste(setdiff(LOW_H2Z, traits), collapse = ","))

N_ukb_vec <- setNames(as.numeric(traitN$ukbb_N), traits)
N_ukb_vec[!is.finite(N_ukb_vec)] <- N_UKB
N_kj_vec <- setNames(as.numeric(traitN$karjalainen_N), traits)
kj_n_missing <- names(N_kj_vec)[!is.finite(N_kj_vec)]
if (length(kj_n_missing)) {
  logmsg("trait_N has empty Karj N for ", paste(kj_n_missing, collapse = ","),
         "; falling back to ", N_KARJ_FALLBACK,
         " (known first-5000 sampler quirk, not a QC failure)")
  N_kj_vec[!is.finite(N_kj_vec)] <- N_KARJ_FALLBACK
}
logmsg("Panel n_traits=", length(traits),
       " includes all 61; glycine absent; low-h2Z kept (flag only): ",
       paste(LOW_H2Z, collapse = ","))
logmsg("Karj N for R estimation: pyruvate=", N_kj_vec[["Pyruvate"]],
       " fallback_applied=", paste(kj_n_missing, collapse = ","),
       " n_fallback=", length(kj_n_missing))

common <- fread(COMMON, header = FALSE)$V1
if (length(common) != 1202568) {
  stop("Expected 1,202,568 common SNPs, got ", length(common))
}
logmsg("Common SNP panel n=", length(common))

if (!snp_map_ok(common)) {
  logmsg("snp_map.qs missing or broken; rebuilding from frozen map into results/exp52_9/cache/ only")
  if (!file_ok(FROZEN_MAP, 100)) stop("Missing frozen SNP map (read-only): ", FROZEN_MAP)
  map52 <- qs::qread(FROZEN_MAP)
  if (!is.data.table(map52)) map52 <- as.data.table(map52)
  stopifnot(all(c("SNP", "CHR", "BP", "A1", "A2") %in% names(map52)))
  map <- map52[match(common, SNP)]
  if (any(is.na(map$SNP))) stop("Frozen SNP map missing common SNPs: ", sum(is.na(map$SNP)))
  assert_write_path(MAP_PATH)
  qs::qsave(map, MAP_PATH)
  logmsg("Rebuilt snp_map.qs n=", nrow(map))
} else {
  logmsg("Reusing existing snp_map.qs (not overwritten)")
}
smap <- qs::qread(MAP_PATH)
if (!is.data.table(smap)) smap <- as.data.table(smap)
if (!identical(smap$SNP, common)) smap <- smap[match(common, SNP)]
stopifnot(identical(smap$SNP, common), all(c("SNP", "CHR", "BP", "A1", "A2") %in% names(smap)))

# --------------------------- 1. Karjalainen 61-trait Z
build_Z_karjalainen <- function() {
  if (valid_kj_z(traits, common)) {
    logmsg("Reuse existing 61-trait Karjalainen Z cache ", Z_KJ_PATH)
    return(invisible(TRUE))
  }
  write_status("Z_KARJALAINEN", "Building 61-trait Karjalainen Z on the 52+9 common SNP panel (reuse frozen 52 columns read-only; add 9 expansion traits).")
  logmsg("Building Karjalainen 61-trait Z cache (does not exist or is invalid)")
  if (!file_ok(FROZEN_Z_KJ, 1e6)) stop("Missing frozen 52-trait Karjalainen Z (read-only): ", FROZEN_Z_KJ)
  logmsg("Loading frozen 52-trait Karjalainen Z (read-only) from ", FROZEN_Z_KJ)
  z52 <- qs::qread(FROZEN_Z_KJ)
  Z52 <- z52$Z
  miss_common <- setdiff(common, rownames(Z52))
  if (length(miss_common)) stop("61-trait common SNPs missing from frozen Karj Z: ", length(miss_common))

  frozen_traits <- intersect(traits, colnames(Z52))
  new_traits <- setdiff(traits, colnames(Z52))
  logmsg("Reusing frozen Karj Z columns n=", length(frozen_traits),
         "; building new expansion columns n=", length(new_traits),
         ": ", paste(new_traits, collapse = ","))
  if (length(frozen_traits) != 52) stop("Expected to reuse 52 frozen Karj traits, got ", length(frozen_traits))
  if (length(new_traits) != 9) stop("Expected 9 expansion traits, got ", length(new_traits),
                                    " ", paste(new_traits, collapse = ","))

  Z <- matrix(NA_real_, length(common), length(traits), dimnames = list(common, traits))
  Z[, frozen_traits] <- Z52[common, frozen_traits, drop = FALSE]
  rm(Z52, z52); gc(verbose = FALSE)

  acc_map <- setNames(character(length(traits)), traits)
  acc_col <- intersect(c("karjalainen_accession", "karj_accession"), names(traitN))
  if (length(acc_col)) acc_map[traitN$trait_code] <- as.character(traitN[[acc_col[1]]])
  if (file_ok(CROSS61, 20)) {
    cr <- fread(CROSS61)
    code_col <- intersect(c("trait_code", "code"), names(cr))[1]
    kj_col <- intersect(c("karjalainen_accession", "karj_accession"), names(cr))[1]
    if (!is.na(code_col) && !is.na(kj_col)) {
      acc_map[as.character(cr[[code_col]])] <- as.character(cr[[kj_col]])
    }
  }
  if (any(!nzchar(acc_map[new_traits]) | is.na(acc_map[new_traits]))) {
    stop("Missing Karjalainen accession for: ",
         paste(new_traits[!nzchar(acc_map[new_traits]) | is.na(acc_map[new_traits])], collapse = ","))
  }

  ref_alleles <- smap[, .(rsid = SNP, ra1 = A1, ra2 = A2)]
  idx <- setNames(seq_along(common), common)
  stats <- data.table(trait = new_traits, accession = acc_map[new_traits], n_present = 0L,
                      n_flip = 0L, n_complement = 0L, n_complement_flip = 0L,
                      n_mismatch = 0L, n_absent = 0L, schema = "", source = "exp9_karjalainen")
  for (j in seq_along(new_traits)) {
    tr <- new_traits[j]
    acc <- acc_map[tr]
    path <- file.path(KARJ9, paste0(acc, ".hm3.tsv.gz"))
    if (!file.exists(path)) path <- file.path(KARJ52, paste0(acc, ".hm3.tsv.gz"))
    if (!file.exists(path)) stop("Missing HM3 for expansion Karj trait ", tr, " acc=", acc)
    logmsg("Z Karj new ", tr, " ", acc, " from ", path)
    got <- read_hm3_z(path)
    stats$schema[j] <- got$kind
    dt <- got$dt[rsid %chin% common]
    dt[, row := idx[rsid]]
    dt <- merge(dt, ref_alleles, by = "rsid", all.x = TRUE)
    dt[, rel := allele_relation(ra1, ra2, a1, a2)]
    n_abs_ref <- sum(is.na(dt$ra1))
    dt <- dt[!is.na(ra1)]
    dt[rel %chin% c("flip", "complement_flip"), z := -z]
    keep <- dt$rel %chin% c("match", "flip", "complement", "complement_flip")
    dt_keep <- dt[keep]
    Z[dt_keep$row, tr] <- dt_keep$z
    stats$n_present[j] <- nrow(dt_keep)
    stats$n_flip[j] <- sum(dt$rel == "flip")
    stats$n_complement[j] <- sum(dt$rel == "complement")
    stats$n_complement_flip[j] <- sum(dt$rel == "complement_flip")
    stats$n_mismatch[j] <- sum(dt$rel == "mismatch") + n_abs_ref
    stats$n_absent[j] <- length(common) - nrow(dt_keep)
    rm(dt, dt_keep, got); gc(verbose = FALSE)
  }
  assert_write_path(ZSTATS_KJ)
  fwrite(stats, ZSTATS_KJ, sep = "\t")
  n_miss <- sum(!is.finite(Z))
  frac <- n_miss / length(Z)
  logmsg("Karj Z missing cells=", n_miss, " / ", length(Z), " frac=", signif(frac, 4))
  if (frac > 0.01) stop("Too many missing Karj Z cells")
  if (n_miss > 0) Z[!is.finite(Z)] <- 0
  assert_write_path(Z_KJ_PATH)
  qs::qsave(list(Z = Z, traits = traits, n_snp = length(common),
                 reused_frozen52 = frozen_traits, new_traits = new_traits,
                 built = now_str()), Z_KJ_PATH)
  assert_write_path(CK_Z_KJ)
  writeLines(paste0("Z_KARJ_OK n_snp=", length(common), " n_trait=", length(traits), " ", now_str()), CK_Z_KJ)
  logmsg("Wrote 61-trait Karjalainen Z cache ", Z_KJ_PATH)
  invisible(TRUE)
}
build_Z_karjalainen()

logmsg("Loading UKB and Karjalainen 61-trait Z caches")
zuk_obj <- qs::qread(Z_UK_PATH)
zkj_obj <- qs::qread(Z_KJ_PATH)
Zuk <- zuk_obj$Z
Zkj <- zkj_obj$Z
stopifnot(is.matrix(Zuk), is.matrix(Zkj))
if (!identical(colnames(Zuk), traits)) {
  if (!all(traits %in% colnames(Zuk))) stop("UKB Z missing traits")
  Zuk <- Zuk[, traits, drop = FALSE]
}
if (!identical(colnames(Zkj), traits)) {
  if (!all(traits %in% colnames(Zkj))) stop("Karj Z missing traits")
  Zkj <- Zkj[, traits, drop = FALSE]
}
if (!identical(rownames(Zuk), common)) {
  if (!all(common %in% rownames(Zuk))) stop("UKB Z missing common SNPs")
  Zuk <- Zuk[common, , drop = FALSE]
}
if (!identical(rownames(Zkj), common)) {
  if (!all(common %in% rownames(Zkj))) stop("Karj Z missing common SNPs")
  Zkj <- Zkj[common, , drop = FALSE]
}
stopifnot(identical(dim(Zuk), c(1202568L, 61L)), identical(dim(Zkj), c(1202568L, 61L)))
if ("Glycine" %in% colnames(Zuk) || "Glycine" %in% colnames(Zkj) ||
    "Gly" %in% colnames(Zuk) || "Gly" %in% colnames(Zkj)) {
  stop("Glycine present in a Z matrix")
}
logmsg("Z UKB dim=", paste(dim(Zuk), collapse = "x"),
       " NA/nonfinite frac=", signif(mean(!is.finite(Zuk)), 4),
       " ; Karj dim=", paste(dim(Zkj), collapse = "x"),
       " NA/nonfinite frac=", signif(mean(!is.finite(Zkj)), 4))

# --------------------------- 2. Karjalainen 61x61 official-style R
build_R_karjalainen <- function() {
  if (valid_kj_r(traits)) {
    logmsg("Reuse existing 61x61 Karjalainen R_ldsc_cond1e3 ", R_KJ_PATH)
    return(invisible(TRUE))
  }
  write_status("R_KARJALAINEN", "Estimating new 61x61 official-style Karjalainen R_ldsc (cond_num=1e3) in results/exp52_9/R/.")
  logmsg("Estimating Karjalainen 61x61 R_ldsc (does not exist or is invalid)")
  if (normalizePath(RDIR, winslash = "/", mustWork = FALSE) ==
      normalizePath(dirname(FROZEN_R_KJ), winslash = "/", mustWork = FALSE)) {
    stop("Refusing to write R into the frozen gls_Rreg directory")
  }
  ldobj <- read_ldscores(LD_DIR)
  ld <- ldobj$ld
  common_snps <- rownames(Zkj)
  ix <- match(common_snps, ld$SNP)
  l2 <- rep(NA_real_, length(common_snps))
  l2[!is.na(ix)] <- ld$L2[ix[!is.na(ix)]]
  mhc <- smap$CHR == MHC_CHR & smap$BP >= MHC_START & smap$BP <= MHC_END
  use <- is.finite(l2) & !mhc
  logmsg("Karj R_ldsc will use ", sum(use), " SNPs with L2 and outside MHC (of ", length(common_snps), ")")
  logmsg("Karj R uses per-trait N from trait_N.tsv (fallback 136016 for empty entries)")
  Rraw <- estimate_R_ldsc(Zkj[use, , drop = FALSE], l2[use], ldobj$M, N_kj_vec, "karjalainen", RDIR)
  obj <- condition_pd(Rraw, COND_NUM, corr = FALSE)
  write_R_matrix(obj$mat, R_KJ_PATH)
  meta_row <- data.table(
    side = "karjalainen", R_source = "R_ldsc", cond_num = COND_NUM,
    n_snp = sum(use), n_trait = 61,
    raw_cond = obj$raw_cond, pd_cond = obj$pd_cond,
    raw_min_ev = obj$raw_min_ev, pd_min_ev = obj$pd_min_ev,
    diag_min = obj$diag_min, diag_max = obj$diag_max,
    note = paste0("new 61x61 matrix; per-trait Karj N; fallback 136016 for ",
                  paste(kj_n_missing, collapse = ","),
                  "; frozen 52x52 R was not copied")
  )
  assert_write_path(R_META)
  if (file_ok(R_META, 20)) {
    oldm <- tryCatch(fread(R_META), error = function(e) NULL)
    if (!is.null(oldm) && "side" %in% names(oldm)) {
      oldm <- oldm[side != "karjalainen"]
      meta_row <- rbind(oldm, meta_row, fill = TRUE)
    }
  }
  fwrite(meta_row, R_META, sep = "\t")
  assert_write_path(CK_R_KJ)
  writeLines(paste0("R_KARJ_OK n_trait=61 cond=", signif(obj$pd_cond, 6), " ", now_str()), CK_R_KJ)
  logmsg("Wrote conditioned 61x61 Karj R ", R_KJ_PATH, " pd_cond=", signif(obj$pd_cond, 6))
  invisible(TRUE)
}
build_R_karjalainen()

Ruk <- read_R_matrix(R_UK_PATH)
Rkj <- read_R_matrix(R_KJ_PATH)
Ruk <- Ruk[traits, traits, drop = FALSE]
Rkj <- Rkj[traits, traits, drop = FALSE]
if (nrow(Ruk) != 61 || ncol(Ruk) != 61) stop("UKB R is not 61x61")
if (nrow(Rkj) != 61 || ncol(Rkj) != 61) stop("Karj R is not 61x61")
logmsg("R UKB dim=", paste(dim(Ruk), collapse = "x"), " cond=", signif(mat_cond(Ruk), 6),
       " ; R Karj dim=", paste(dim(Rkj), collapse = "x"), " cond=", signif(mat_cond(Rkj), 6))

# --------------------------- 3. Load frozen GFA seed 101 fit (F_hat * scale)
logmsg("Loading frozen GFA fit RDS (read-only): ", FIT_RDS)
fit <- readRDS(FIT_RDS)
if (is.null(fit$F_hat) || is.null(fit$scale)) {
  stop("Frozen fit RDS missing F_hat or scale; refusing to project from F_hat alone")
}
if (is.null(rownames(fit$F_hat))) rownames(fit$F_hat) <- traits
if (!identical(rownames(fit$F_hat), traits)) {
  if (!all(traits %in% rownames(fit$F_hat))) stop("Frozen F_hat is missing traits")
  logmsg("Reordering F_hat rows to trait_N order")
  ord <- match(traits, rownames(fit$F_hat))
  fit$F_hat <- fit$F_hat[ord, , drop = FALSE]
  if (!is.null(names(fit$scale)) && all(traits %in% names(fit$scale))) {
    fit$scale <- fit$scale[traits]
  } else if (length(fit$scale) == length(ord)) {
    fit$scale <- fit$scale[ord]
  }
}
if ("Glycine" %in% rownames(fit$F_hat) || "Gly" %in% rownames(fit$F_hat)) {
  stop("Glycine present in frozen F_hat")
}
kfac <- ncol(fit$F_hat)
logmsg("Frozen F k=", kfac, " n_trait=", nrow(fit$F_hat),
       " scale_length=", length(fit$scale),
       " method=", if (!is.null(fit$method)) fit$method else "NA")
if (kfac != EXPECTED_K) {
  stop("Expected frozen seed 101 k=", EXPECTED_K, ", got ", kfac)
}

# --------------------------- 4. Genome-wide GLS, both sides
do_gls_side <- function(side, Zfull, Ruse, out_path, ck_path) {
  if (file_ok(ck_path, 10) && file_ok(out_path, 1e6)) {
    logmsg("Reuse existing GLS ", side, " ", out_path)
    return(invisible(TRUE))
  }
  write_status(paste("GLS", side),
               paste0("Genome-wide gls_Rreg, Sigma=cohort R_ldsc_cond1e3, plus_resid=FALSE, k=", kfac, ", side=", side, "."))
  logmsg("GLS genome-wide gls_Rreg ", side, " n_snp=", nrow(Zfull), " k=", kfac,
         " plus_resid=FALSE")
  t0 <- proc.time()[3]
  ret <- run_gls(Zfull, fit, Ruse, plus_resid = FALSE)
  write_fgwas(ret, rownames(Zfull), out_path)
  lam <- sapply(seq_len(kfac), function(j) lambda_gc(ret$L[, j] / ret$S[, j]))
  logmsg("  ", side, " lambda F1-Fk: ", paste(signif(lam, 4), collapse = " / "),
         "  Sigma_cond=", signif(ret$Sigma_cond, 4),
         "  elapsed_s=", round(proc.time()[3] - t0))
  assert_write_path(ck_path)
  writeLines(paste0("GLS_OK ", side, " k=", kfac, " ", now_str(),
                    " lambda=", paste(signif(lam, 4), collapse = ",")), ck_path)
  rm(ret); gc(verbose = FALSE)
  invisible(TRUE)
}

if (any(!is.finite(Zuk))) Zuk[!is.finite(Zuk)] <- 0
if (any(!is.finite(Zkj))) Zkj[!is.finite(Zkj)] <- 0

do_gls_side("ukbb_eur", Zuk, Ruk, UK_FG, CK_GLS_UK)
do_gls_side("karjalainen", Zkj, Rkj, KJ_FG, CK_GLS_KJ)

# --------------------------- 5. QC and one-to-one cross-cohort (no retain/drop)
do_qc_cross <- function() {
  if (file_ok(CK_QC, 10) && file_ok(QC_TAB, 200) && file_ok(SUM_OUT, 200) &&
      file_ok(LEAD_OUT, 200) && file_ok(REPORT, 200)) {
    logmsg("Reuse existing QC / cross-cohort / report")
    return(invisible(TRUE))
  }
  write_status("QC_CROSS", "Computing lambda_GC, LDSC intercept/h2, and one-to-one cross-cohort metrics. No retain/drop decision.")
  logmsg("Reading EUR HM3 LD scores from ", LD_DIR)
  ldpack <- read_ldscores(LD_DIR)
  logmsg("LD scores n_snp=", ldpack$n_ld_snps, " M_5_50=", ldpack$M)
  logmsg("Factor-level LDSC uses Step 2 constants N_UKB=", N_UKB, " N_KARJ=", N_KARJ)

  pos <- merge(
    smap[, .(SNP, CHR_map = CHR, BP_map = BP)],
    ldpack$ld[, .(SNP, CHR_ld = CHR, BP_ld = BP, L2)],
    by = "SNP", all.x = TRUE, sort = FALSE
  )
  pos[, `:=`(
    CHR = fifelse(is.finite(CHR_ld), CHR_ld, CHR_map),
    BP  = fifelse(is.finite(BP_ld), BP_ld, BP_map)
  )]

  ldsc_one_file <- function(path, side, method, scope, Nuse) {
    outp <- file.path(OUT, sprintf("qc_ldsc_%s_%s_%s.tsv", side, method, scope))
    if (file_ok(outp, 200)) {
      logmsg("Skip existing LDSC QC ", method, " ", side, " ", scope)
      return(fread(outp))
    }
    logmsg("QC/LDSC ", method, " ", side, " ", scope, " N=", Nuse)
    dt <- fread(path)
    zcols <- paste0("Z_F", seq_len(kfac))
    if (!all(c("SNP", zcols) %in% names(dt))) {
      stop("Unexpected F-GWAS columns in ", path, ": ", paste(names(dt), collapse = ","))
    }
    dt <- merge(dt, pos[, .(SNP, CHR, BP, L2)], by = "SNP", all.x = TRUE, sort = FALSE)
    rows <- list()
    for (j in seq_len(kfac)) {
      z <- dt[[paste0("Z_F", j)]]
      row <- data.table(
        side = side, method = method, scope = scope, factor = paste0("F", j),
        n_snps = nrow(dt),
        lambda_gc = lambda_gc(z),
        mean_chisq = mean(z[is.finite(z)]^2),
        se_median = median(dt[[paste0("SE_F", j)]], na.rm = TRUE)
      )
      if (sum(is.finite(dt$L2) & is.finite(z)) >= 2000) {
        uni1 <- tryCatch(
          ldsc_univar(dt$SNP, z, dt$CHR, dt$BP, dt$L2, dt$L2, Nuse, ldpack$M,
                      drop_mhc = TRUE, two_step = FALSE),
          error = function(e) { logmsg("ldsc fail ", method, " ", side, " ", j, ": ", e$message); NULL }
        )
        uni2 <- tryCatch(
          ldsc_univar(dt$SNP, z, dt$CHR, dt$BP, dt$L2, dt$L2, Nuse, ldpack$M,
                      drop_mhc = TRUE, two_step = TRUE),
          error = function(e) NULL
        )
        row[, `:=`(
          n_ldsc = if (is.null(uni1)) NA else uni1$n_snp,
          ldsc_intercept = if (is.null(uni1)) NA else uni1$intercept,
          ldsc_intercept_se = if (is.null(uni1)) NA else uni1$intercept_se,
          ldsc_h2 = if (is.null(uni1)) NA else uni1$h2,
          ldsc_h2_se = if (is.null(uni1)) NA else uni1$h2_se,
          ldsc_mean_chisq = if (is.null(uni1)) NA else uni1$mean_chisq,
          ldsc_twostep_intercept = if (is.null(uni2)) NA else uni2$intercept,
          ldsc_twostep_h2 = if (is.null(uni2)) NA else uni2$h2
        )]
      } else {
        row[, `:=`(n_ldsc = NA_real_, ldsc_intercept = NA_real_,
                   ldsc_intercept_se = NA_real_, ldsc_h2 = NA_real_,
                   ldsc_h2_se = NA_real_, ldsc_mean_chisq = NA_real_,
                   ldsc_twostep_intercept = NA_real_, ldsc_twostep_h2 = NA_real_)]
      }
      rows[[j]] <- row
      logmsg("  ", side, " ", method, " ", scope, " F", j,
             " lambda=", signif(row$lambda_gc, 4),
             " intercept=", signif(row$ldsc_intercept, 4),
             " h2=", signif(row$ldsc_h2, 4))
    }
    out <- rbindlist(rows, fill = TRUE)
    assert_write_path(outp)
    fwrite(out, outp, sep = "\t")
    out
  }

  qc_new <- rbind(
    ldsc_one_file(UK_FG, "ukbb_eur", METHOD, "genomewide", N_UKB),
    ldsc_one_file(KJ_FG, "karjalainen", METHOD, "genomewide", N_KARJ)
  )
  assert_write_path(QC_TAB)
  fwrite(qc_new, QC_TAB, sep = "\t")

  h2_of <- function(sd, fac) {
    v <- qc_new[side == sd & method == METHOD & scope == "genomewide" & factor == fac, ldsc_h2]
    if (length(v) != 1 || !is.finite(v[1])) stop("Missing genome-wide ldsc_h2 for ", sd, " ", fac)
    as.numeric(v[1])
  }
  read_z <- function(path) {
    zcols <- paste0("Z_F", seq_len(kfac))
    hdr <- names(fread(path, nrows = 0))
    if (!all(c("SNP", zcols) %in% hdr)) stop("Unexpected F-GWAS columns in ", path, ": ", paste(hdr, collapse = ","))
    dt <- fread(path, select = c("SNP", zcols))
    if (!identical(dt$SNP, common)) {
      logmsg("SNP order differs from common in ", path, "; realigning")
      dt <- dt[match(common, SNP)]
      if (any(is.na(dt$SNP))) stop("F-GWAS missing common SNPs: ", path)
    }
    as.matrix(dt[, ..zcols])
  }
  logmsg("Reading new UKB F-GWAS ", UK_FG)
  Zu <- read_z(UK_FG)
  logmsg("Reading new Karj F-GWAS ", KJ_FG)
  Zk <- read_z(KJ_FG)

  sum_rows <- list()
  lead_rows <- list()
  for (j in seq_len(kfac)) {
    fac <- paste0("F", j)
    ps <- file.path(OUT, "partial", sprintf("cross_%s_%s.tsv", METHOD, fac))
    pl <- file.path(OUT, "partial", sprintf("leads_%s_%s.tsv", METHOD, fac))
    if (file.exists(ps) && file.info(ps)$size > 50 && file.exists(pl) && file.info(pl)$size > 50) {
      logmsg("skip cross ", fac, " (partial exists)")
      sum_rows[[j]] <- fread(ps)
      lead_rows[[j]] <- fread(pl)
      next
    }
    h2u <- h2_of("ukbb_eur", fac)
    h2k <- h2_of("karjalainen", fac)
    logmsg("metrics ", fac, " h2_ukb=", signif(h2u, 4), " h2_karj=", signif(h2k, 4))
    write_status(paste("CROSS", fac), paste0("Pearson / sign / leads / LDSC rg for ", fac, "."))
    res <- cross_metrics(common, Zu[, j], Zk[, j], smap, ldpack$ld, ldpack$M, h2u, h2k, fac, METHOD)
    assert_write_path(ps)
    assert_write_path(pl)
    fwrite(res$summary, ps, sep = "\t")
    fwrite(res$leads, pl, sep = "\t")
    sum_rows[[j]] <- res$summary
    lead_rows[[j]] <- res$leads
    logmsg("  ", fac, " pearson=", signif(res$summary$pearson_z, 4),
           " sign_|Z|>=8=", signif(res$summary$sign_concord_zukb_ge8, 4),
           " rg=", signif(res$summary$ldsc_rg, 4),
           " rg_se=", signif(res$summary$ldsc_rg_se, 4))
  }
  cross_dt <- rbindlist(sum_rows, fill = TRUE)
  lead_dt <- rbindlist(lead_rows, fill = TRUE)
  assert_write_path(SUM_OUT)
  assert_write_path(LEAD_OUT)
  fwrite(cross_dt, SUM_OUT, sep = "\t")
  fwrite(lead_dt, LEAD_OUT, sep = "\t")

  lines <- c(
    "# Independent 52+9 (61-trait) genome-wide gls_Rreg projection",
    "",
    paste0("Finished: ", now_str(), " (", TZ_LABEL, ")"),
    "",
    "## Scope",
    "- Independent analysis of the 61-trait panel (frozen 52 + 9 expansion). Glycine was not included.",
    "- Frozen GFA seed 101 was used as-is (k=30; confirmed identical in subspace to seed 202). GFA was not refit; kmax was not locked.",
    "- Official-style `gls_Rreg`: Sigma = cohort residual `R_ldsc` conditioned at 1e3, `plus_resid=FALSE`. Loadings came from the fit RDS as `F_hat * scale`, not from the TSV alone.",
    "- UKB 61-trait Z and 61x61 R were reused from `results/exp52_9/cache/` and `results/exp52_9/R/`.",
    "- Karjalainen 61-trait Z was built in this job (or reused if already valid): frozen 52-trait Karj Z read-only, plus 9 expansion traits from `inputs/sumstats/exp52_9/karjalainen/` with the same allele alignment as the UKB Z builder.",
    "- Karjalainen 61x61 `R_ldsc` was estimated in this job (or reused if already valid) with per-trait Karjalainen N from `trait_N.tsv`. Empty Karj N for Ala, Gln, and/or His (known first-5000 sampler quirk, not a QC failure) fell back to 136016 and was logged.",
    "- All 61 traits were kept, including Omega_3, DHA, XL_HDL_TG, XS_VLDL_L, and XS_VLDL_FC, which remain flagged because UKB univariate LDSC h2 Z did not pass the suggested Z>4 gate. No retain/drop factor-QC decision was made.",
    "- No GLEANR GRID and no disease F-GWAS / MR.",
    "- Frozen matched-52 / week1 cache / gls_Rreg / step2_ldclump / main GFA objects were not overwritten.",
    "",
    "## Projection",
    "",
    paste0("- Seed: ", SEED_PRIMARY, "; selected k: **", kfac, "**."),
    paste0("- SNPs: **", length(common), "** (52+9 common HM3 panel)."),
    paste0("- Traits: **61**."),
    paste0("- Method: ", METHOD, " (`GFA:::loadings_gls`, Sigma = cohort R, plus_resid=FALSE)."),
    "",
    "## Calibration and cross-cohort constants",
    "",
    "Factor-level LDSC (univariate h2/intercept and bivariate rg) used the same constants as matched-52 Step 2: **N_UKB=413897**, **N_KARJ=136016**, LD scores supplied by `LDSC_EUR_LD_DIR`, MHC chr6 25–34 Mb dropped, 1 Mb UKB distance-clumped leads (`CLUMP_KB_MAIN=1000`), LDSC rg with chr-block jackknife, h2 from the new genome-wide univar LDSC. Per-trait Karjalainen N was used only when estimating the 61x61 residual R, not for factor-level LDSC.",
    "",
    "## Calibration (genome-wide lambda_GC and LDSC intercept / h2)",
    "",
    "| Factor | UKB lambda_GC | UKB intercept | UKB h2 | Karj lambda_GC | Karj intercept | Karj h2 |",
    "| --- | --- | --- | --- | --- | --- | --- |"
  )
  for (j in seq_len(kfac)) {
    fac <- paste0("F", j)
    u <- qc_new[side == "ukbb_eur" & factor == fac][1]
    k <- qc_new[side == "karjalainen" & factor == fac][1]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s |",
      fac, fmt(u$lambda_gc), fmt(u$ldsc_intercept), fmt(u$ldsc_h2),
      fmt(k$lambda_gc), fmt(k$ldsc_intercept), fmt(k$ldsc_h2)
    ))
  }
  lines <- c(
    lines,
    "",
    "## Cross-cohort (same metrics as matched-52 Step 2)",
    "",
    "| Factor | Pearson Z | Sign \\|Z_UKB\\|≥8 | LDSC rg (SE) | Lead \\|Z\\|≥8 n (sign) |",
    "| --- | --- | --- | --- | --- |"
  )
  for (j in seq_len(kfac)) {
    fac <- paste0("F", j)
    nw <- cross_dt[factor == fac][1]
    ld <- lead_dt[factor == fac & lead_set == "ukb_clumped_1mb_zge8"][1]
    nlead <- if (nrow(ld) && !is.na(ld$n_leads[1])) as.character(as.integer(ld$n_leads[1])) else "NA"
    slead <- if (nrow(ld)) fmt(ld$sign_concord[1]) else "NA"
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s (SE %s) | %s (%s) |",
      fac, fmt(nw$pearson_z), fmt(nw$sign_concord_zukb_ge8),
      fmt(nw$ldsc_rg), fmt(nw$ldsc_rg_se), nlead, slead
    ))
  }
  lines <- c(
    lines,
    "",
    "## Full QC table",
    "",
    "```",
    paste(capture.output(print(qc_new)), collapse = "\n"),
    "```",
    "",
    "## Full cross-cohort summary",
    "",
    "```",
    paste(capture.output(print(cross_dt)), collapse = "\n"),
    "```",
    "",
    "## Leads",
    "",
    "```",
    paste(capture.output(print(lead_dt)), collapse = "\n"),
    "```",
    "",
    "## Paths",
    "",
    "- Frozen F / fit (read-only): `results/exp52_9/gfa/ukbb_ldclump_seed101_F_hat.tsv`, `.../ukbb_ldclump_seed101.rds`",
    "- Karj Z: `results/exp52_9/cache/Z_karjalainen.qs`",
    "- Karj R: `results/exp52_9/R/R_ldsc_cond1e3_karjalainen.tsv`",
    "- F-GWAS: `results/exp52_9/gls_Rreg/fgwas_<side>_gls_Rreg_genomewide.tsv.gz`",
    "- QC: `results/exp52_9/gls_Rreg/qc_lambda_ldsc_gls_Rreg.tsv`",
    "- Cross-cohort: `results/exp52_9/gls_Rreg/cross_cohort_gls_Rreg.tsv`",
    "- Leads: `results/exp52_9/gls_Rreg/cross_cohort_leads_gls_Rreg.tsv`",
    "",
    "## Not done in this job",
    "",
    "- No retain/drop factor-QC decision.",
    "- GLEANR GRID: not started.",
    "- Disease F-GWAS / MR: not started.",
    "- Frozen matched-52 objects were not overwritten."
  )
  assert_write_path(REPORT)
  writeLines(lines, REPORT)
  assert_write_path(CK_QC)
  writeLines(paste0("QC_CROSS_OK k=", kfac, " ", now_str()), CK_QC)
  invisible(TRUE)
}

do_qc_cross()

assert_write_path(DONE_FLAG)
writeLines(paste0(
  "DONE ", now_str(), "\n",
  "k=", kfac, "\n",
  "n_snp=", length(common), "\n",
  "n_trait=61\n",
  "method=gls_Rreg\n",
  "plus_resid=FALSE\n",
  "N_UKB_factor_ldsc=", N_UKB, "\n",
  "N_KARJ_factor_ldsc=", N_KARJ, "\n",
  "report=", REPORT, "\n",
  "factor_qc_decision=not_made\n",
  "gleanr=not_started\n",
  "disease_MR=not_started\n",
  "frozen52=not_overwritten\n"
), DONE_FLAG)
write_status("DONE", paste0("Wrote ", REPORT, " and ", DONE_FLAG, ". No retain/drop decision; disease MR and GLEANR were not started."))
logmsg("DONE 52+9 genome-wide gls_Rreg. k=", kfac, " report=", REPORT)
quit(save = "no", status = 0)
