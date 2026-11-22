#!/usr/bin/env Rscript
# Independent 52+9 (61-trait) UKB GFA + FactorGo fit.
# Outputs ONLY under results/exp52_9/. Does NOT overwrite frozen matched-52
# objects, week1 cache, step2_ldclump, gls_Rreg, or main GFA / F-GWAS.
# Does NOT start disease F-GWAS / MR, and does NOT start a GLEANR GRID.
#
# Pipeline (resumable):
#   1. Build 61-trait UKB Z on the 52+9 common SNP panel (reuse frozen 52
#      columns by SNP subset; add 9 expansion traits with the same allele
#      alignment). Glycine is not included. Five UKB univariate LDSC
#      h2 Z<=4 traits are kept and only flagged.
#   2. Estimate a new 61x61 official-style R_ldsc and condition at 1e3
#      in this folder. Do not copy the frozen 52x52 R.
#   3. plink --indep-pairwise 1000 kb 100 0.1 on the 61-trait common SNPs
#      (MHC dropped; same target window as Step 2).
#   4. gfa_fit on UKB, regularized R, seed 101, kmax left at package default
#      (not locked at 5).
#   5. FactorGo two seeds (101, 202), official -k = n_traits = 61, SVD init.

suppressPackageStartupMessages({
  library(GFA)
  library(data.table)
  library(qs)
  library(Matrix)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
OUT  <- file.path(PROJ, "results/exp52_9")
CACHE <- file.path(OUT, "cache")
RDIR <- file.path(OUT, "R")
PRUNE_DIR <- file.path(OUT, "prune")
GFA_DIR <- file.path(OUT, "gfa")
FG_DIR <- file.path(OUT, "factorgo")

FROZEN_Z <- file.path(PROJ, "results/week1/fgwas/cache/Z_ukbb_eur.qs")
FROZEN_MAP <- file.path(PROJ, "results/week1/fgwas/cache/snp_map.qs")
FROZEN_R <- file.path(PROJ, "results/week1/fgwas/gls_Rreg/R_ldsc_cond1e3_ukbb_eur.tsv")
COMMON <- file.path(OUT, "common_snps.txt")
TRAIT_N <- file.path(OUT, "trait_N.tsv")
LDSC_SUM <- file.path(OUT, "ldsc_univariate/UKBB_EUR_UNIVAR_LDSC_SUMMARY.tsv")
UK52 <- file.path(PROJ, "inputs/sumstats/ukbb_eur")
UK9  <- file.path(PROJ, "inputs/sumstats/exp52_9/ukbb_eur")
CROSS61 <- file.path(PROJ, "inputs/crosswalk/exp52_9_trait_gcst_list.tsv")

LD_DIR <- Sys.getenv("LDSC_EUR_LD_DIR", unset = file.path(PROJ, "external/ldsc/eur_w_ld_chr"))
PLINK <- Sys.getenv("PLINK_EXECUTABLE", unset = "plink")
BFILE <- Sys.getenv("PLINK_EUR_BFILE", unset = file.path(PROJ, "external/1000g/g1000_eur"))
FG_EXE <- Sys.getenv("FACTORGO_EXECUTABLE", unset = "factorgo")

LOW_H2Z <- c("Omega_3", "DHA", "XL_HDL_TG", "XS_VLDL_L", "XS_VLDL_FC")
N_UKB_DEFAULT <- 413897
MHC_CHR <- 6L
MHC_START <- 25e6
MHC_END <- 34e6
COND_NUM <- 1000
SEED_PRIMARY <- 101L
SEED_STABILITY <- 202L
PRUNE_WINDOW_KB <- 1000
PRUNE_STEP <- 100
PRUNE_R2_START <- 0.1
TARGET_N_MIN <- 15000
TARGET_N_MAX <- 80000
TZ_LABEL <- "Asia/Shanghai"

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PRUNE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GFA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FG_DIR, recursive = TRUE, showWarnings = FALSE)

Z_PATH <- file.path(CACHE, "Z_ukbb_eur.qs")
MAP_PATH <- file.path(CACHE, "snp_map.qs")
ZSTATS <- file.path(CACHE, "Z_build_stats.tsv")
FLAG_PATH <- file.path(OUT, "LOW_H2Z_FLAG.tsv")
R_PAIR <- file.path(RDIR, "R_ldsc_pairs_ukbb_eur.tsv")
R_RAW <- file.path(RDIR, "R_ldsc_raw_ukbb_eur.tsv")
R_COND <- file.path(RDIR, "R_ldsc_cond1e3_ukbb_eur.tsv")
R_META <- file.path(RDIR, "R_meta.tsv")
PRUNE_IN <- file.path(PRUNE_DIR, "g1000_indep.prune.in")
PRUNE_META <- file.path(PRUNE_DIR, "prune_meta.tsv")
PRUNE_EXTRACT <- file.path(PRUNE_DIR, "extract_universe.txt")
F_101 <- file.path(GFA_DIR, "ukbb_ldclump_seed101_F_hat.tsv")
RDS_101 <- file.path(GFA_DIR, "ukbb_ldclump_seed101.rds")
F_202 <- file.path(GFA_DIR, "ukbb_ldclump_seed202_F_hat.tsv")
RDS_202 <- file.path(GFA_DIR, "ukbb_ldclump_seed202.rds")
FIT_META <- file.path(GFA_DIR, "fit_meta.tsv")
R_USED <- file.path(GFA_DIR, "R_used_for_gfa_fit_seed101.tsv")
STATUS <- file.path(OUT, "status.md")
REPORT <- file.path(OUT, "EXP52_9_FIT_REPORT.md")
DONE_FLAG <- file.path(OUT, "FIT_DONE.txt")
PIDFILE <- file.path(OUT, "exp52_9_fit.pid")
CK_Z <- file.path(OUT, "checkpoint_Z.txt")
CK_R <- file.path(OUT, "checkpoint_R.txt")
CK_PRUNE <- file.path(OUT, "checkpoint_prune.txt")
CK_GFA <- file.path(OUT, "checkpoint_gfa_seed101.txt")
CK_FG <- file.path(OUT, "checkpoint_factorgo.txt")

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
write_status <- function(step, detail = "") {
  lines <- c(
    "# 52+9 UKB GFA + FactorGo — live status",
    "",
    paste0("Updated: ", now_str(), " (", TZ_LABEL, ")"),
    "",
    paste0("**Step:** ", step),
    "",
    detail,
    "",
    "## Policy",
    "- Independent 52+9 analysis. Outputs only under `results/exp52_9/`.",
    "- Frozen matched-52 objects, week1 cache, step2_ldclump, gls_Rreg, and main GFA / F-GWAS are not overwritten.",
    "- 61 traits are all kept, including five UKB univariate LDSC h2 Z that did not pass Z>4; they are flagged, not dropped.",
    "- Glycine is not included.",
    "- k is not locked at 5.",
    "- Residual R is a new 61x61 official-style R_ldsc (cond_num=1e3), not the frozen 52x52 R.",
    "- No GLEANR GRID, no disease F-GWAS / MR, no genome-wide metabolite GLS in this job.",
    "",
    "- Script: `results/exp52_9/run_exp52_9_gfa_factorgo.R`",
    "- Log: `results/exp52_9/exp52_9_rscript.log`"
  )
  con <- file(STATUS, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}
mat_cond <- function(M) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  max(ev) / max(min(ev), .Machine$double.eps)
}
prepare_R_for_gfa <- function(R, target_cond = 1000) {
  cn0 <- mat_cond(R)
  R_use <- R
  if (cn0 > target_cond) {
    logmsg("R cond=", signif(cn0, 6), " exceeds target ", target_cond,
           "; re-applying Matrix::nearPD(posd.tol=1/", target_cond, ")")
    pd <- Matrix::nearPD(R, corr = FALSE, keepDiag = TRUE, posd.tol = 1 / target_cond)
    R_use <- as.matrix(pd$mat)
    dimnames(R_use) <- dimnames(R)
    cn1 <- mat_cond(R_use)
    logmsg("  after nearPD cond=", signif(cn1, 6))
  } else {
    cn1 <- cn0
  }
  cond_param <- target_cond
  if (cn1 > target_cond) {
    cond_param <- cn1 * 1.01
    logmsg("  nearPD still slightly above ", target_cond,
           "; setting params$cond_num=", signif(cond_param, 6),
           " only so check_R accepts the official-style 1e3 R")
  }
  list(R = R_use, cond_param = cond_param, cond_before = cn0, cond_after = mat_cond(R_use))
}
write_R_matrix <- function(M, path) {
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
run_plink_prune <- function(extract_path, r2, out_prefix) {
  log_path <- paste0(out_prefix, ".plink.log")
  args <- c(
    "--bfile", BFILE,
    "--extract", extract_path,
    "--indep-pairwise", as.character(PRUNE_WINDOW_KB), "kb", as.character(PRUNE_STEP),
    sprintf("%g", r2),
    "--out", out_prefix
  )
  logmsg("plink ", paste(args, collapse = " "))
  rc <- system2(PLINK, args, stdout = log_path, stderr = log_path)
  prune_in <- paste0(out_prefix, ".prune.in")
  list(rc = rc, log = log_path, prune_in = prune_in, r2 = r2, args = args)
}
parse_plink_counts <- function(log_path) {
  if (!file.exists(log_path)) return(list(loaded = NA, after_extract = NA, log_tail = ""))
  lines <- readLines(log_path, warn = FALSE)
  loaded <- NA_integer_
  after <- NA_integer_
  m1 <- grep("variants loaded from .bim file", lines, value = TRUE, ignore.case = TRUE)
  if (length(m1)) loaded <- as.integer(gsub("[^0-9]", "", sub("variants.*", "", m1[1])))
  m2 <- grep("variants remaining after main filters|variants remaining after --extract",
             lines, value = TRUE, ignore.case = TRUE)
  if (length(m2)) after <- as.integer(gsub("[^0-9]", "", sub("variants.*", "", m2[1])))
  list(loaded = loaded, after_extract = after, log_tail = paste(tail(lines, 20), collapse = " | "))
}

# =============================================================================
cat(sprintf("pid=%s\nstarted=%s\nscript=results/exp52_9/run_exp52_9_gfa_factorgo.R\n",
            Sys.getpid(), now_str()), file = PIDFILE)
logmsg("START 52+9 UKB GFA + FactorGo; pid=", Sys.getpid())
logmsg("Checkpoints present: Z=", file_ok(CK_Z, 10), " R=", file_ok(CK_R, 10),
       " prune=", file_ok(CK_PRUNE, 10), " gfa=", file_ok(CK_GFA, 10))
write_status("START", "Preflight and policy checks.")

# Refuse to write into frozen locations
frozen_block <- c(FROZEN_Z, FROZEN_MAP, FROZEN_R,
                  file.path(PROJ, "results/week1/fgwas/step2_ldclump/frozen_ukbb_ldclump_seed101_F_hat.tsv"),
                  file.path(PROJ, "results/week1/fgwas/step2_ldclump/R_used_for_gfa_fit_seed101.tsv"))
logmsg("Frozen paths will not be written: ", paste(basename(frozen_block), collapse = ", "))
if (file.exists(DONE_FLAG) && file_ok(F_101, 100) && file_ok(file.path(FG_DIR, "FactorGo_ukbb_ldclump_seed101.Zm.tsv.gz"), 100) &&
    file_ok(file.path(FG_DIR, "FactorGo_ukbb_ldclump_seed202.Zm.tsv.gz"), 100)) {
  logmsg("FIT_DONE already present and GFA/FactorGo outputs exist; exiting idempotently")
  quit(status = 0)
}

traitN <- fread(TRAIT_N)
stopifnot("trait_code" %in% names(traitN))
if ("Glycine" %in% traitN$trait_code || "Gly" %in% traitN$trait_code) {
  stop("Glycine found in trait_N.tsv; this 52+9 run must exclude it")
}
traits <- traitN$trait_code
if (length(traits) != 61) stop("Expected 61 traits, got ", length(traits))
if (!all(LOW_H2Z %in% traits)) stop("Low-h2Z flag traits missing from panel: ",
                                    paste(setdiff(LOW_H2Z, traits), collapse = ","))
Nvec <- setNames(as.numeric(traitN$ukbb_N), traits)
Nvec[!is.finite(Nvec)] <- N_UKB_DEFAULT
logmsg("Panel n_traits=", length(traits), " includes all 61; glycine absent; low-h2Z kept: ",
       paste(LOW_H2Z, collapse = ","))

# Flag the five UKB LDSC h2 Z that did not pass Z>4; do not drop them.
if (file_ok(LDSC_SUM, 50)) {
  sumtab <- fread(LDSC_SUM)
  keep_cols <- intersect(c("trait_code", "h2", "h2_se", "h2_z", "h2z_ge_4", "lambda_gc", "intercept"), names(sumtab))
  flag <- sumtab[trait_code %in% LOW_H2Z, ..keep_cols]
  flag[, action := "kept_not_dropped"]
  flag[, reason := "UKB univariate LDSC h2 Z did not pass suggested Z>4"]
} else {
  flag <- data.table(trait_code = LOW_H2Z, action = "kept_not_dropped",
                     reason = "UKB univariate LDSC h2 Z did not pass suggested Z>4; summary tsv not found at runtime")
}
fwrite(flag, FLAG_PATH, sep = "\t")

common <- fread(COMMON, header = FALSE)$V1
if (length(common) < 200000) stop("common SNP panel too small: ", length(common))
logmsg("Common SNP panel n=", length(common))

# --------------------------- 1. Z cache (UKB, 61 traits)
build_Z <- function() {
  if (file_ok(CK_Z, 10) && file_ok(Z_PATH, 1e6) && file_ok(MAP_PATH, 1000)) {
    logmsg("Reuse existing 61-trait Z cache")
    return(invisible(TRUE))
  }
  write_status("Z", "Building 61-trait UKB Z on the 52+9 common SNP panel.")
  if (!file_ok(FROZEN_Z, 1e6)) stop("Missing frozen 52-trait Z (read-only): ", FROZEN_Z)
  if (!file_ok(FROZEN_MAP, 100)) stop("Missing frozen SNP map (read-only): ", FROZEN_MAP)
  logmsg("Loading frozen 52-trait Z and SNP map (read-only)")
  z52 <- qs::qread(FROZEN_Z)
  Z52 <- z52$Z
  map52 <- qs::qread(FROZEN_MAP)
  if (!is.data.table(map52)) map52 <- as.data.table(map52)
  stopifnot(all(c("SNP", "CHR", "BP", "A1", "A2") %in% names(map52)))
  miss_common <- setdiff(common, rownames(Z52))
  if (length(miss_common)) stop("61-trait common SNPs missing from frozen Z: ", length(miss_common))
  map <- map52[match(common, SNP)]
  if (any(is.na(map$SNP))) stop("SNP map missing common SNPs: ", sum(is.na(map$SNP)))
  qs::qsave(map, MAP_PATH)

  frozen_traits <- intersect(traits, colnames(Z52))
  new_traits <- setdiff(traits, colnames(Z52))
  logmsg("Reusing frozen Z columns n=", length(frozen_traits),
         "; building new expansion columns n=", length(new_traits),
         ": ", paste(new_traits, collapse = ","))
  if (length(frozen_traits) != 52) stop("Expected to reuse 52 frozen traits, got ", length(frozen_traits))
  if (length(new_traits) != 9) stop("Expected 9 expansion traits, got ", length(new_traits),
                                    " ", paste(new_traits, collapse = ","))

  Z <- matrix(NA_real_, length(common), length(traits), dimnames = list(common, traits))
  Z[, frozen_traits] <- Z52[common, frozen_traits, drop = FALSE]
  rm(Z52, z52); gc(verbose = FALSE)

  acc_map <- setNames(character(length(traits)), traits)
  # prefer explicit accessions from trait_N if present
  acc_col <- intersect(c("ukbb_eur_accession", "ukbb_accession"), names(traitN))
  if (length(acc_col)) acc_map[traitN$trait_code] <- as.character(traitN[[acc_col[1]]])
  if (file_ok(CROSS61, 20)) {
    cr <- fread(CROSS61)
    code_col <- intersect(c("trait_code", "code"), names(cr))[1]
    uk_col <- intersect(c("ukbb_eur_accession", "ukbb_accession"), names(cr))[1]
    if (!is.na(code_col) && !is.na(uk_col)) {
      acc_map[as.character(cr[[code_col]])] <- as.character(cr[[uk_col]])
    }
  }
  if (any(!nzchar(acc_map[new_traits]))) stop("Missing UKB accession for: ",
                                              paste(new_traits[!nzchar(acc_map[new_traits])], collapse = ","))

  ref_alleles <- map[, .(rsid = SNP, ra1 = A1, ra2 = A2)]
  idx <- setNames(seq_along(common), common)
  stats <- data.table(trait = new_traits, accession = acc_map[new_traits], n_present = 0L,
                      n_flip = 0L, n_complement = 0L, n_complement_flip = 0L,
                      n_mismatch = 0L, n_absent = 0L, schema = "", source = "exp9")
  for (j in seq_along(new_traits)) {
    tr <- new_traits[j]
    acc <- acc_map[tr]
    path <- file.path(UK9, paste0(acc, ".hm3.tsv.gz"))
    if (!file.exists(path)) path <- file.path(UK52, paste0(acc, ".hm3.tsv.gz"))
    if (!file.exists(path)) stop("Missing HM3 for expansion trait ", tr, " acc=", acc)
    logmsg("Z new ", tr, " ", acc, " from ", path)
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
  fwrite(stats, ZSTATS, sep = "\t")
  n_miss <- sum(!is.finite(Z))
  frac <- n_miss / length(Z)
  logmsg("Z missing cells=", n_miss, " / ", length(Z), " frac=", signif(frac, 4))
  if (frac > 0.01) stop("Too many missing Z cells")
  if (n_miss > 0) Z[!is.finite(Z)] <- 0
  qs::qsave(list(Z = Z, traits = traits, n_snp = length(common),
                 reused_frozen52 = frozen_traits, new_traits = new_traits,
                 built = now_str()), Z_PATH)
  writeLines(paste0("Z_OK n_snp=", length(common), " n_trait=", length(traits), " ", now_str()), CK_Z)
  logmsg("Wrote 61-trait Z cache ", Z_PATH)
  invisible(TRUE)
}
build_Z()
zobj <- qs::qread(Z_PATH)
Zuk <- zobj$Z
stopifnot(is.matrix(Zuk), ncol(Zuk) == 61, all(colnames(Zuk) == traits) || setequal(colnames(Zuk), traits))
if (!identical(colnames(Zuk), traits)) Zuk <- Zuk[, traits, drop = FALSE]
smap <- qs::qread(MAP_PATH)
if (!is.data.table(smap)) smap <- as.data.table(smap)
stopifnot(identical(smap$SNP, rownames(Zuk)) || setequal(smap$SNP, rownames(Zuk)))
if (!identical(smap$SNP, rownames(Zuk))) smap <- smap[match(rownames(Zuk), SNP)]
logmsg("Z dim=", paste(dim(Zuk), collapse = "x"), " NA/nonfinite frac=", signif(mean(!is.finite(Zuk)), 4))

# --------------------------- 2. 61x61 official-style R
build_R <- function() {
  if (file_ok(CK_R, 10) && file_ok(R_COND, 200) && file_ok(R_RAW, 200)) {
    Rchk <- read_R_matrix(R_COND)
    if (nrow(Rchk) == 61 && ncol(Rchk) == 61 && setequal(rownames(Rchk), traits)) {
      logmsg("Reuse existing 61x61 R_ldsc_cond1e3")
      return(invisible(TRUE))
    }
    logmsg("Existing R is not 61x61 or trait set differs; recomputing in the new path")
  }
  write_status("R", "Estimating new 61x61 official-style R_ldsc (cond_num=1e3) in results/exp52_9/R/.")
  # Safety: never write the frozen 52 R path; refuse if someone pointed RDIR there
  if (normalizePath(RDIR, mustWork = FALSE) == normalizePath(dirname(FROZEN_R), mustWork = FALSE)) {
    stop("Refusing to write R into the frozen gls_Rreg directory")
  }
  ldobj <- read_ldscores(LD_DIR)
  ld <- ldobj$ld
  common_snps <- rownames(Zuk)
  ix <- match(common_snps, ld$SNP)
  l2 <- rep(NA_real_, length(common_snps))
  l2[!is.na(ix)] <- ld$L2[ix[!is.na(ix)]]
  # drop MHC for R estimation, same as Step 1/2
  mhc <- smap$CHR == MHC_CHR & smap$BP >= MHC_START & smap$BP <= MHC_END
  use <- is.finite(l2) & !mhc
  logmsg("R_ldsc will use ", sum(use), " SNPs with L2 and outside MHC (of ", length(common_snps), ")")
  Rraw <- estimate_R_ldsc(Zuk[use, , drop = FALSE], l2[use], ldobj$M, Nvec, "ukbb_eur", RDIR)
  obj <- condition_pd(Rraw, COND_NUM, corr = FALSE)
  write_R_matrix(obj$mat, R_COND)
  meta <- data.table(
    side = "ukbb_eur", R_source = "R_ldsc", cond_num = COND_NUM,
    n_snp = sum(use), n_trait = 61,
    raw_cond = obj$raw_cond, pd_cond = obj$pd_cond,
    raw_min_ev = obj$raw_min_ev, pd_min_ev = obj$pd_min_ev,
    diag_min = obj$diag_min, diag_max = obj$diag_max,
    note = "new 61x61 matrix; frozen 52x52 R was not copied"
  )
  fwrite(meta, R_META, sep = "\t")
  writeLines(paste0("R_OK n_trait=61 cond=", signif(obj$pd_cond, 6), " ", now_str()), CK_R)
  logmsg("Wrote conditioned 61x61 R ", R_COND, " pd_cond=", signif(obj$pd_cond, 6))
  invisible(TRUE)
}
build_R()
Ruk <- read_R_matrix(R_COND)
Ruk <- Ruk[traits, traits, drop = FALSE]
if (nrow(Ruk) != 61) stop("R is not 61x61")
logmsg("R dim=", paste(dim(Ruk), collapse = "x"), " cond=", signif(mat_cond(Ruk), 6))

# --------------------------- 3. LD prune
do_prune <- function() {
  if (file_ok(CK_PRUNE, 10) && file_ok(PRUNE_IN, 100) && file_ok(PRUNE_META, 20)) {
    meta <- fread(PRUNE_META)
    n <- if ("n_pruned" %in% names(meta)) as.integer(meta$n_pruned[1]) else length(readLines(PRUNE_IN))
    logmsg("Reuse existing prune: n=", n)
    return(invisible(TRUE))
  }
  write_status("PRUNE", "plink --indep-pairwise 1000 kb 100 0.1 on 61-trait common SNPs (not p-value clump).")
  universe_ids <- rownames(Zuk)
  mhc_map <- smap$CHR == MHC_CHR & smap$BP >= MHC_START & smap$BP <= MHC_END
  after_mhc <- setdiff(universe_ids, smap$SNP[mhc_map])
  logmsg("Universe Z=", length(universe_ids), " after snp_map MHC drop=", length(after_mhc),
         " MHC dropped=", sum(smap$SNP %in% universe_ids & mhc_map))
  logmsg("Reading g1000_eur.bim for intersection and MHC")
  bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1, 2, 4),
               col.names = c("CHR", "SNP", "BP"))
  bim[, CHR := as.integer(gsub("^chr", "", as.character(CHR), ignore.case = TRUE))]
  bim_mhc <- bim[CHR == MHC_CHR & BP >= MHC_START & BP <= MHC_END, SNP]
  in_bim <- after_mhc[after_mhc %in% bim$SNP]
  in_bim <- setdiff(in_bim, bim_mhc)
  logmsg("After intersect g1000_eur and bim-MHC drop: n=", length(in_bim))
  if (length(in_bim) < TARGET_N_MIN) stop("Prune universe too small after MHC/g1000 filters: ", length(in_bim))
  fwrite(data.table(SNP = in_bim), PRUNE_EXTRACT, col.names = FALSE)
  rm(bim); gc(verbose = FALSE)
  run_one <- function(r2) {
    prefix <- file.path(PRUNE_DIR, sprintf("g1000_indep_r2_%s", gsub("\\.", "p", sprintf("%g", r2))))
    res <- run_plink_prune(PRUNE_EXTRACT, r2, prefix)
    counts <- parse_plink_counts(res$log)
    if (!file.exists(res$prune_in)) {
      stop("plink prune failed (no prune.in). rc=", res$rc, " tail: ", counts$log_tail)
    }
    snps <- fread(res$prune_in, header = FALSE)$V1
    logmsg("prune r2=", r2, " n=", length(snps),
           " plink_loaded=", counts$loaded, " after_extract=", counts$after_extract, " rc=", res$rc)
    list(snps = snps, r2 = r2, res = res, counts = counts)
  }
  used <- run_one(PRUNE_R2_START)
  flag_note <- sprintf("start --indep-pairwise %d kb %d %g", PRUNE_WINDOW_KB, PRUNE_STEP, PRUNE_R2_START)
  if (length(used$snps) > TARGET_N_MAX) {
    logmsg("n_pruned=", length(used$snps), " > ", TARGET_N_MAX, "; tightening to r2=0.05")
    used <- run_one(0.05)
    flag_note <- paste0(flag_note, sprintf(" -> tightened --indep-pairwise %d kb %d 0.05", PRUNE_WINDOW_KB, PRUNE_STEP))
  } else if (length(used$snps) < TARGET_N_MIN) {
    logmsg("n_pruned=", length(used$snps), " < ", TARGET_N_MIN, "; loosening to r2=0.2")
    used <- run_one(0.2)
    flag_note <- paste0(flag_note, sprintf(" -> loosened --indep-pairwise %d kb %d 0.2", PRUNE_WINDOW_KB, PRUNE_STEP))
  }
  n <- length(used$snps)
  if (n < 1000) stop("Pruned SNP count implausibly small: ", n)
  if (n < TARGET_N_MIN || n > TARGET_N_MAX) {
    logmsg("WARNING: pruned n=", n, " is outside the 15k-80k target after adjustment; proceeding with this set")
  }
  fwrite(data.table(SNP = used$snps), PRUNE_IN, col.names = FALSE)
  meta <- data.table(
    window_kb = PRUNE_WINDOW_KB, step = PRUNE_STEP, r2 = used$r2,
    flags = sprintf("--indep-pairwise %d kb %d %g", PRUNE_WINDOW_KB, PRUNE_STEP, used$r2),
    flag_note = flag_note, n_pruned = n,
    plink_loaded = used$counts$loaded, plink_after_extract = used$counts$after_extract,
    plink_log = used$res$log, extract_path = PRUNE_EXTRACT, prune_in = PRUNE_IN,
    panel = "exp52_9_61trait"
  )
  fwrite(meta, PRUNE_META, sep = "\t")
  writeLines(paste0("PRUNE_OK n=", n, " r2=", used$r2, " ", now_str()), CK_PRUNE)
  logmsg("PRUNE done. flags=", meta$flags[1], " n=", n)
  invisible(TRUE)
}
logmsg("BEGIN LD prune (plink --indep-pairwise 1000 kb 100 0.1)")
do_prune()
if (!file_ok(PRUNE_IN, 100)) stop("LD prune did not produce a usable prune.in at ", PRUNE_IN)
prune_snps <- unique(as.character(fread(PRUNE_IN, header = FALSE)$V1))
prune_snps <- prune_snps[prune_snps %in% rownames(Zuk)]
if (length(prune_snps) < 1000) stop("prune_snps too small after intersecting with Z: ", length(prune_snps))
logmsg("Using pruned SNPs n=", length(prune_snps))

# --------------------------- 4. GFA seed 101
do_fit <- function(seed, f_path, rds_path, ck_path = NULL, r_used_path = NULL) {
  if (file_ok(f_path, 100) && file_ok(rds_path, 1000)) {
    Fh <- as.matrix(read.delim(f_path, check.names = FALSE, row.names = 1))
    logmsg("Reuse existing fit seed=", seed, " k=", ncol(Fh), " path=", f_path)
    return(invisible(TRUE))
  }
  write_status("GFA", paste0("gfa_fit on UKB pruned SNPs, 61 traits, seed=", seed,
                             ", new 61x61 R_ldsc_cond1e3, kmax left at package default."))
  Zfit <- Zuk[prune_snps, traits, drop = FALSE]
  params <- GFA:::gfa_default_parameters()
  if (!is.null(params$kmax) && params$kmax <= 5) {
    logmsg("default kmax=", params$kmax, " is too low; setting kmax=20")
    params$kmax <- 20
  }
  logmsg("Leaving params$kmax=", if (is.null(params$kmax)) "NULL (package will set 2*p=122)" else as.character(params$kmax))
  prepR <- prepare_R_for_gfa(Ruk, target_cond = if (is.null(params$cond_num)) 1000 else params$cond_num)
  R_fit <- prepR$R
  params$cond_num <- prepR$cond_param
  if (!is.null(r_used_path)) write_R_matrix(R_fit, r_used_path)
  logmsg("R for gfa_fit written; cond_param=", signif(params$cond_num, 6),
         " cond_after=", signif(prepR$cond_after, 6))
  logmsg("Calling gfa_fit seed=", seed, " Z=", paste(dim(Zfit), collapse = "x"), " ...")
  set.seed(seed)
  fit <- gfa_fit(Z_hat = Zfit, N = Nvec, R = R_fit, params = params)
  Fh <- fit$F_hat
  if (is.null(Fh) || ncol(Fh) < 1) stop("GFA returned no factors for seed=", seed)
  if (is.null(rownames(Fh))) rownames(Fh) <- traits
  Fh <- Fh[traits, , drop = FALSE]
  colnames(Fh) <- paste0("F", seq_len(ncol(Fh)))
  write.table(Fh, f_path, sep = "\t", quote = FALSE, col.names = NA)
  saveRDS(fit, rds_path)
  k <- ncol(Fh)
  logmsg("gfa_fit seed=", seed, " method=", if (!is.null(fit$method)) fit$method else "NA",
         " k=", k)
  kmax_used <- if (!is.null(fit$params$kmax)) fit$params$kmax else NA
  meta <- data.table(
    seed = seed, k = k, n_snp = nrow(Zfit), n_trait = ncol(Zfit),
    kmax_param = if (is.null(params$kmax)) "NULL->2p" else as.character(params$kmax),
    kmax_used = kmax_used, cond_num = params$cond_num,
    method = if (!is.null(fit$method)) fit$method else NA_character_,
    finished = now_str()
  )
  if (file_ok(FIT_META, 20)) {
    oldm <- tryCatch(fread(FIT_META), error = function(e) NULL)
    if (!is.null(oldm)) {
      oldm <- oldm[as.integer(seed) != as.integer(meta$seed)]
      meta <- rbind(oldm, meta, fill = TRUE)
    }
  }
  fwrite(meta, FIT_META, sep = "\t")
  if (!is.null(ck_path)) writeLines(paste0("FIT_OK seed=", seed, " k=", k, " ", now_str()), ck_path)
  invisible(TRUE)
}
do_fit(SEED_PRIMARY, F_101, RDS_101, CK_GFA, R_USED)
Fh <- as.matrix(read.delim(F_101, check.names = FALSE, row.names = 1))
logmsg("Primary GFA k=", ncol(Fh))

# --------------------------- 5. FactorGo two seeds
if (!exists("prune_snps", inherits = FALSE) || length(prune_snps) < 1000) {
  stop("Refusing to write FactorGo inputs: prune_snps is not ready")
}
fg_z <- file.path(FG_DIR, "FactorGo_ukbb_ldclump.Z.tsv")
fg_n <- file.path(FG_DIR, "FactorGo_ukbb_ldclump.N.tsv")
if (!(file_ok(fg_z, 1e6) && file_ok(fg_n, 20))) {
  logmsg("Writing FactorGo Z/N inputs after prune n=", length(prune_snps))
  Zfit <- Zuk[prune_snps, traits, drop = FALSE]
  fwrite(data.table(SNP = rownames(Zfit), Zfit), fg_z, sep = "\t")
  writeLines(c("N", as.character(as.numeric(Nvec[traits]))), fg_n)
} else {
  logmsg("Reuse existing FactorGo Z/N inputs")
  Zfit <- Zuk[prune_snps, traits, drop = FALSE]
}
FG_KMAX <- length(traits)
run_factorgo <- function(seed) {
  tag <- sprintf("FactorGo_ukbb_ldclump_seed%d", seed)
  prefix <- file.path(FG_DIR, tag)
  zm_path <- paste0(prefix, ".Zm.tsv.gz")
  fac_path <- paste0(prefix, ".factor.tsv.gz")
  if (file_ok(zm_path, 100) && file_ok(fac_path, 20)) {
    logmsg("Reuse FactorGo seed=", seed)
    return(invisible(TRUE))
  }
  write_status(paste("FACTORGO seed", seed),
               paste0("factorgo -k ", FG_KMAX, " on n_snp=", nrow(Zfit), ", seed=", seed))
  logmsg("RUN factorgo seed=", seed, " -k ", FG_KMAX)
  t0 <- Sys.time()
  rc <- system2(FG_EXE, args = c(
    fg_z, fg_n,
    "-k", as.character(FG_KMAX),
    "-s", as.character(seed),
    "--init-factor", "svd",
    "--max-iter", "10000",
    "--elbo-tol", "1e-3",
    "--rate", "100",
    "-o", prefix
  ), stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt"))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  writeLines(c(paste("factorgo", fg_z, fg_n, "-k", FG_KMAX, "-s", seed, "--init-factor svd"),
               paste("rc", rc, "elapsed_min", signif(elapsed, 4))),
             paste0(prefix, ".command.txt"))
  if (!identical(as.integer(rc), 0L) || !file_ok(zm_path, 100)) {
    err <- if (file.exists(paste0(prefix, ".stderr.txt")))
      paste(readLines(paste0(prefix, ".stderr.txt"), warn = FALSE), collapse = "\n") else ""
    stop("FactorGo failed seed=", seed, " rc=", rc, " stderr_tail=", substr(err, max(1, nchar(err) - 1500), nchar(err)))
  }
  logmsg("OK FactorGo seed=", seed, " elapsed_min=", signif(elapsed, 4))
  invisible(TRUE)
}
run_factorgo(SEED_PRIMARY)
run_factorgo(SEED_STABILITY)
writeLines(paste0("FACTORGO_OK seeds=101,202 kmax=", FG_KMAX, " ", now_str()), CK_FG)

# Optional GFA seed 202 after FactorGo, so the primary path is ready first.
tryCatch({
  do_fit(SEED_STABILITY, F_202, RDS_202, NULL, NULL)
}, error = function(e) logmsg("Optional GFA seed 202 failed or was skipped: ", e$message))

pmeta <- if (file_ok(PRUNE_META)) fread(PRUNE_META) else NULL
fmeta <- if (file_ok(FIT_META)) fread(FIT_META)[seed == SEED_PRIMARY][1] else NULL
report <- c(
  "# Independent 52+9 (61-trait) UKB GFA + FactorGo",
  "",
  paste0("Finished: ", now_str(), " (", TZ_LABEL, ")"),
  "",
  "## Scope",
  "- Independent analysis of the 61-trait panel (frozen 52 + 9 expansion). Glycine was not included.",
  "- All 61 traits were kept, including Omega_3, DHA, XL_HDL_TG, XS_VLDL_L, and XS_VLDL_FC, which are flagged in `LOW_H2Z_FLAG.tsv` because UKB univariate LDSC h2 Z did not pass the suggested Z>4 gate.",
  "- SNP panel: 61-trait common HM3 SNPs, then `plink --indep-pairwise 1000 kb 100 0.1` (MHC dropped).",
  "- Residual R: new 61x61 official-style `R_ldsc` conditioned at 1e3, written under `results/exp52_9/R/`. The frozen 52x52 R was not copied.",
  "- GFA: seed 101, kmax left at the package default (2*p = 122), not locked at 5.",
  "- FactorGo: official `-k` = 61, SVD init, seeds 101 and 202.",
  "- No GLEANR GRID, no disease F-GWAS / MR, and no genome-wide metabolite GLS in this job.",
  "- Frozen matched-52 objects were not overwritten.",
  "",
  "## Results",
  if (!is.null(pmeta)) paste0("- Prune flags: `", pmeta$flags[1], "`; n_pruned = **", pmeta$n_pruned[1], "**.") else "- Prune metadata missing.",
  paste0("- GFA seed 101 k = **", ncol(Fh), "**."),
  paste0("- FactorGo requested k = **", FG_KMAX, "**; outputs under `results/exp52_9/factorgo/`."),
  "",
  "## Files",
  "- Z cache: `results/exp52_9/cache/Z_ukbb_eur.qs`",
  "- R: `results/exp52_9/R/R_ldsc_cond1e3_ukbb_eur.tsv`",
  "- GFA F: `results/exp52_9/gfa/ukbb_ldclump_seed101_F_hat.tsv`",
  "- FactorGo: `results/exp52_9/factorgo/FactorGo_ukbb_ldclump_seed101` and `..._seed202`",
  "- Low-h2Z flag: `results/exp52_9/LOW_H2Z_FLAG.tsv`"
)
con <- file(REPORT, open = "w", encoding = "UTF-8")
writeLines(report, con)
close(con)
writeLines(c(
  paste0("done=", now_str()),
  "panel=exp52_9_61trait",
  paste0("gfa_k=", ncol(Fh)),
  paste0("factorgo_k_requested=", FG_KMAX),
  paste0("n_pruned=", length(prune_snps)),
  "gleanr=not_started",
  "disease_MR=not_started",
  "frozen52=not_overwritten"
), DONE_FLAG)
write_status("DONE", paste0("GFA k=", ncol(Fh), "; FactorGo two seeds finished; n_pruned=", length(prune_snps), "."))
logmsg("DONE 52+9 GFA + FactorGo. GFA k=", ncol(Fh), " FactorGo kmax=", FG_KMAX)
quit(status = 0)
