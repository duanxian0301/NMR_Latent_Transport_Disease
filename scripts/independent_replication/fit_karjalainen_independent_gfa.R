# Phase S1: Karjalainen independent GFA
# Writes ONLY under INDEPENDENT_GFA_REPLICATION_20261119
# prepare_R_for_gfa() and mat_cond() are copied VERBATIM from
#   scripts/core/run_exp52_9_gfa.R (public path; copied from the official project script)

options(stringsAsFactors = FALSE, warn = 1)
TZ_LABEL <- "Asia/Shanghai"

args <- commandArgs(trailingOnly = TRUE)
PREFLIGHT_ONLY <- "--preflight-only" %in% args

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
NEWDIR <- file.path(PROJ, "INDEPENDENT_GFA_REPLICATION_20261119")
OUTDIR <- file.path(NEWDIR, "s1_fit")

# Official inputs: READ ONLY
TRAIT_N <- file.path(PROJ, "results/exp52_9/trait_N.tsv")
Z_PATH <- file.path(PROJ, "results/exp52_9/cache/Z_karjalainen.qs")
PRUNE_IN <- file.path(PROJ, "results/exp52_9/prune/g1000_indep.prune.in")
R_OFFICIAL <- file.path(PROJ, "results/exp52_9/R/R_ldsc_cond1e3_karjalainen.tsv")

# Write paths: new directory only
N_PATH <- file.path(NEWDIR, "S1_KARJALAINEN_N.tsv")
R_USED <- file.path(OUTDIR, "R_used_for_karj_independent_gfa.tsv")
R_META <- file.path(OUTDIR, "R_PREPARE_CONDITION.tsv")
PREFLIGHT_PATH <- file.path(OUTDIR, "S1_PREFLIGHT.txt")
STATUS_PATH <- file.path(NEWDIR, "S1_FIT_STATUS.md")
LOG_PATH <- file.path(OUTDIR, "s1_karjalainen_independent_gfa.log")
RDS_PATH <- file.path(OUTDIR, "karjalainen_independent_gfa_seed101.rds")
FHAT_PATH <- file.path(OUTDIR, "karjalainen_independent_F_hat.tsv")
FHAT_SCALED_PATH <- file.path(OUTDIR, "karjalainen_independent_F_hat_times_scale.tsv")
SUMMARY_PATH <- file.path(NEWDIR, "S1_FIT_SUMMARY.tsv")
PID_PATH <- file.path(OUTDIR, "s1_fit.pid")

P3_N <- c(Ala = 136016, Gln = 131686, His = 136016)
EXPECTED_GFA_VERSION <- "1.0.0.0470"
OUTER_SEED <- 101L

FORBIDDEN_WRITE_SUBSTR <- c(
  "results/exp52_9/gfa/",
  "results/exp52_9/gls_Rreg",
  "results/exp52_9/gfa_stability_snp_split",
  "results/exp52_9/cache/",
  "results/exp52_9/R/",
  "results/exp52_9/prune/",
  "results/exp52_9/manuscript",
  "results/exp52_9/trait_N.tsv",
  "results/exp52_9/run_exp52_9_gfa_factorgo.R",
  "gls_commonR_ukb",
  "P108_METHOD_GATE",
  "AGING_EXTENSION_20261119",
  "audit_followup_20261117",
  "PHASE2E_INDICATOR_NDD_LAVA",
  "Cross_Organ_Genetic_Programs_NDD"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

norm_path <- function(path) {
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE), error = function(e) path)
  gsub("\\\\", "/", p)
}

assert_write_paths_safe <- function(paths) {
  ok_root <- norm_path(NEWDIR)
  for (path in paths) {
    pn <- norm_path(path)
    if (!startsWith(tolower(pn), tolower(ok_root))) {
      stop("Refusing write path outside the new directory: ", path)
    }
    rel <- substring(pn, nchar(ok_root) + 1L)
    # After the new-dir prefix, still block accidental official subpaths
    for (bad in FORBIDDEN_WRITE_SUBSTR) {
      if (grepl(bad, pn, fixed = TRUE) && !grepl("INDEPENDENT_GFA_REPLICATION_20261119", pn, fixed = TRUE)) {
        stop("Refusing write path that looks official: ", path)
      }
    }
    if (grepl("INDEPENDENT_GFA_REPLICATION_20261119", pn, fixed = TRUE)) next
  }
  invisible(TRUE)
}

WRITE_PATHS <- c(
  N_PATH, R_USED, R_META, PREFLIGHT_PATH, STATUS_PATH, LOG_PATH,
  RDS_PATH, FHAT_PATH, FHAT_SCALED_PATH, SUMMARY_PATH, PID_PATH
)
assert_write_paths_safe(WRITE_PATHS)

# Open log (append so preflight and fit share one file)
LOG_CON <- file(LOG_PATH, open = "at", encoding = "UTF-8")

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")

logmsg <- function(...) {
  msg <- sprintf("[%s] %s\n", now_str(), paste0(..., collapse = ""))
  cat(msg)
  flush(stdout())
  try(flush.console(), silent = TRUE)
  if (exists("LOG_CON", inherits = FALSE) && isOpen(LOG_CON)) {
    cat(msg, file = LOG_CON)
    flush(LOG_CON)
  }
}

# --- verbatim helpers from official run_exp52_9_gfa_factorgo.R ---
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
# --- end verbatim helpers ---

write_utf8_lines <- function(path, lines) {
  assert_write_paths_safe(path)
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, useBytes = FALSE)
}

write_status <- function(state, extra = character()) {
  pid <- Sys.getpid()
  lines <- c(
    "# S1 Karjalainen independent GFA status",
    "",
    paste0("Updated: ", now_str(), " (", TZ_LABEL, ")"),
    "",
    paste0("**State:** ", state),
    "",
    paste0("- PID (Sys.getpid): ", pid),
    paste0("- Command: Rscript --vanilla ",
           "s1_fit/s1_karjalainen_independent_gfa.R",
           if (PREFLIGHT_ONLY) " --preflight-only" else ""),
    paste0("- Log: `", LOG_PATH, "`"),
    paste0("- Script: `", file.path(OUTDIR, "s1_karjalainen_independent_gfa.R"), "`"),
    paste0("- Protocol: `", file.path(NEWDIR, "S1_INDEPENDENT_GFA_PROTOCOL.md"), "`"),
    "",
    extra,
    "",
    "## Policy",
    "- Outputs only under INDEPENDENT_GFA_REPLICATION_20261119.",
    "- Official UKB seed 101/202, official R, official Z, trait_N.tsv, GLS, P108, Aging, and manuscript are not modified.",
    "- K is chosen by the package rule; K=30 sensitivity is NOT_RUN this launch.",
    "- S2 and scientific judgments 1-3 are NOT_RUN this launch.",
    "- set.seed(101) is for official call-pattern comparability; default init_fn contains set.seed(666)."
  )
  write_utf8_lines(STATUS_PATH, lines)
  writeLines(as.character(pid), PID_PATH)
}

close_log <- function() {
  if (exists("LOG_CON", inherits = FALSE) && isOpen(LOG_CON)) {
    try(close(LOG_CON), silent = TRUE)
  }
}

on_fatal <- function(e) {
  msg <- paste0("FAILED: ", conditionMessage(e))
  try(logmsg(msg), silent = TRUE)
  try(logmsg("Traceback follows"), silent = TRUE)
  try({
    tb <- paste(capture.output(traceback(e, 1L)), collapse = "\n")
    logmsg(tb)
  }, silent = TRUE)
  try(write_status(if (PREFLIGHT_ONLY) "PREFLIGHT_FAILED" else "FAILED", c(
    paste0("- Failure time: ", now_str(), " (", TZ_LABEL, ")"),
    paste0("- Error: ", conditionMessage(e)),
    "- gfa_fit was not treated as successful. No fabricated results.",
    if (PREFLIGHT_ONLY) "- This was preflight-only; gfa_fit was not launched." else NULL
  )), silent = TRUE)
  try(close_log(), silent = TRUE)
  quit(save = "no", status = 1)
}

run_preflight <- function() {
  logmsg("==== S1 preflight start ====")
  logmsg("PREFLIGHT_ONLY=", PREFLIGHT_ONLY, " pid=", Sys.getpid())

  logmsg("Loading packages")
  suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(qs)
    library(GFA)
  })
  ver_obj <- packageVersion("GFA")
  ver <- as.character(ver_obj)
  desc_ver <- as.character(packageDescription("GFA")$Version)
  logmsg("packageVersion(GFA)=", ver, " DESCRIPTION Version=", desc_ver)
  # packageVersion() prints 1.0.0.470 because numeric components drop the leading zero;
  # DESCRIPTION Version is 1.0.0.0470. These compare equal as package_version objects.
  version_ok <- isTRUE(ver_obj == package_version(EXPECTED_GFA_VERSION)) ||
    identical(desc_ver, EXPECTED_GFA_VERSION) ||
    identical(ver, EXPECTED_GFA_VERSION)
  if (!version_ok) {
    stop("GFA version is packageVersion=", ver, " DESCRIPTION=", desc_ver,
         " but required version is ", EXPECTED_GFA_VERSION)
  }
  if (!identical(desc_ver, EXPECTED_GFA_VERSION)) {
    stop("GFA DESCRIPTION Version is ", desc_ver, " but required is ", EXPECTED_GFA_VERSION)
  }

  if (!file.exists(TRAIT_N)) stop("Missing official trait_N.tsv: ", TRAIT_N)
  if (!file.exists(Z_PATH)) stop("Missing official Z cache: ", Z_PATH)
  if (!file.exists(PRUNE_IN)) stop("Missing official prune.in: ", PRUNE_IN)
  if (!file.exists(R_OFFICIAL)) stop("Missing official Karjalainen R: ", R_OFFICIAL)
  if (grepl("ukbb|ukb", basename(R_OFFICIAL), ignore.case = TRUE)) {
    stop("R_OFFICIAL path looks like a UKB matrix: ", R_OFFICIAL)
  }
  if (!grepl("karjalainen", basename(R_OFFICIAL), ignore.case = TRUE)) {
    stop("R_OFFICIAL basename is not the official Karjalainen cond R: ", R_OFFICIAL)
  }

  traitN <- fread(TRAIT_N)
  if (!("trait_code" %in% names(traitN))) stop("trait_N.tsv missing trait_code")
  if (!("karjalainen_N" %in% names(traitN))) stop("trait_N.tsv missing karjalainen_N")
  traits <- as.character(traitN$trait_code)
  if (length(traits) != 61L) stop("Expected 61 traits, got ", length(traits))
  if (anyDuplicated(traits)) stop("Duplicate trait codes in trait_N.tsv")
  logmsg("Trait order n=", length(traits), " first=", traits[1], " last=", traits[length(traits)])

  official_N <- setNames(as.numeric(traitN$karjalainen_N), traits)
  Nvec <- official_N
  for (tr in names(P3_N)) {
    if (!tr %in% traits) stop("P3 trait missing from official trait_N: ", tr)
    Nvec[[tr]] <- unname(P3_N[[tr]])
  }
  others <- setdiff(traits, names(P3_N))
  if (any(!is.finite(official_N[others]))) {
    stop("Official Karjalainen N missing/non-finite for: ",
         paste(others[!is.finite(official_N[others])], collapse = ","))
  }
  if (!isTRUE(all.equal(as.numeric(Nvec[others]), as.numeric(official_N[others]), tolerance = 0))) {
    stop("S1 N for non-Ala/Gln/His traits does not match official trait_N.tsv")
  }
  if (!isTRUE(all.equal(unname(Nvec["Ala"]), 136016)) ||
      !isTRUE(all.equal(unname(Nvec["Gln"]), 131686)) ||
      !isTRUE(all.equal(unname(Nvec["His"]), 136016))) {
    stop("Ala/Gln/His N are not the required P3 recovered values 136016/131686/136016")
  }
  if (length(Nvec) != 61L || any(!is.finite(Nvec)) || any(Nvec <= 0)) {
    stop("Nvec failed finite/positive/length-61 check")
  }
  mode_n <- as.numeric(names(sort(table(Nvec), decreasing = TRUE)[1]))
  if (isTRUE(all(Nvec == 136016))) {
    stop("Refusing N vector that sets every trait to 136016")
  }
  n_tab <- data.table(
    trait_order = seq_along(traits),
    trait_code = traits,
    official_trait_N_karjalainen_N = as.numeric(official_N[traits]),
    s1_N_used = as.numeric(Nvec[traits]),
    n_source = ifelse(traits %in% names(P3_N), "P3_recovered", "official_trait_N"),
    scale_sqrt_N = sqrt(as.numeric(Nvec[traits]))
  )
  assert_write_paths_safe(N_PATH)
  fwrite(n_tab, N_PATH, sep = "\t")
  logmsg("Wrote N table ", N_PATH, " Gln=", Nvec[["Gln"]], " mode_N=", mode_n,
         " min_N=", min(Nvec), " max_N=", max(Nvec))

  logmsg("Reading official prune.in (file order, not alphabetical sort)")
  prune_fread <- unique(as.character(fread(PRUNE_IN, header = FALSE)$V1))
  prune_lines <- unique(as.character(readLines(PRUNE_IN, warn = FALSE)))
  prune_lines <- prune_lines[nzchar(prune_lines)]
  if (!identical(prune_fread, prune_lines)) {
    stop("prune.in file order differs between fread and readLines; refusing to guess")
  }
  if (length(prune_fread) != 69975L) {
    stop("Expected 69975 prune SNPs, got ", length(prune_fread))
  }
  if (isTRUE(all(prune_fread == sort(prune_fread))) && !identical(prune_fread, sort(prune_fread))) {
    # no-op; just document if already alphabetical
  }
  logmsg("prune.in unique SNPs=", length(prune_fread),
         " identical_to_alpha_sort=", identical(prune_fread, sort(prune_fread)))

  logmsg("Reading official Karjalainen Z (read-only): ", Z_PATH)
  zobj <- qs::qread(Z_PATH)
  if (!is.list(zobj) || is.null(zobj$Z)) stop("Z_karjalainen.qs does not contain $Z")
  Z <- zobj$Z
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (!identical(as.integer(dim(Z)), c(1202568L, 61L))) {
    stop("Expected Z dim 1202568 x 61, got ", paste(dim(Z), collapse = "x"))
  }
  if (!identical(colnames(Z), traits)) {
    if (!all(traits %in% colnames(Z))) {
      stop("Z is missing trait columns: ", paste(setdiff(traits, colnames(Z)), collapse = ","))
    }
    logmsg("Reordering Z columns to official trait_N order")
    Z <- Z[, traits, drop = FALSE]
  }
  if (!identical(colnames(Z), traits)) stop("Z colnames still do not match trait_N order")

  missing_snps <- setdiff(prune_fread, rownames(Z))
  if (length(missing_snps) > 0) {
    stop("Not all prune.in SNPs are present in Z: missing n=", length(missing_snps),
         " e.g. ", paste(head(missing_snps, 5), collapse = ","))
  }
  prune_snps <- prune_fread
  Zfit <- Z[prune_snps, traits, drop = FALSE]
  if (!identical(as.integer(dim(Zfit)), c(69975L, 61L))) {
    stop("Expected Zfit dim 69975 x 61, got ", paste(dim(Zfit), collapse = "x"))
  }
  if (!identical(rownames(Zfit), prune_snps)) stop("Zfit row order is not prune.in file order")
  if (!identical(colnames(Zfit), traits)) stop("Zfit colnames are not official trait order")

  n_exact_zero <- sum(Zfit == 0, na.rm = TRUE)
  n_nonfinite <- sum(!is.finite(Zfit))
  zero_frac <- n_exact_zero / length(Zfit)
  logmsg("Zfit exact zeros=", n_exact_zero, " frac=", signif(zero_frac, 6),
         " nonfinite=", n_nonfinite, " (S0 reported 3124 zeros, ~0.073%)")

  rm(Z, zobj)
  gc()

  logmsg("Reading official Karjalainen cond R (read-only): ", R_OFFICIAL)
  Rkj <- read_R_matrix(R_OFFICIAL)
  if (!identical(as.integer(dim(Rkj)), c(61L, 61L))) {
    stop("Expected official R 61 x 61, got ", paste(dim(Rkj), collapse = "x"))
  }
  if (!identical(rownames(Rkj), traits) || !identical(colnames(Rkj), traits)) {
    if (!all(traits %in% rownames(Rkj)) || !all(traits %in% colnames(Rkj))) {
      stop("Official Karjalainen R does not contain all 61 traits")
    }
    Rkj <- Rkj[traits, traits, drop = FALSE]
  }
  if (!identical(rownames(Rkj), traits) || !identical(colnames(Rkj), traits)) {
    stop("Official Karjalainen R trait order could not be aligned")
  }
  if (max(abs(Rkj - t(Rkj))) > 1e-8) stop("Official Karjalainen R is not symmetric")

  params <- GFA:::gfa_default_parameters()
  if (!is.null(params$kmax)) {
    stop("params$kmax is not NULL at default_parameters(); refusing to change the official pattern. kmax=",
         params$kmax)
  }
  target_cond <- if (is.null(params$cond_num)) 1000 else params$cond_num
  logmsg("gfa_default_parameters: kmax=NULL cond_num=", target_cond)

  prepR <- prepare_R_for_gfa(Rkj, target_cond = target_cond)
  R_fit <- prepR$R
  params$cond_num <- prepR$cond_param
  if (!is.null(params$kmax)) {
    stop("params$kmax was altered during R preparation; refusing to continue")
  }
  if (!identical(rownames(R_fit), traits) || !identical(colnames(R_fit), traits)) {
    stop("Prepared R dimnames do not match official trait order")
  }
  if (max(abs(R_fit - t(R_fit))) > 1e-8) stop("Prepared R is not symmetric")
  ev_fit <- eigen(R_fit, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev_fit <= 0)) {
    stop("Prepared R is not positive definite; min ev=", min(ev_fit))
  }

  write_R_matrix(R_fit, R_USED)
  r_meta <- data.table(
    source_R = R_OFFICIAL,
    prepared_R = R_USED,
    n_trait = 61L,
    trait_order_source = "results/exp52_9/trait_N.tsv",
    cond_before = prepR$cond_before,
    cond_after = prepR$cond_after,
    target_cond = target_cond,
    params_cond_num = params$cond_num,
    params_kmax_before_gfa_fit = "NULL",
    min_ev_after = min(ev_fit),
    max_ev_after = max(ev_fit),
    is_symmetric = TRUE,
    is_positive_definite = TRUE,
    used_UKB_R = FALSE,
    official_R_used_Gln_fallback_N_136016 = TRUE,
    s1_N_Gln_for_scale = 131686,
    n_corrected_R_sensitivity = "NOT_RUN",
    prepared_at = now_str()
  )
  fwrite(r_meta, R_META, sep = "\t")
  logmsg("Prepared R written. cond_before=", signif(prepR$cond_before, 8),
         " cond_after=", signif(prepR$cond_after, 8),
         " params$cond_num=", signif(params$cond_num, 8))

  pf_lines <- c(
    "S1_PREFLIGHT=PASS",
    paste0("time=", now_str(), " (", TZ_LABEL, ")"),
    paste0("GFA_version=", ver),
    paste0("Z_dim=", paste(c(1202568, 61), collapse = "x")),
    paste0("Zfit_dim=", paste(dim(Zfit), collapse = "x")),
    paste0("n_train_snp=", length(prune_snps)),
    paste0("Zfit_row_order=prune.in_file_order"),
    paste0("n_exact_zero_Zfit=", n_exact_zero),
    paste0("zero_frac_Zfit=", signif(zero_frac, 6)),
    paste0("n_nonfinite_Zfit=", n_nonfinite),
    paste0("N_Ala=", Nvec[["Ala"]]),
    paste0("N_Gln=", Nvec[["Gln"]]),
    paste0("N_His=", Nvec[["His"]]),
    paste0("N_min=", min(Nvec)),
    paste0("N_max=", max(Nvec)),
    paste0("R_source=", R_OFFICIAL),
    paste0("R_cond_before=", prepR$cond_before),
    paste0("R_cond_after=", prepR$cond_after),
    paste0("params_cond_num=", params$cond_num),
    "params_kmax=NULL",
    paste0("N_table=", N_PATH),
    paste0("R_used=", R_USED)
  )
  write_utf8_lines(PREFLIGHT_PATH, pf_lines)
  logmsg("==== S1 preflight PASS ====")

  list(
    traits = traits,
    Nvec = Nvec,
    prune_snps = prune_snps,
    Zfit = Zfit,
    R_fit = R_fit,
    params = params,
    prepR = prepR,
    n_exact_zero = n_exact_zero,
    zero_frac = zero_frac,
    n_nonfinite = n_nonfinite,
    gfa_version = ver
  )
}

save_success <- function(pf, fit, t0, t1, warns) {
  traits <- pf$traits
  Fh <- fit$F_hat
  if (is.null(Fh) || !is.matrix(Fh) || ncol(Fh) < 1) {
    stop("GFA returned no factors")
  }
  if (is.null(rownames(Fh))) rownames(Fh) <- traits
  Fh <- Fh[traits, , drop = FALSE]
  k <- ncol(Fh)
  colnames(Fh) <- paste0("Karj_F", seq_len(k))
  scale_vec <- fit$scale
  if (is.null(scale_vec)) scale_vec <- sqrt(as.numeric(pf$Nvec[traits]))
  if (is.null(names(scale_vec))) names(scale_vec) <- traits
  scale_vec <- scale_vec[traits]
  Fh_scaled <- Fh * as.numeric(scale_vec)
  dimnames(Fh_scaled) <- dimnames(Fh)

  write.table(Fh, FHAT_PATH, sep = "\t", quote = FALSE, col.names = NA)
  write.table(Fh_scaled, FHAT_SCALED_PATH, sep = "\t", quote = FALSE, col.names = NA)

  elapsed_hours <- as.numeric(difftime(t1, t0, units = "hours"))
  method <- if (!is.null(fit$method)) fit$method else NA_character_
  kmax_used <- if (!is.null(fit$params$kmax)) fit$params$kmax else NA
  pve <- if (!is.null(fit$gfa_pve)) fit$gfa_pve else NA

  obj <- list(
    fit = fit,
    method = method,
    L_hat = if (!is.null(fit$L_hat)) fit$L_hat else NULL,
    F_hat = Fh,
    F_hat_times_scale = Fh_scaled,
    scale = scale_vec,
    params_input = pf$params,
    params_fit = if (!is.null(fit$params)) fit$params else NULL,
    R = pf$R_fit,
    gfa_pve = if (!is.null(fit$gfa_pve)) fit$gfa_pve else NULL,
    sessionInfo = sessionInfo(),
    seed_notes = paste(
      "Outer set.seed(101) was called immediately before gfa_fit for official",
      "call-pattern comparability, NOT as a claim of an independent random start.",
      "Default init_fn contains set.seed(666), so the outer seed may not randomize",
      "greedy initialization. This is a recorded limitation."
    ),
    N = pf$Nvec,
    training_snps = pf$prune_snps,
    package_version = pf$gfa_version,
    n_exact_zero_Zfit = pf$n_exact_zero,
    zero_frac_Zfit = pf$zero_frac,
    n_nonfinite_Zfit = pf$n_nonfinite,
    R_condition = list(
      source = R_OFFICIAL,
      cond_before = pf$prepR$cond_before,
      cond_after = pf$prepR$cond_after,
      params_cond_num = pf$params$cond_num,
      official_R_used_Gln_fallback_N_136016 = TRUE
    ),
    warnings = warns,
    started = format(t0, "%Y-%m-%d %H:%M:%S %z"),
    finished = format(t1, "%Y-%m-%d %H:%M:%S %z"),
    elapsed_hours = elapsed_hours,
    k = k,
    kmax_param = "NULL",
    kmax_used = kmax_used,
    tz = TZ_LABEL
  )
  saveRDS(obj, RDS_PATH)

  summary <- data.table(
    state = "COMPLETE",
    started = obj$started,
    finished = obj$finished,
    elapsed_hours = elapsed_hours,
    K = k,
    n_snp = nrow(pf$Zfit),
    n_trait = ncol(pf$Zfit),
    method = method,
    kmax_param = "NULL",
    kmax_used = kmax_used,
    cond_before = pf$prepR$cond_before,
    cond_after = pf$prepR$cond_after,
    params_cond_num = pf$params$cond_num,
    package_version = pf$gfa_version,
    N_Gln = unname(pf$Nvec[["Gln"]]),
    n_exact_zero_Zfit = pf$n_exact_zero,
    n_warning = length(warns),
    rds = RDS_PATH,
    F_hat = FHAT_PATH,
    F_hat_times_scale = FHAT_SCALED_PATH,
    log = LOG_PATH
  )
  fwrite(summary, SUMMARY_PATH, sep = "\t")

  logmsg("gfa_fit complete method=", method, " K=", k, " elapsed_hours=", signif(elapsed_hours, 6))
  if (length(warns)) logmsg("Warnings n=", length(warns), " : ", paste(warns, collapse = " | "))
  if (!is.null(fit$gfa_pve)) {
    logmsg("gfa_pve available; length=", length(fit$gfa_pve))
  }

  write_status("COMPLETE", c(
    paste0("- Start: ", obj$started, " (", TZ_LABEL, ")"),
    paste0("- Finish: ", obj$finished, " (", TZ_LABEL, ")"),
    paste0("- Elapsed hours: ", signif(elapsed_hours, 6)),
    paste0("- K (package rule, not forced): ", k),
    paste0("- method: ", method),
    paste0("- RDS: `", RDS_PATH, "`"),
    paste0("- F_hat: `", FHAT_PATH, "`"),
    paste0("- F_hat * scale: `", FHAT_SCALED_PATH, "`"),
    paste0("- Summary: `", SUMMARY_PATH, "`"),
    "- K=30 sensitivity, S2, and scientific judgments 1-3 remain NOT_RUN."
  ))
}

# ---------------- main ----------------
tryCatch({
  logmsg("S1 script start cwd=", getwd(), " R=", R.version.string)
  pf <- run_preflight()
  if (PREFLIGHT_ONLY) {
    logmsg("Preflight-only mode: not calling gfa_fit, not writing IN_PROGRESS fit status")
    close_log()
    quit(save = "no", status = 0)
  }

  write_status("IN_PROGRESS", c(
    paste0("- Start: ", now_str(), " (", TZ_LABEL, ")"),
    "- Preflight: PASS",
    paste0("- Zfit: ", paste(dim(pf$Zfit), collapse = "x")),
    paste0("- N_Gln (scale): ", pf$Nvec[["Gln"]]),
    paste0("- R cond_before/after: ", signif(pf$prepR$cond_before, 8), " / ", signif(pf$prepR$cond_after, 8)),
    paste0("- params$cond_num: ", signif(pf$params$cond_num, 8)),
    "- params$kmax: NULL (package will set 2*p)",
    "- Official files are not write targets.",
    "- Fit is running; this file will be updated to COMPLETE or FAILED."
  ))

  logmsg("Seed note: calling set.seed(101) for official call-pattern comparability.")
  logmsg("Limitation: default init_fn contains set.seed(666); outer seed may not randomize greedy initialization.")
  logmsg("Calling gfa_fit Z=", paste(dim(pf$Zfit), collapse = "x"),
         " kmax=NULL cond_num=", signif(pf$params$cond_num, 8))
  if (!is.null(pf$params$kmax)) stop("kmax is not NULL immediately before gfa_fit")
  if (!identical(names(pf$Nvec), pf$traits) && !identical(names(pf$Nvec), NULL)) {
    # names should be traits; gfa_fit accepts a vector of length p
  }
  N_for_fit <- as.numeric(pf$Nvec[pf$traits])
  names(N_for_fit) <- pf$traits

  warns <- character()
  t0 <- Sys.time()
  set.seed(OUTER_SEED)
  fit <- withCallingHandlers(
    gfa_fit(Z_hat = pf$Zfit, N = N_for_fit, R = pf$R_fit, params = pf$params),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      logmsg("WARNING: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  t1 <- Sys.time()
  save_success(pf, fit, t0, t1, warns)
  close_log()
  quit(save = "no", status = 0)
}, error = function(e) on_fatal(e))
