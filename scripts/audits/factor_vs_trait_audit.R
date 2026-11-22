#!/usr/bin/env Rscript
# P2C: dominance, factor vs dominant-trait, factor-lead vs 61-trait locus audit.
# Reuses existing GLS factor GWAS and Z caches. Does not overwrite GFA/GLS/LAVA/Z/R/manuscript.
# Language rule: never call "factor significant, trait not" a new independent causal locus.

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
EXP  <- file.path(PROJ, "results/exp52_9")
OUT  <- file.path(PROJ, "audit_followup_20261117/P2_FINAL")
STAGE <- Sys.getenv("NMR_GFA_AUDIT_MIRROR", unset = OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)

RETAIN <- c("F9","F11","F13","F14","F15","F17","F18","F20","F22","F23","F24","F25","F26")
CLUMP_KB_MAIN <- 1000
P_GWAS <- 5e-8
Z_GWAS <- qnorm(P_GWAS / 2, lower = FALSE)  # ~5.45131
SHADOW_Z <- 4  # pre-specified dominant-shadow |Z| cutoff
MHC_CHR <- 6; MHC_START <- 25e6; MHC_END <- 34e6

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
      paste0(..., collapse = ""), "\n", sep = "")
  flush.console()
}
stage <- function(path) {
  if (file.exists(path)) file.copy(path, file.path(STAGE, basename(path)), overwrite = TRUE)
  path
}
p_from_z <- function(z) 2 * pnorm(-abs(z))

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
  is.finite(chr) & is.finite(bp) & chr == MHC_CHR & bp >= MHC_START & bp <= MHC_END
}

# ---- 1. Loadings (GLS matrix) ----
logmsg("START P2C")
Fpath <- file.path(OUT, "FULL_GFA_LOADINGS_61x30.tsv")
if (!file.exists(Fpath)) stop("Missing ", Fpath, "; run export first")
Fdt <- fread(Fpath)
traits <- Fdt$trait
Fgls <- as.matrix(Fdt[, -1, with = FALSE])
rownames(Fgls) <- traits
stopifnot(all(RETAIN %in% colnames(Fgls)), length(traits) == 61)

sc <- fread(file.path(OUT, "FULL_GFA_SCALE_VECTOR.tsv"))
scale_nuniq <- uniqueN(round(sc$scale, 10))
scale_const <- scale_nuniq == 1
logmsg("scale unique (rounded 1e-10)=", scale_nuniq, " constant=", scale_const,
       " median_scale=", median(sc$scale))

# F_hat for |loading|>=0.10 secondary count (S4a-style threshold is only meaningful here)
Fh_path <- file.path(OUT, "FULL_GFA_LOADINGS_61x30_Fhat.tsv")
Fh <- as.matrix(fread(Fh_path)[, -1, with = FALSE])
rownames(Fh) <- traits

base_rows <- list()
for (fac in RETAIN) {
  v <- as.numeric(Fgls[, fac])
  ss <- sum(v^2)
  pshare <- (v^2) / ss
  names(pshare) <- traits
  i_dom <- which.max(abs(v))
  vh <- as.numeric(Fh[, fac])
  base_rows[[fac]] <- data.table(
    factor = fac,
    loading_matrix = "F_hat_times_scale",
    scale_constant_across_traits = scale_const,
    dominant_trait = traits[i_dom],
    dominant_loading_scaled = v[i_dom],
    dominant_loading_Fhat = vh[i_dom],
    max_loading_sq_share = pshare[i_dom],
    concentration_HHI = sum(pshare^2),
    effective_n_traits_participation_ratio = 1 / sum(pshare^2),
    n_traits_abs_ge_0.10_Fhat = sum(abs(vh) >= 0.10),
    n_traits_abs_ge_0.10_scaled = sum(abs(v) >= 0.10),
    n_traits_share_ge_0.01 = sum(pshare >= 0.01),
    top3_share = sum(sort(pshare, decreasing = TRUE)[1:3]),
    formula_share = "p_j = loading_j^2 / sum(loading^2); effective_n = 1/sum(p_j^2); concentration_HHI = sum(p_j^2)",
    note_0.10 = "Primary share metrics use F_hat*scale (GLS). The 0.10 threshold is reported on raw F_hat (S4a-style) and on the scaled matrix; if scale is ~sqrt(N) and nearly constant, shares are invariant to scale but the 0.10 count is not."
  )
}
base0 <- rbindlist(base_rows)
logmsg("Dominant traits: ", paste(base0$factor, base0$dominant_trait, sep = "=", collapse = ", "))

# ---- 2. Maps, Z, factor GWAS ----
logmsg("Loading snp_map and Z caches (read-only)")
smap <- qs::qread(file.path(EXP, "cache/snp_map.qs"))
if (is.data.frame(smap)) smap <- as.data.table(smap)
# normalize names
nm <- names(smap)
if (!"SNP" %in% nm) {
  cand <- intersect(c("rsid", "snp", "SNPID", "id"), nm)
  if (length(cand)) setnames(smap, cand[1], "SNP")
}
if (!"CHR" %in% names(smap)) {
  cand <- intersect(c("chr", "chrom", "CHR"), names(smap))
  if (length(cand)) setnames(smap, cand[1], "CHR")
}
if (!"BP" %in% names(smap)) {
  cand <- intersect(c("bp", "pos", "POS", "position", "BP"), names(smap))
  if (length(cand)) setnames(smap, cand[1], "BP")
}
smap[, CHR := as.integer(CHR)]
smap[, BP := as.integer(BP)]
logmsg("snp_map n=", nrow(smap), " cols=", paste(names(smap), collapse = ","))

Zuk <- qs::qread(file.path(EXP, "cache/Z_ukbb_eur.qs"))
Zkj <- qs::qread(file.path(EXP, "cache/Z_karjalainen.qs"))
# accept matrix or list
if (is.list(Zuk) && !is.matrix(Zuk) && !is.null(Zuk$Z)) Zuk <- Zuk$Z
if (is.list(Zkj) && !is.matrix(Zkj) && !is.null(Zkj$Z)) Zkj <- Zkj$Z
Zuk <- as.matrix(Zuk); Zkj <- as.matrix(Zkj)
logmsg("Zuk dim=", paste(dim(Zuk), collapse="x"),
       " Zkj dim=", paste(dim(Zkj), collapse="x"),
       " col_match_traits=", isTRUE(all(traits %in% colnames(Zuk))))
if (!all(traits %in% colnames(Zuk))) {
  logmsg("Zuk cols sample=", paste(colnames(Zuk)[1:min(5, ncol(Zuk))], collapse=","))
  stop("Z cache columns do not contain the 61 traits")
}
Zuk <- Zuk[, traits, drop = FALSE]
if (all(traits %in% colnames(Zkj))) Zkj <- Zkj[, traits, drop = FALSE]

snp_z <- rownames(Zuk)
if (is.null(snp_z)) stop("Zuk has no rownames")

UK_FG <- file.path(EXP, "gls_Rreg/fgwas_ukbb_eur_gls_Rreg_genomewide.tsv.gz")
KJ_FG <- file.path(EXP, "gls_Rreg/fgwas_karjalainen_gls_Rreg_genomewide.tsv.gz")
zcols <- paste0("Z_", RETAIN)
pcols <- paste0("P_", RETAIN)
logmsg("Reading UKB factor GWAS (selected columns)")
fg_uk <- fread(UK_FG, select = c("SNP", zcols, pcols))
logmsg("Reading Karj factor GWAS (selected columns)")
fg_kj <- fread(KJ_FG, select = c("SNP", zcols, pcols))
setnames(fg_uk, zcols, paste0("UKB_", RETAIN))
setnames(fg_uk, pcols, paste0("UKB_P_", RETAIN))
setnames(fg_kj, zcols, paste0("KARJ_", RETAIN))
setnames(fg_kj, pcols, paste0("KARJ_P_", RETAIN))

# Align SNPs: intersection of factor GWAS, Z, map
common <- Reduce(intersect, list(fg_uk$SNP, rownames(Zuk), smap$SNP))
logmsg("common SNPs factorGWAS x Zuk x map = ", length(common))
fg_uk <- fg_uk[SNP %in% common]
smap_u <- unique(smap[SNP %in% common, .(SNP, CHR, BP)], by = "SNP")
fg_uk <- merge(fg_uk, smap_u, by = "SNP", all.x = TRUE)
# keep Z row order index
z_index <- setNames(seq_len(nrow(Zuk)), rownames(Zuk))

# GLS weights from official formula (plus_resid=FALSE): try loadings_gls source + reconstruct
Rpath <- file.path(EXP, "R/R_ldsc_cond1e3_ukbb_eur.tsv")
if (!file.exists(Rpath)) {
  # find
  cands <- list.files(file.path(EXP, "R"), pattern = "R_ldsc_cond1e3_ukbb", full.names = TRUE)
  logmsg("R candidates: ", paste(cands, collapse = " | "))
  if (length(cands)) Rpath <- cands[1]
}
W <- NULL
ldo_ok <- FALSE
gls_note <- ""
if (file.exists(Rpath)) {
  Ruk <- as.matrix(fread(Rpath), rownames = 1)
  # fread rownames via first col
  if (is.null(rownames(Ruk)) || !all(traits %in% rownames(Ruk))) {
    Rdt <- fread(Rpath)
    rn <- Rdt[[1]]
    Ruk <- as.matrix(Rdt[, -1, with = FALSE])
    rownames(Ruk) <- rn
    colnames(Ruk) <- names(Rdt)[-1]
  }
  Ruk <- Ruk[traits, traits, drop = FALSE]
  ev <- eigen(Ruk, symmetric = TRUE, only.values = TRUE)$values
  logmsg("Ruk dim=", paste(dim(Ruk), collapse="x"), " min_ev=", min(ev))
  # dump loadings_gls
  gls_src <- tryCatch(paste(deparse(GFA:::loadings_gls), collapse = "\n"), error = function(e) e$message)
  writeLines(c("GFA:::loadings_gls", gls_src), file.path(OUT, "P2C_LOADINGS_GLS_SOURCE.txt"))
  stage(file.path(OUT, "P2C_LOADINGS_GLS_SOURCE.txt"))
  # Official project used plus_resid=FALSE, Sigma=R, S=1
  # Standard GLS: score_j = t(F_j) %*% solve(R) %*% z  / sqrt( t(F_j) %*% solve(R) %*% F_j )
  Ri <- tryCatch(solve(Ruk), error = function(e) NULL)
  if (!is.null(Ri)) {
    W <- Ri %*% Fgls[, RETAIN, drop = FALSE]   # 61 x 13, un-normalized
    den <- sqrt(colSums(Fgls[, RETAIN, drop = FALSE] * W))
    Wn <- sweep(W, 2, den, "/")
    colnames(Wn) <- RETAIN
    rownames(Wn) <- traits
    W <- Wn
    ldo_ok <- TRUE
    gls_note <- "Weights are official-style GLS: w = R^{-1} F / sqrt(F' R^{-1} F), plus_resid=FALSE, F=F_hat*scale, R=UKB R_ldsc_cond1e3. Confirm against GFA:::loadings_gls source dumped in P2C_LOADINGS_GLS_SOURCE.txt."
    # QC: reconstruct Z_f approx and correlate with official factor Z on a subset
    # Use only SNPs present in both
    common2 <- intersect(fg_uk$SNP, names(z_index))
    # sample or all — 1.2e6 * 61 may be heavy but we already have Z
    logmsg("QC reconstructing GLS scores on ", length(common2), " SNPs")
    Zc <- Zuk[z_index[common2], traits, drop = FALSE]
    Zc[!is.finite(Zc)] <- 0
    Zrec <- Zc %*% W
    qc_cor <- sapply(RETAIN, function(fac) {
      zf <- fg_uk[[paste0("UKB_", fac)]][match(common2, fg_uk$SNP)]
      zr <- Zrec[, fac]
      ok <- is.finite(zf) & is.finite(zr)
      if (sum(ok) < 1000) return(NA_real_)
      cor(zf[ok], zr[ok])
    })
    logmsg("GLS reconstruct vs official Z cor: ", paste(sprintf("%s=%.6f", RETAIN, qc_cor), collapse = ", "))
    fwrite(data.table(factor = RETAIN, cor_reconstructed_vs_official_Z = as.numeric(qc_cor),
                      note = gls_note),
           file.path(OUT, "P2C_GLS_WEIGHT_QC.tsv"), sep = "\t")
    stage(file.path(OUT, "P2C_GLS_WEIGHT_QC.tsv"))
    if (any(!is.finite(qc_cor) | abs(qc_cor) < 0.999)) {
      logmsg("WARNING: reconstructed GLS Z does not match official factor Z at |r|>=0.999; LDO will be flagged PARTIAL")
      ldo_ok <- all(is.finite(qc_cor) & abs(qc_cor) >= 0.99)
    }
    rm(Zc, Zrec); gc(verbose = FALSE)
  } else {
    gls_note <- "solve(R) failed; leave-dominant-out GLS blocked"
  }
} else {
  gls_note <- "UKB residual R not found; leave-dominant-out GLS blocked"
  logmsg(gls_note)
}

# ---- 3. Factor vs dominant trait genome-wide + leads ----
# Z dominant from Zuk; factor Z from fg_uk
fg_kj_u <- fg_kj[SNP %in% fg_uk$SNP]
setkey(fg_kj_u, SNP)

cmp_rows <- list()
lead_rows <- list()   # factor p5e8 leads with trait-locus audit
cond_rows <- list()

# Pre-index map by chr for window searches
smap_u[, row_z := z_index[SNP]]
smap_u <- smap_u[is.finite(row_z)]
setkey(smap_u, CHR, BP)

# Precompute per-trait genome-wide significant row indices (UKB)
trait_sig_index <- lapply(traits, function(tr) {
  z <- Zuk[, tr]
  which(is.finite(z) & abs(z) >= Z_GWAS)
})
names(trait_sig_index) <- traits
trait_shadow_index <- lapply(traits, function(tr) {
  z <- Zuk[, tr]
  which(is.finite(z) & abs(z) >= SHADOW_Z)
})
names(trait_shadow_index) <- traits

for (fac in RETAIN) {
  logmsg("Factor ", fac)
  dom <- base0[factor == fac, dominant_trait]
  zf <- fg_uk[[paste0("UKB_", fac)]]
  pf <- fg_uk[[paste0("UKB_P_", fac)]]
  if (all(is.na(pf))) pf <- p_from_z(zf)
  # dominant trait Z on same SNP order as fg_uk
  ix <- z_index[fg_uk$SNP]
  zd <- Zuk[ix, dom]
  zd[!is.finite(zd)] <- NA_real_
  ok <- is.finite(zf) & is.finite(zd)
  pear <- cor(zf[ok], zd[ok], method = "pearson")
  spear <- cor(zf[ok], zd[ok], method = "spearman")
  # sign at |Z_f|>=8 and |Z_dom|>=8
  w8f <- ok & abs(zf) >= 8
  w8d <- ok & abs(zd) >= 8
  sign8f <- if (sum(w8f) >= 10) mean(sign(zf[w8f]) == sign(zd[w8f])) else NA_real_
  sign8d <- if (sum(w8d) >= 10) mean(sign(zf[w8d]) == sign(zd[w8d])) else NA_real_
  # OLS residual Zf ~ Zd
  fit_ols <- lm(zf[ok] ~ zd[ok])
  b <- coef(fit_ols)
  r2 <- summary(fit_ols)$r.squared
  z_resid <- rep(NA_real_, length(zf))
  z_resid[ok] <- as.numeric(residuals(fit_ols))
  # LDO weights
  z_ldo <- rep(NA_real_, length(zf))
  ldo_status <- "BLOCKED"
  if (isTRUE(ldo_ok) && !is.null(W) && fac %in% colnames(W)) {
    w <- W[, fac]
    w_ldo <- w
    w_ldo[dom] <- 0
    # renormalize in the GLS metric: denom = sqrt(w_ldo' R w_ldo) = sqrt(F_ldo' R^{-1} F_ldo)
    # w already equals R^{-1} F / s, so F_ldo = F with F[dom]=0, recompute
    Fldo <- Fgls[, fac]
    Fldo[dom] <- 0
    if (file.exists(Rpath)) {
      Ri_dom_den <- as.numeric(sqrt(t(Fldo) %*% solve(Ruk, Fldo)))
      if (is.finite(Ri_dom_den) && Ri_dom_den > 0) {
        w_ldo2 <- as.numeric(solve(Ruk, Fldo) / Ri_dom_den)
        names(w_ldo2) <- traits
        Zc <- Zuk[ix, traits, drop = FALSE]
        Zc[!is.finite(Zc)] <- 0
        z_ldo <- as.numeric(Zc %*% w_ldo2)
        ldo_status <- "COMPLETE"
        rm(Zc)
      }
    }
  }

  # clump factor and dominant at p<5e-8
  cl_f <- clump_distance(fg_uk$CHR, fg_uk$BP, zf, kb = CLUMP_KB_MAIN) & ok & (p_from_z(zf) < P_GWAS)
  cl_d <- clump_distance(fg_uk$CHR, fg_uk$BP, zd, kb = CLUMP_KB_MAIN) & ok & (p_from_z(zd) < P_GWAS)
  n_f <- sum(cl_f); n_d <- sum(cl_d)
  # overlap: factor lead has a dominant lead within 1 Mb
  f_leads <- fg_uk[cl_f, .(SNP, CHR, BP, zf = zf[cl_f], pf = p_from_z(zf[cl_f]))]
  d_leads <- data.table(SNP = fg_uk$SNP[cl_d], CHR = fg_uk$CHR[cl_d], BP = fg_uk$BP[cl_d],
                        zd = zd[cl_d], pd = p_from_z(zd[cl_d]))
  n_overlap_f <- 0L
  if (n_f && n_d) {
    f_leads[, overlap_dom := FALSE]
    for (i in seq_len(nrow(f_leads))) {
      hit <- any(d_leads$CHR == f_leads$CHR[i] & abs(d_leads$BP - f_leads$BP[i]) <= CLUMP_KB_MAIN * 1000)
      f_leads$overlap_dom[i] <- hit
    }
    n_overlap_f <- sum(f_leads$overlap_dom)
  } else if (n_f) f_leads[, overlap_dom := FALSE]

  n_resid_sig <- sum(cl_f & is.finite(z_resid) & p_from_z(z_resid) < P_GWAS, na.rm = TRUE)
  n_ldo_sig <- sum(cl_f & is.finite(z_ldo) & p_from_z(z_ldo) < P_GWAS, na.rm = TRUE)

  cmp_rows[[fac]] <- cbind(
    base0[factor == fac],
    data.table(
      n_common_snps_used = sum(ok),
      pearson_Z_factor_vs_dominant = pear,
      spearman_Z_factor_vs_dominant = spear,
      sign_concord_factor_zge8 = sign8f,
      n_snps_factor_zge8 = sum(w8f),
      sign_concord_dominant_zge8 = sign8d,
      n_snps_dominant_zge8 = sum(w8d),
      ols_slope_Zf_on_Zdom = unname(b[2]),
      ols_intercept = unname(b[1]),
      ols_r2 = r2,
      n_factor_leads_p5e8_1Mb = n_f,
      n_dominant_leads_p5e8_1Mb = n_d,
      n_factor_leads_overlapping_dominant_lead_1Mb = n_overlap_f,
      frac_factor_leads_overlapping_dominant_lead = if (n_f) n_overlap_f / n_f else NA_real_,
      n_factor_leads_with_OLS_residual_p5e8 = n_resid_sig,
      n_factor_leads_with_LDO_p5e8 = if (ldo_status == "COMPLETE") n_ldo_sig else NA_integer_,
      ldo_status = ldo_status,
      ldsc_rg_factor_vs_dominant = NA_real_,
      ldsc_rg_status = "NOT_RUN",
      ldsc_rg_note = "Bivariate LDSC rg of factor vs dominant trait was not pre-computed and was not launched in this package (concurrent P2A/LAVA load). Genome-wide Z correlation and 1 Mb lead overlap are the primary comparisons."
    )
  )

  # ---- 4. Per-lead audit vs all 61 traits ----
  if (!n_f) next
  # also mergeto karj factor Z
  f_leads <- merge(f_leads, fg_kj_u[, c("SNP", paste0("KARJ_", fac), paste0("KARJ_P_", fac)), with = FALSE],
                   by = "SNP", all.x = TRUE)
  setnames(f_leads, paste0("KARJ_", fac), "Z_factor_karj", skip_absent = TRUE)
  setnames(f_leads, paste0("KARJ_P_", fac), "P_factor_karj", skip_absent = TRUE)
  if (!"Z_factor_karj" %in% names(f_leads)) f_leads[, Z_factor_karj := NA_real_]
  if (!"P_factor_karj" %in% names(f_leads)) f_leads[, P_factor_karj := NA_real_]

  for (i in seq_len(nrow(f_leads))) {
    chr <- f_leads$CHR[i]; bp <- f_leads$BP[i]; snp <- f_leads$SNP[i]
    win <- smap_u[CHR == chr & BP >= (bp - CLUMP_KB_MAIN * 1000) & BP <= (bp + CLUMP_KB_MAIN * 1000)]
    rows_w <- win$row_z
    # max |Z| per trait in window
    if (length(rows_w)) {
      Zw <- Zuk[rows_w, traits, drop = FALSE]
      maxabs <- apply(abs(Zw), 2, function(x) if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE))
    } else {
      maxabs <- setNames(rep(NA_real_, length(traits)), traits)
    }
    n_traits_gwas_sig_locus <- sum(is.finite(maxabs) & maxabs >= Z_GWAS)
    n_traits_shadow_locus <- sum(is.finite(maxabs) & maxabs >= SHADOW_Z)
    traits_sig <- names(maxabs)[is.finite(maxabs) & maxabs >= Z_GWAS]
    # values at exact lead SNP
    i_lead <- z_index[[snp]]
    z_lead_all <- if (!is.na(i_lead)) Zuk[i_lead, traits] else setNames(rep(NA_real_, length(traits)), traits)
    z_dom_lead <- as.numeric(z_lead_all[[dom]])
    z_dom_max_locus <- as.numeric(maxabs[[dom]])
    max_trait <- if (any(is.finite(maxabs))) names(which.max(replace(maxabs, !is.finite(maxabs), -Inf))) else NA_character_
    max_trait_z <- if (!is.na(max_trait)) as.numeric(maxabs[[max_trait]]) else NA_real_
    # contributing at lead: |w_j Z_j|
    contrib <- NULL
    if (!is.null(W) && fac %in% colnames(W) && !is.na(i_lead)) {
      contrib <- abs(W[, fac] * replace(z_lead_all, !is.finite(z_lead_all), 0))
    }
    top_contrib <- if (!is.null(contrib)) {
      o <- order(contrib, decreasing = TRUE)
      paste(sprintf("%s:%.3f", traits[o[1:min(5, length(o))]], contrib[o[1:min(5, length(o))]]), collapse = ";")
    } else ""

    already <- n_traits_gwas_sig_locus >= 1
    shadow <- (!already) && (is.finite(z_dom_lead) && abs(z_dom_lead) >= SHADOW_Z ||
                               is.finite(z_dom_max_locus) && abs(z_dom_max_locus) >= SHADOW_Z)
    if (already) {
      cls <- "already_discovered"
    } else if (shadow) {
      cls <- "significance_gain_dominant_shadow"
    } else {
      cls <- "candidate_incremental_threshold"
    }
    cls_note <- "All classes remain inside the linear span of the 61 trait Z scores at each SNP (GLS factor Z is a linear combination). candidate_incremental_threshold is a significance-threshold discovery, not a new independent causal variant."

    # Karj same SNP for factor and dominant
    z_dom_kj <- if (!is.null(colnames(Zkj)) && dom %in% colnames(Zkj) && !is.na(i_lead) &&
                    !is.null(rownames(Zkj)) && snp %in% rownames(Zkj)) {
      as.numeric(Zkj[snp, dom])
    } else if (!is.null(colnames(Zkj)) && dom %in% colnames(Zkj) && !is.na(i_lead) &&
               !is.null(rownames(Zkj)) && identical(rownames(Zkj), rownames(Zuk))) {
      as.numeric(Zkj[i_lead, dom])
    } else NA_real_

    zf_kj <- f_leads$Z_factor_karj[i]
    pf_kj <- f_leads$P_factor_karj[i]
    karj_flag <- "ok"
    if (!is.finite(zf_kj)) karj_flag <- "missing_or_unaligned_factor_Z_karj"

    lead_rows[[length(lead_rows) + 1L]] <- data.table(
      factor = fac,
      dominant_trait = dom,
      lead_SNP = snp,
      CHR = chr,
      BP = bp,
      in_mhc = in_mhc(chr, bp),
      clump_kb = CLUMP_KB_MAIN,
      p_threshold = P_GWAS,
      Z_factor_ukb = f_leads$zf[i],
      P_factor_ukb = f_leads$pf[i],
      Z_dominant_ukb_at_lead = z_dom_lead,
      P_dominant_ukb_at_lead = p_from_z(z_dom_lead),
      Z_dominant_ukb_max_in_1Mb = z_dom_max_locus,
      n_traits_p5e8_in_1Mb = n_traits_gwas_sig_locus,
      n_traits_absZ_ge4_in_1Mb = n_traits_shadow_locus,
      traits_p5e8_in_1Mb = paste(traits_sig, collapse = ";"),
      trait_with_max_absZ_in_1Mb = max_trait,
      max_absZ_any_trait_in_1Mb = max_trait_z,
      overlaps_dominant_trait_lead = if ("overlap_dom" %in% names(f_leads)) isTRUE(f_leads$overlap_dom[i]) else NA,
      top_GLS_contributions_at_lead = top_contrib,
      Z_factor_OLS_residual_after_dominant = z_resid[match(snp, fg_uk$SNP)],
      P_factor_OLS_residual_after_dominant = p_from_z(z_resid[match(snp, fg_uk$SNP)]),
      Z_factor_LDO = z_ldo[match(snp, fg_uk$SNP)],
      P_factor_LDO = p_from_z(z_ldo[match(snp, fg_uk$SNP)]),
      ldo_status = ldo_status,
      Z_factor_karj = zf_kj,
      P_factor_karj = pf_kj,
      sign_concord_factor_ukb_karj = if (is.finite(zf_kj) && is.finite(f_leads$zf[i])) as.integer(sign(zf_kj) == sign(f_leads$zf[i])) else NA_integer_,
      Z_dominant_karj_at_lead = z_dom_kj,
      P_dominant_karj_at_lead = p_from_z(z_dom_kj),
      karj_qc_flag = karj_flag,
      locus_class = cls,
      class_note = cls_note
    )
  }
  gc(verbose = FALSE)
}

cmp <- rbindlist(cmp_rows, fill = TRUE)
leads <- rbindlist(lead_rows, fill = TRUE)

# Summarize novel locus audit at factor level
audit_sum <- leads[, .(
  n_factor_leads_p5e8 = .N,
  n_already_discovered = sum(locus_class == "already_discovered"),
  n_significance_gain_dominant_shadow = sum(locus_class == "significance_gain_dominant_shadow"),
  n_candidate_incremental_threshold = sum(locus_class == "candidate_incremental_threshold"),
  n_mhc = sum(in_mhc),
  n_already_discovered_nonMHC = sum(locus_class == "already_discovered" & !in_mhc),
  n_candidate_incremental_nonMHC = sum(locus_class == "candidate_incremental_threshold" & !in_mhc),
  n_candidate_with_karj_factor_p5e8 = sum(locus_class == "candidate_incremental_threshold" & is.finite(P_factor_karj) & P_factor_karj < P_GWAS),
  n_candidate_sign_concordant_karj = sum(locus_class == "candidate_incremental_threshold" & sign_concord_factor_ukb_karj == 1, na.rm = TRUE)
), by = factor]

fwrite(cmp, file.path(OUT, "FACTOR_VS_TRAIT_BASELINE.tsv"), sep = "\t")
fwrite(leads, file.path(OUT, "FACTOR_NOVEL_LOCUS_CATALOG.tsv"), sep = "\t")
fwrite(audit_sum, file.path(OUT, "FACTOR_NOVEL_LOCUS_AUDIT.tsv"), sep = "\t")
# candidate-only subset for easier review
fwrite(leads[locus_class == "candidate_incremental_threshold"],
       file.path(OUT, "FACTOR_NOVEL_LOCUS_CANDIDATES.tsv"), sep = "\t")

for (fn in c("FACTOR_VS_TRAIT_BASELINE.tsv", "FACTOR_NOVEL_LOCUS_CATALOG.tsv",
             "FACTOR_NOVEL_LOCUS_AUDIT.tsv", "FACTOR_NOVEL_LOCUS_CANDIDATES.tsv")) {
  stage(file.path(OUT, fn))
}

logmsg("LEAD CLASS TOTALS already=", sum(leads$locus_class=="already_discovered"),
       " shadow=", sum(leads$locus_class=="significance_gain_dominant_shadow"),
       " candidate=", sum(leads$locus_class=="candidate_incremental_threshold"),
       " n_leads=", nrow(leads))
print(audit_sum)
logmsg("DONE P2C")
