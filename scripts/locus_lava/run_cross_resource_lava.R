#!/usr/bin/env Rscript
# Same-factor LAVA local rg: UKB vs Karjalainen, 13 main retained factors only.
# Requires complete UKB EUR v1.1 LD reference for chr1-22.
# Writes only under results/exp52_9/local_genetic/lava/.
# Does not run cross-factor or factor-disease LAVA.

options(stringsAsFactors = FALSE)
.libPaths(c(Sys.getenv("LAVA_R_LIBRARY", unset = file.path("external", "lava", "R_library_0.1.5")), .libPaths()))
suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
OUT  <- file.path(PROJ, "results/exp52_9/local_genetic")
LAVA_DIR <- file.path(OUT, "lava")
SUM_DIR <- file.path(OUT, "lava_sumstats")
LOC_FILE <- Sys.getenv("LAVA_LOCUS_FILE", unset = file.path(PROJ, "external/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
REF_PREFIX <- Sys.getenv("LAVA_UKB_REFERENCE_PREFIX", unset = file.path(PROJ, "external/lava/reference_ukb_v1.1/selected/lava-ukb-v1.1"))
INFO_IN <- file.path(OUT, "lava_input.info.txt")
UNIV_THRESH <- 0.05
MHC_START <- 25000000
MHC_STOP <- 34000000
MAIN_FACTORS <- c("F9","F11","F13","F14","F15","F17","F18","F20","F22","F23","F24","F25","F26")
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
       ". Do not launch LAVA until chr1-22 are present. Do not restart existing curl/aria2 downloads from this script.")
}
if (!file.exists(LOC_FILE)) stop("Missing locus file: ", LOC_FILE)
if (!file.exists(INFO_IN)) stop("Missing lava_input.info.txt; run run_exp52_9_local_genetic.R first")

loci_all <- read.loci(LOC_FILE)
loci_all <- loci_all[loci_all$CHR %in% 1:22, ]
loci_all$mhc <- loci_all$CHR == 6 & loci_all$START < MHC_STOP & loci_all$STOP > MHC_START
logmsg("loci=", nrow(loci_all), " MHC loci=", sum(loci_all$mhc))

info_all <- fread(INFO_IN)
safe_df <- function(x) {
  if (is.null(x)) return(NULL)
  as.data.frame(x, stringsAsFactors = FALSE)
}
append_table <- function(x, file) {
  if (is.null(x) || nrow(as.data.table(x)) == 0) return(invisible(NULL))
  fwrite(as.data.table(x), file, sep = "\t", append = file.exists(file), col.names = !file.exists(file))
}

univ_file <- file.path(LAVA_DIR, "lava_univ_ukb_karj_main.tsv")
bivar_file <- file.path(LAVA_DIR, "lava_bivar_ukb_karj_main.tsv")
status_file <- file.path(LAVA_DIR, "lava_locus_status.tsv")
for (f in c(univ_file, bivar_file, status_file)) if (file.exists(f)) file.remove(f)

for (fac in MAIN_FACTORS) {
  ph_u <- paste0("UKB_", fac)
  ph_k <- paste0("KARJ_", fac)
  info <- info_all[phenotype %in% c(ph_u, ph_k)]
  if (nrow(info) != 2L) stop("input.info missing pair for ", fac)
  if (any(!file.exists(info$filename))) stop("Missing sumstats for ", fac, ": ", paste(info$filename[!file.exists(info$filename)], collapse = ", "))

  work <- file.path(LAVA_DIR, "work", fac)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(work, "input.info.txt")
  ov_path <- file.path(work, "sample.overlap.txt")
  # Quantitative GWAS: write explicit NA. LAVA get.input.info() uses read.table()
  # with default whitespace, so empty cases/controls collapse and fail.
  info[, cases := NA_real_]
  info[, controls := NA_real_]
  write.table(as.data.frame(info[, .(phenotype, cases, controls, filename)]),
              file = info_path, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
  ov <- diag(2)
  dimnames(ov) <- list(c(ph_u, ph_k), c(ph_u, ph_k))
  write.table(ov, ov_path, quote = FALSE, sep = "\t", col.names = NA)

  logmsg("process.input ", fac)
  input <- process.input(
    input.info.file = info_path,
    sample.overlap.file = ov_path,
    ref.prefix = REF_PREFIX,
    phenos = c(ph_u, ph_k),
    input.dir = NULL
  )

  for (i in seq_len(nrow(loci_all))) {
    loc_meta <- loci_all[i, , drop = FALSE]
    loc_id <- loc_meta$LOC
    chr <- loc_meta$CHR
    locus <- tryCatch(
      process.locus(loc_meta, input, phenos = c(ph_u, ph_k)),
      error = function(e) {
        append_table(data.table(factor = fac, chr = chr, locus = loc_id,
                                start = loc_meta$START, stop = loc_meta$STOP, mhc = loc_meta$mhc,
                                status = "process_locus_error", detail = conditionMessage(e)), status_file)
        NULL
      }
    )
    if (is.null(locus)) next
    univ <- tryCatch(safe_df(run.univ(locus, phenos = c(ph_u, ph_k))), error = function(e) {
      append_table(data.table(factor = fac, chr = chr, locus = loc_id,
                              start = loc_meta$START, stop = loc_meta$STOP, mhc = loc_meta$mhc,
                              status = "run_univ_error", detail = conditionMessage(e)), status_file)
      NULL
    })
    if (!is.null(univ) && nrow(univ) > 0) {
      univ$factor <- fac
      univ$chr <- chr
      univ$locus <- loc_id
      univ$start <- loc_meta$START
      univ$stop <- loc_meta$STOP
      univ$mhc <- loc_meta$mhc
      append_table(univ, univ_file)
    }
    if (is.null(univ)) next
    phen_col <- if ("phen" %in% names(univ)) "phen" else if ("phenotype" %in% names(univ)) "phenotype" else NA
    if (is.na(phen_col)) next
    p_u <- univ[["p"]][univ[[phen_col]] == ph_u]
    p_k <- univ[["p"]][univ[[phen_col]] == ph_k]
    if (length(p_u) != 1 || length(p_k) != 1) next
    if (!is.finite(p_u) || !is.finite(p_k) || p_u >= UNIV_THRESH || p_k >= UNIV_THRESH) next
    b <- tryCatch(safe_df(run.bivar(locus, phenos = c(ph_u, ph_k))), error = function(e) {
      append_table(data.table(factor = fac, chr = chr, locus = loc_id,
                              status = "run_bivar_error", detail = conditionMessage(e)), status_file)
      NULL
    })
    if (is.null(b) || nrow(b) == 0) next
    b$factor <- fac
    b$chr <- chr
    b$locus <- loc_id
    b$start <- loc_meta$START
    b$stop <- loc_meta$STOP
    b$mhc <- loc_meta$mhc
    b$univ_p_ukb <- p_u
    b$univ_p_karj <- p_k
    append_table(b, bivar_file)
  }
  logmsg("finished loci for ", fac)
}

# FDR: primary = within-factor, non-MHC bivariate tests; secondary = all non-MHC tests pooled.
if (file.exists(bivar_file)) {
  bv <- fread(bivar_file)
  pcol <- if ("p" %in% names(bv)) "p" else if ("pval" %in% names(bv)) "pval" else NA
  if (!is.na(pcol)) {
    bv[, fdr_within_factor_nonMHC := NA_real_]
    bv[mhc == FALSE, fdr_within_factor_nonMHC := p.adjust(get(pcol), method = "BH"), by = factor]
    nonmhc <- bv$mhc == FALSE
    bv[, fdr_pooled_nonMHC := NA_real_]
    if (any(nonmhc)) bv[nonmhc, fdr_pooled_nonMHC := p.adjust(get(pcol), method = "BH")]
    fwrite(bv, file.path(LAVA_DIR, "lava_bivar_ukb_karj_main_fdr.tsv"), sep = "\t")
    summ <- bv[, .(
      n_bivar = .N,
      n_bivar_nonMHC = sum(mhc == FALSE),
      n_fdr05_within_nonMHC = sum(mhc == FALSE & is.finite(fdr_within_factor_nonMHC) & fdr_within_factor_nonMHC < 0.05),
      n_fdr05_pooled_nonMHC = sum(mhc == FALSE & is.finite(fdr_pooled_nonMHC) & fdr_pooled_nonMHC < 0.05),
      median_rg_nonMHC = median(if ("rho" %in% names(bv)) rho[mhc == FALSE] else if ("rg" %in% names(bv)) rg[mhc == FALSE] else NA_real_, na.rm = TRUE)
    ), by = factor]
    fwrite(summ, file.path(LAVA_DIR, "lava_bivar_summary_main.tsv"), sep = "\t")
  }
}

writeLines(c(paste0("done=", now_str(), " ", TZ_LABEL), paste0("univ.thresh=", UNIV_THRESH)), file.path(LAVA_DIR, "LAVA_DONE.txt"))
logmsg("ALL LAVA DONE")
