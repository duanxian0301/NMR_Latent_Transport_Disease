#!/usr/bin/env Rscript
# Export complete 61x30 loadings from seed-101 RDS (and seed-202 comparison).
# Does not overwrite GFA RDS/TSV, Z, R, LAVA, or manuscript files.
# Writes only under audit_followup_20261117/P2_FINAL/.

suppressPackageStartupMessages({
  library(data.table)
})

PROJ <- normalizePath(Sys.getenv("NMR_GFA_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)
GFA  <- file.path(PROJ, "results/exp52_9/gfa")
OUT  <- file.path(PROJ, "audit_followup_20261117/P2_FINAL")
STAGE <- Sys.getenv("NMR_GFA_AUDIT_MIRROR", unset = OUT)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "scripts"), showWarnings = FALSE, recursive = TRUE)
dir.create(STAGE, showWarnings = FALSE, recursive = TRUE)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
      paste0(..., collapse = ""), "\n", sep = "")
  flush.console()
}

copy_also_stage <- function(path) {
  if (file.exists(path) && dir.exists(STAGE)) {
    file.copy(path, file.path(STAGE, basename(path)), overwrite = TRUE)
  }
}

write_loading_matrix <- function(M, path) {
  dt <- data.table(trait = rownames(M))
  dt <- cbind(dt, as.data.table(M))
  fwrite(dt, path, sep = "\t")
}

count_zeros <- function(M, tol = 0) {
  v <- as.numeric(M)
  list(
    n = length(v),
    n_exact_zero = sum(v == 0),
    n_abs_le_1e16 = sum(abs(v) <= 1e-16),
    n_abs_le_1e12 = sum(abs(v) <= 1e-12),
    n_abs_le_1e8  = sum(abs(v) <= 1e-8),
    n_na = sum(!is.finite(v)),
    min_abs_nonzero = if (any(v != 0 & is.finite(v))) min(abs(v[v != 0 & is.finite(v)])) else NA_real_,
    min = min(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE)
  )
}

RDS_101 <- file.path(GFA, "ukbb_ldclump_seed101.rds")
RDS_202 <- file.path(GFA, "ukbb_ldclump_seed202.rds")
TSV_101 <- file.path(GFA, "ukbb_ldclump_seed101_F_hat.tsv")
TRAIT_N <- file.path(PROJ, "results/exp52_9/trait_N.tsv")

logmsg("START export loadings and seed-identity inspect")
logmsg("Reading seed 101 RDS: ", RDS_101, " size=", file.info(RDS_101)$size)
fit101 <- readRDS(RDS_101)
logmsg("seed101 class=", paste(class(fit101), collapse = ","),
       " names=", paste(names(fit101), collapse = ", "))

if (is.null(fit101$F_hat)) stop("fit101 missing F_hat")
if (is.null(fit101$scale)) stop("fit101 missing scale")

Fh <- as.matrix(fit101$F_hat)
sc <- fit101$scale
logmsg("F_hat dim=", paste(dim(Fh), collapse = "x"),
       " rownames_n=", length(rownames(Fh)),
       " colnames=", paste(colnames(Fh), collapse = ","),
       " scale_length=", length(sc),
       " scale_names=", if (is.null(names(sc))) "NULL" else paste(names(sc), collapse = ","))

# Align scale to rows of F_hat (per-trait)
if (!is.null(names(sc)) && !is.null(rownames(Fh)) && all(rownames(Fh) %in% names(sc))) {
  sc <- sc[rownames(Fh)]
} else if (length(sc) != nrow(Fh)) {
  stop("Cannot align scale (len=", length(sc), ") to F_hat rows (", nrow(Fh), ")")
}
Fgls <- Fh * as.numeric(sc)  # traits x factors; scale recycled down rows in R if len==nrow
rownames(Fgls) <- rownames(Fh)
colnames(Fgls) <- colnames(Fh)
if (is.null(colnames(Fgls)) || !grepl("^F", colnames(Fgls)[1])) {
  colnames(Fh) <- colnames(Fgls) <- paste0("F", seq_len(ncol(Fh)))
}

# Optional trait_N order check (do not reorder unless needed; keep RDS row order)
if (file.exists(TRAIT_N)) {
  tn <- fread(TRAIT_N)
  logmsg("trait_N traits=", nrow(tn),
         " identical_order_to_F=", identical(tn$trait_code, rownames(Fh)),
         " setequal=", isTRUE(all(sort(tn$trait_code) == sort(rownames(Fh)))))
}

z_fh <- count_zeros(Fh)
z_gls <- count_zeros(Fgls)
logmsg("F_hat exact zeros=", z_fh$n_exact_zero, "/", z_fh$n,
       " min_abs_nonzero=", format(z_fh$min_abs_nonzero, digits = 8))
logmsg("F_hat*scale exact zeros=", z_gls$n_exact_zero, "/", z_gls$n,
       " min_abs_nonzero=", format(z_gls$min_abs_nonzero, digits = 8))

# Compare to exported TSV (read-only)
if (file.exists(TSV_101)) {
  Ft <- as.matrix(read.delim(TSV_101, check.names = FALSE, row.names = 1))
  logmsg("exported TSV dim=", paste(dim(Ft), collapse = "x"),
         " exact zeros=", count_zeros(Ft)$n_exact_zero,
         " max_abs_diff_vs_RDS_Fhat=",
         if (identical(dim(Ft), dim(Fh))) max(abs(Ft - Fh)) else "DIM_MISMATCH")
}

path_fhat <- file.path(OUT, "FULL_GFA_LOADINGS_61x30_Fhat.tsv")
path_gls  <- file.path(OUT, "FULL_GFA_LOADINGS_61x30.tsv")
write_loading_matrix(Fh, path_fhat)
write_loading_matrix(Fgls, path_gls)

# scale vector
tr_sc <- if (!is.null(names(sc))) names(sc) else rownames(Fh)
sc_dt <- data.table(trait = tr_sc, scale = as.numeric(sc))
fwrite(sc_dt, file.path(OUT, "FULL_GFA_SCALE_VECTOR.tsv"), sep = "\t")

readme <- c(
  "# Full GFA 61x30 loading matrices (seed 101)",
  "",
  "Source: `results/exp52_9/gfa/ukbb_ldclump_seed101.rds` (not the pre-exported TSV, not S4a filtered loadings).",
  paste0("Exported: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  "",
  "## Files",
  "",
  "- `FULL_GFA_LOADINGS_61x30.tsv` — **this is `F_hat * scale`**, the matrix that entered official GLS (`GFA:::loadings_gls`). Use this for dominance / concentration / factor-vs-trait work unless a note says otherwise.",
  "- `FULL_GFA_LOADINGS_61x30_Fhat.tsv` — raw `fit$F_hat` from the RDS.",
  "- `FULL_GFA_SCALE_VECTOR.tsv` — the length-61 per-trait `fit$scale` used in the product.",
  "",
  "## Zero counts (exact `== 0`)",
  "",
  paste0("- F_hat: ", z_fh$n_exact_zero, " / ", z_fh$n,
         "; min |nonzero| = ", format(z_fh$min_abs_nonzero, digits = 8)),
  paste0("- F_hat * scale: ", z_gls$n_exact_zero, " / ", z_gls$n,
         "; min |nonzero| = ", format(z_gls$min_abs_nonzero, digits = 8)),
  paste0("- Additional near-zero counts for F_hat: |x|<=1e-16 ", z_fh$n_abs_le_1e16,
         "; |x|<=1e-12 ", z_fh$n_abs_le_1e12,
         "; |x|<=1e-8 ", z_fh$n_abs_le_1e8),
  "",
  if (z_fh$n_exact_zero == 0)
    "The RDS loading matrix is dense: there are **no exact structural zeros**. A complete 61x30 export is therefore dense. S4a / `factor_qc_loadings.tsv` remains a filtered subset and was not used as a substitute."
  else
    paste0("The RDS contains exact zeros (n=", z_fh$n_exact_zero, "). These are included in the complete matrix."),
  "",
  paste0("Dimensions: ", nrow(Fh), " traits x ", ncol(Fh), " factors. Row order is the RDS `F_hat` row order."),
  paste0("Method in fit object: ", if (!is.null(fit101$method)) fit101$method else "NA"),
  ""
)
writeLines(readme, file.path(OUT, "FULL_GFA_LOADINGS_README.md"))

# ------------------- seed 101 vs 202 identity -------------------
logmsg("Reading seed 202 RDS for slot-wise comparison")
fit202 <- readRDS(RDS_202)

cmp_rows <- list()
slot_cmp <- function(name, a, b) {
  row <- data.table(
    slot = name,
    class_101 = paste(class(a), collapse = ","),
    class_202 = paste(class(b), collapse = ","),
    identical = NA,
    max_abs_diff = NA_real_,
    note = ""
  )
  if (is.null(a) && is.null(b)) {
    row$identical <- TRUE
    row$note <- "both NULL"
  } else if (is.numeric(a) || is.matrix(a) || is.logical(a)) {
    if (!is.null(dim(a)) || !is.null(dim(b))) {
      ma <- as.matrix(a); mb <- as.matrix(b)
      if (!identical(dim(ma), dim(mb))) {
        row$identical <- FALSE
        row$note <- paste("dim", paste(dim(ma), collapse = "x"), "vs", paste(dim(mb), collapse = "x"))
      } else {
        d <- max(abs(as.numeric(ma) - as.numeric(mb)), na.rm = TRUE)
        row$max_abs_diff <- d
        row$identical <- isTRUE(all.equal(ma, mb, tolerance = 0, check.attributes = FALSE))
        row$note <- paste0("attr_identical=", identical(dimnames(ma), dimnames(mb)))
      }
    } else if (length(a) == length(b)) {
      d <- max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE)
      row$max_abs_diff <- d
      row$identical <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0, check.attributes = FALSE))
    } else {
      row$identical <- FALSE
      row$note <- paste("length", length(a), "vs", length(b))
    }
  } else if (is.character(a) || is.factor(a) || is.list(a) || is.numeric(a) == FALSE) {
    row$identical <- identical(a, b)
    if (!isTRUE(row$identical)) {
      row$note <- paste0("typeof=", typeof(a), "/", typeof(b), " len=", length(a), "/", length(b))
    }
  }
  row
}

nms <- union(names(fit101), names(fit202))
for (nm in nms) {
  cmp_rows[[nm]] <- slot_cmp(nm, fit101[[nm]], fit202[[nm]])
}
cmp_dt <- rbindlist(cmp_rows, fill = TRUE)
fwrite(cmp_dt, file.path(OUT, "GFA_SEED101_VS_202_SLOT_COMPARE.tsv"), sep = "\t")
logmsg("slot compare written; n_slots=", nrow(cmp_dt),
       " n_identical=", sum(cmp_dt$identical, na.rm = TRUE))

# GFA package source: does gfa_fit / fixed_factors use RNG?
gfa_pkg <- tryCatch(find.package("GFA", quiet = TRUE), error = function(e) NA_character_)
logmsg("GFA package path: ", paste(gfa_pkg, collapse = ","))
src_out <- file.path(OUT, "GFA_PACKAGE_RNG_INSPECT.txt")
sink(src_out)
cat("GFA package:", paste(gfa_pkg, collapse = ","), "\n")
print(sessionInfo())
cat("\n===== names(fit101) =====\n")
print(names(fit101))
if (!is.null(fit101$method)) cat("method=", fit101$method, "\n")
if (!is.null(fit101$params)) {
  cat("\n===== fit101$params =====\n")
  str(fit101$params)
}
cat("\n===== getAnywhere / function dumps =====\n")
fns <- c("gfa_fit", "gfa_fixed_factors", "fixed_factors", "gfa_fit_fixed",
         "choose_k", "gfa_select", "initialize_factors")
ns <- tryCatch(asNamespace("GFA"), error = function(e) NULL)
if (!is.null(ns)) {
  cat("GFA exported:", paste(getNamespaceExports("GFA"), collapse = ", "), "\n")
  cat("GFA internal ls:\n")
  print(ls(ns, all.names = TRUE))
  for (fn in ls(ns, all.names = TRUE)) {
    obj <- ns[[fn]]
    if (!is.function(obj)) next
    txt <- paste(deparse(obj), collapse = "\n")
    hits <- grepl("set\\.seed|runif|rnorm|rgamma|sample\\(|\\.Random\\.seed|rbeta|rpois", txt)
    # always dump likely fit functions; dump others only if RNG mentioned
    keep <- grepl("fit|factor|init|fixed|gfa", fn, ignore.case = TRUE) || hits
    if (!keep) next
    cat("\n----- function:", fn, " RNG_like=", hits, " -----\n")
    cat(txt, "\n")
  }
} else {
  cat("Could not open GFA namespace\n")
}
sink()
logmsg("Wrote GFA inspect ", src_out)

# file times / sizes
fi <- rbind(
  data.table(file = "ukbb_ldclump_seed101.rds", as.list(file.info(RDS_101))),
  data.table(file = "ukbb_ldclump_seed202.rds", as.list(file.info(RDS_202))),
  data.table(file = "ukbb_ldclump_seed101_F_hat.tsv", as.list(file.info(TSV_101))),
  data.table(file = "ukbb_ldclump_seed202_F_hat.tsv", as.list(file.info(file.path(GFA, "ukbb_ldclump_seed202_F_hat.tsv"))))
)
fwrite(fi, file.path(OUT, "GFA_SEED_FILEINFO.tsv"), sep = "\t")

# stage copies of small outputs
for (fn in c(
  "FULL_GFA_LOADINGS_61x30.tsv",
  "FULL_GFA_LOADINGS_61x30_Fhat.tsv",
  "FULL_GFA_SCALE_VECTOR.tsv",
  "FULL_GFA_LOADINGS_README.md",
  "GFA_SEED101_VS_202_SLOT_COMPARE.tsv",
  "GFA_PACKAGE_RNG_INSPECT.txt",
  "GFA_SEED_FILEINFO.tsv"
)) {
  copy_also_stage(file.path(OUT, fn))
}

logmsg("DONE export + seed inspect")
