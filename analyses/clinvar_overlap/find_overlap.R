#!/usr/bin/env Rscript
options(echo = TRUE, error = function() { traceback(2); quit(status = 1) })

suppressPackageStartupMessages({
  library(config)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(readr)
  library(glue)
  library(fs)
  library(future.apply)
})

# -------------------------------
# 1. Load configuration
# -------------------------------
cfg <- config::get()
required_keys <- c("vcf_file", "bed_files", "output_folder", "log_file", "ignore_strand")
missing <- setdiff(required_keys, names(cfg))
if (length(missing) > 0) stop(glue("Missing config entries: {paste(missing, collapse = ', ')}"))

dir_create(cfg$output_folder)
log_file <- file.path(cfg$output_folder, cfg$log_file)
log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- glue("{timestamp} [{level}] {msg}\n")
  cat(full_msg, file = log_con, append = TRUE)
  message(glue("[{level}] {msg}"))
}

log_msg("Starting ClinVar-BED overlap analysis with variant-type stratification")

# -------------------------------
# 2. Load ClinVar VCF
# -------------------------------
log_msg(glue("Loading ClinVar VCF: {cfg$vcf_file}"))
vcf <- readVcf(cfg$vcf_file, param = ScanVcfParam(info = "CLNSIG"))
gr <- rowRanges(vcf)
info <- info(vcf)

if (!"CLNSIG" %in% colnames(info)) stop("VCF missing INFO/CLNSIG field")

# Flatten CLNSIG field
clnsig <- vapply(info(vcf)$CLNSIG, function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  paste0(unique(na.omit(as.character(x))), collapse = ",")
}, character(1L))

# Unique variant QC
variant_id <- sprintf("%s:%d_%s>%s", as.character(seqnames(gr)), start(gr),
                      as.character(ref(vcf)), as.character(unlist(alt(vcf))))
if (length(variant_id) != data.table::uniqueN(variant_id)) {
  log_msg("Duplicate variants detected in ClinVar VCF", level = "ERROR")
  stop("Please deduplicate the VCF before proceeding.")
}
log_msg(glue("QC PASS: {length(variant_id)} unique variants confirmed"))

# -------------------------------
# 3. Classify variant types (SNP vs INDEL ≤50bp vs LARGE_VAR)
# -------------------------------
ref_alleles <- as.character(ref(vcf))
alt_alleles <- as.character(unlist(alt(vcf)))
allele_diffs <- abs(nchar(ref_alleles) - nchar(alt_alleles))

variant_type <- ifelse(
  nchar(ref_alleles) == 1 & nchar(alt_alleles) == 1,
  "SNP",
  ifelse(allele_diffs < 50, "INDEL", "LARGE_VAR")
)

log_msg(glue("Variant type breakdown: {paste(names(table(variant_type)), table(variant_type), collapse = '; ')}"))

# -------------------------------
# 4. Convert GRanges → data.table
# -------------------------------
gr_dt <- data.table(
  VAR_ID   = sprintf("%s:%d_%s>%s",
                     as.character(seqnames(gr)), start(gr),
                     as.character(ref(vcf)), as.character(unlist(alt(vcf)))),
  chr       = as.character(seqnames(gr)),
  start     = start(gr),
  end       = end(gr),
  CLNSIG    = clnsig,
  VAR_TYPE  = variant_type
)

# Normalize chromosome names
gr_dt[, chr := sub("^chr", "", chr)]
setkey(gr_dt, chr, start, end)

# -------------------------------
# 5. CLNSIG and VAR_TYPE levels
# -------------------------------
clnsig_levels <- sort(unique(gr_dt$CLNSIG))
vartype_levels <- c("SNP", "INDEL", "LARGE_VAR")

log_msg(glue("CLNSIG values considered: {paste(clnsig_levels, collapse = ', ')}"))
log_msg(glue("Variant types considered: {paste(vartype_levels, collapse = ', ')}"))

# Total variants per CLNSIG × VAR_TYPE
total_per_group <- gr_dt[, .N, by = .(CLNSIG, VAR_TYPE)]
setnames(total_per_group, "N", "Total_variants")

log_msg(glue("Converted GRanges to data.table ({nrow(gr_dt)} variants)"))

# -------------------------------
# 6. Compute summary per BED file
# -------------------------------
compute_summary_fast <- function(bed_path, gr_dt, total_per_group) {
  bed_name <- gsub("\\.bed(\\.gz)?$", "", basename(bed_path))
  log_msg(glue("Processing BED: {bed_name}"))
  start_time <- Sys.time()

  # Read BED or BED.gz
  cmd_line <- if (grepl("\\.gz$", bed_path)) glue("zcat {bed_path}") else glue("cat {bed_path}")
  bed <- tryCatch({
    tmp <- fread(cmd = cmd_line, data.table = TRUE)
    tmp <- tmp[!grepl("^#", tmp[[1]])]
    bed <- tmp[, .(chr = tmp[[1]], start = tmp[[2]], end = tmp[[3]])]
    bed[, chr := sub("^chr", "", chr)]
    
    # Filter invalid regions
    bed <- bed[end > start]

  }, error = function(e) {
    log_msg(glue("  ERROR reading {bed_name}: {e$message}"), level = "ERROR")
    return(NULL)
  })

  if (is.null(bed) || nrow(bed) == 0L) {
    log_msg(glue("  Warning: {bed_name} has 0 regions, skipping."), level = "WARN")
    return(data.frame(
      BED_file = bed_name,
      CLNSIG = total_per_group$CLNSIG,
      VAR_TYPE = total_per_group$VAR_TYPE,
      Overlap_variants = 0,
      Non_overlap_variants = total_per_group$Total_variants,
      Total_variants = total_per_group$Total_variants,
      Overlap_proportion = 0
    ))
  }

  setkey(bed, chr, start, end)
  ov <- foverlaps(gr_dt, bed, nomatch = 0L)

  # Count unique variants per CLNSIG × VAR_TYPE
  overlap_counts <- ov[, .(Overlap_variants = uniqueN(VAR_ID)), by = .(CLNSIG, VAR_TYPE)]

  summary_dt <- merge(total_per_group, overlap_counts,
                      by = c("CLNSIG", "VAR_TYPE"), all.x = TRUE)
  summary_dt[is.na(Overlap_variants), Overlap_variants := 0L]
  summary_dt[, Non_overlap_variants := Total_variants - Overlap_variants]
  summary_dt[, Overlap_proportion := round(Overlap_variants / Total_variants, 4)]
  summary_dt[, BED_file := bed_name]

  total_hits <- sum(summary_dt$Overlap_variants)
  total_prop <- round(total_hits / sum(summary_dt$Total_variants) * 100, 2)
  duration <- round(difftime(Sys.time(), start_time, units = "secs"), 1)
  log_msg(glue("  Overlapping variants: {total_hits}/{sum(summary_dt$Total_variants)} ({total_prop}%) in {duration}s"))

  summary_dt[, .(BED_file, VAR_TYPE, CLNSIG, Overlap_variants,
                 Non_overlap_variants, Total_variants, Overlap_proportion)]
}

# -------------------------------
# 7. Run in parallel or sequential
# -------------------------------
parallelize <- if (!is.null(cfg$parallel)) cfg$parallel else FALSE
if (parallelize) {
  plan(multisession)
  log_msg("Running overlap analysis in parallel mode")
  summary_list <- future_lapply(cfg$bed_files, compute_summary_fast,
                                gr_dt = gr_dt, total_per_group = total_per_group)
} else {
  log_msg("Running overlap analysis in sequential mode")
  summary_list <- lapply(cfg$bed_files, compute_summary_fast,
                         gr_dt = gr_dt, total_per_group = total_per_group)
}

summary_df <- rbindlist(summary_list, fill = TRUE)

# -------------------------------
# 8. Order and output
# -------------------------------
summary_df[, CLNSIG := factor(CLNSIG, levels = clnsig_levels)]
summary_df[, VAR_TYPE := factor(VAR_TYPE, levels = vartype_levels)]
summary_df <- summary_df[order(BED_file, VAR_TYPE, CLNSIG)]

output_file <- file.path(cfg$output_folder, "clinvar_overlap_by_CLNSIG_and_VAR_TYPE.csv")
write_csv(summary_df, output_file)
log_msg(glue("Saved ClinVar overlap summary: {output_file}"))

# -------------------------------
# 9. Optional per-variant-type summary
# -------------------------------
type_summary <- summary_df[, .(
  Total_variants = sum(Total_variants),
  Overlap_variants = sum(Overlap_variants),
  Overlap_proportion = round(sum(Overlap_variants) / sum(Total_variants), 4)
), by = .(BED_file, VAR_TYPE)]

write_csv(type_summary, file.path(cfg$output_folder, "clinvar_overlap_by_VAR_TYPE.csv"))
log_msg("Saved variant-type summary table")

log_msg(glue("Output folder: {cfg$output_folder}"))
log_msg("Done.")
