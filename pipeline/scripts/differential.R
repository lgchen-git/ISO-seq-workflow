#!/usr/bin/env Rscript

stopf <- function(...) stop(sprintf(...), call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) stopf("Arguments must be supplied as --key value pairs")
opts <- as.list(args[seq(2L, length(args), by = 2L)])
names(opts) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])
required <- c("gene-counts", "transcript-counts", "catalog", "samples", "contrasts", "output-dir")
missing <- setdiff(required, names(opts))
if (length(missing)) stopf("Missing arguments: %s", paste(missing, collapse = ", "))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                      na.strings = character(), quote = "", comment.char = "")
write_tsv <- function(x, path, columns = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(columns)) columns <- names(x)
  if (is.null(x) || nrow(x) == 0L) x <- as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns))
  write.table(x[, columns, drop = FALSE], path, sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA")
}
extract_ids <- function(x, preferred = NULL) {
  if (is.null(x) || !length(x)) return(character())
  if (is.data.frame(x) || is.matrix(x)) {
    nms <- colnames(x)
    if (!is.null(preferred) && !is.null(nms) && preferred %in% nms) {
      x <- x[, preferred, drop = TRUE]
    } else if (ncol(x) > 0L) {
      x <- x[, 1L, drop = TRUE]
    } else {
      return(character())
    }
  }
  if (is.null(x) || !length(x)) character() else as.character(x)
}
as_count_matrix <- function(path, feature_col, sample_ids) {
  data <- read_tsv(path)
  if (!feature_col %in% names(data) || !all(sample_ids %in% names(data))) {
    stopf("Count matrix %s is missing required feature/sample columns", path)
  }
  ids <- data[[feature_col]]
  values <- as.matrix(data[, sample_ids, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) stopf("Invalid counts in %s", path)
  if (any(abs(values - round(values)) > 1e-8)) stopf("Non-integer counts in %s; edgeR/DRIMSeq require read counts", path)
  storage.mode(values) <- "integer"
  rownames(values) <- ids
  values
}

dir.create(opts[["output-dir"]], recursive = TRUE, showWarnings = FALSE)
samples <- read_tsv(opts$samples)
contrasts <- read_tsv(opts$contrasts)
catalog <- read_tsv(opts$catalog)
sample_ids <- samples$sample_id
gene_counts <- as_count_matrix(opts[["gene-counts"]], "gene_id", sample_ids)
tx_counts <- as_count_matrix(opts[["transcript-counts"]], "transcript_id", sample_ids)

gene_output <- list()
tx_output <- list()
dtu_gene_output <- list()
dtu_tx_output <- list()
contrast_summary <- list()

empty_gene <- c("contrast_id", "feature_id", "logFC", "logCPM", "F", "PValue", "FDR")
empty_dtu_gene <- c("contrast_id", "gene_id", "lr", "df", "pvalue", "adj_pvalue", "stageR_gene", "stageR_significant")
empty_dtu_tx <- c("contrast_id", "gene_id", "feature_id", "logit_change", "lr", "df", "pvalue", "adj_pvalue", "stageR_gene", "stageR_transcript", "stageR_significant")

run_edger <- function(counts, samples_sub, contrast_id, feature_type) {
  if (nrow(counts) == 0L) return(data.frame())
  y <- edgeR::DGEList(counts = counts)
  design <- stats::model.matrix(if ("batch" %in% names(samples_sub) && length(unique(samples_sub$batch)) > 1L) {
    ~ batch + group
  } else {
    ~ group
  }, data = samples_sub)
  keep <- edgeR::filterByExpr(y, design = design)
  if (!any(keep)) return(data.frame())
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coefficient <- grep("^group", colnames(design))
  if (length(coefficient) != 1L) stopf("Contrast %s: expected one group coefficient in edgeR design", contrast_id)
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  out <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  out$feature_id <- rownames(out)
  out$contrast_id <- contrast_id
  out[, c("contrast_id", "feature_id", "logFC", "logCPM", "F", "PValue", "FDR"), drop = FALSE]
}

if (nrow(contrasts)) {
  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("DRIMSeq", quietly = TRUE) ||
      !requireNamespace("stageR", quietly = TRUE)) {
    stopf("Statistical dependencies edgeR, DRIMSeq and stageR must be installed in the workflow image")
  }
  for (i in seq_len(nrow(contrasts))) {
    contrast <- contrasts[i, , drop = FALSE]
    cid <- contrast$contrast_id[[1L]]
    ref <- contrast$reference_group[[1L]]
    test_group <- contrast$test_group[[1L]]
    samples_sub <- samples[samples$group %in% c(ref, test_group), , drop = FALSE]
    samples_sub$group <- factor(samples_sub$group, levels = c(ref, test_group))
    if ("batch" %in% names(samples_sub) && length(unique(samples_sub$batch)) > 1L) {
      samples_sub$batch <- factor(samples_sub$batch)
    }
    rownames(samples_sub) <- samples_sub$sample_id
    columns <- samples_sub$sample_id
    contrast_gene_counts <- gene_counts[, columns, drop = FALSE]
    contrast_tx_counts <- tx_counts[, columns, drop = FALSE]

    gene_result <- run_edger(contrast_gene_counts, samples_sub, cid, "gene")
    tx_result <- run_edger(contrast_tx_counts, samples_sub, cid, "transcript")
    gene_output[[length(gene_output) + 1L]] <- gene_result
    tx_output[[length(tx_output) + 1L]] <- tx_result

    catalog_sub <- catalog[match(rownames(contrast_tx_counts), catalog$transcript_id), , drop = FALSE]
    keep <- !is.na(catalog_sub$gene_id) & catalog_sub$gene_id != "" & catalog_sub$gene_id != "."
    long_counts <- data.frame(
      feature_id = rownames(contrast_tx_counts)[keep],
      gene_id = catalog_sub$gene_id[keep],
      contrast_tx_counts[keep, , drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    genes_before <- length(unique(long_counts$gene_id))
    dtu_gene <- data.frame()
    dtu_tx <- data.frame()
    if (nrow(long_counts)) {
      d <- DRIMSeq::dmDSdata(counts = long_counts, samples = samples_sub[, c("sample_id", "group", intersect("batch", names(samples_sub))), drop = FALSE])
      d <- DRIMSeq::dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_expr = 2,
                            min_gene_expr = 10, min_feature_expr = 3)
      if (nrow(DRIMSeq::counts(d)) > 0L) {
        design <- stats::model.matrix(if ("batch" %in% names(samples_sub) && length(unique(samples_sub$batch)) > 1L) {
          ~ batch + group
        } else {
          ~ group
        }, data = DRIMSeq::samples(d))
        set.seed(1108)
        d <- DRIMSeq::dmPrecision(d, design = design)
        d <- DRIMSeq::dmFit(d, design = design, verbose = 0)
        coefficient <- grep("^group", colnames(design))
        if (length(coefficient) != 1L) stopf("Contrast %s: expected one group coefficient in DRIMSeq design", cid)
        d <- DRIMSeq::dmTest(d, coef = coefficient, verbose = 0)
        screen <- DRIMSeq::results(d, level = "gene")
        confirm <- DRIMSeq::results(d, level = "feature")
        valid_confirm <- is.finite(confirm$pvalue)
        valid_screen <- is.finite(screen$pvalue) & screen$gene_id %in% confirm$gene_id[valid_confirm]
        screen <- screen[valid_screen, , drop = FALSE]
        valid_confirm <- valid_confirm & confirm$gene_id %in% screen$gene_id
        confirm <- confirm[valid_confirm, , drop = FALSE]
        if (nrow(screen) && nrow(confirm)) {
          p_screen <- screen$pvalue
          names(p_screen) <- screen$gene_id
          p_confirm <- matrix(confirm$pvalue, ncol = 1L,
                             dimnames = list(confirm$feature_id, "condition"))
          stage <- stageR::stageRTx(
            pScreen = p_screen,
            pConfirmation = p_confirm,
            pScreenAdjusted = FALSE,
            tx2gene = unique(confirm[, c("feature_id", "gene_id")])
          )
          stage <- stageR::stageWiseAdjustment(stage, method = "dtu", alpha = 0.05)
          stage_adjusted <- stageR::getAdjustedPValues(stage, order = FALSE, onlySignificantGenes = FALSE)
          stage_frame <- as.data.frame(stage_adjusted, check.names = FALSE)
          if ("txID" %in% names(stage_frame)) {
            stage_frame$feature_id <- stage_frame$txID
          } else if ("transcript" %in% names(stage_frame)) {
            stage_frame$feature_id <- stage_frame$transcript
          } else if ("feature_id" %in% names(stage_frame)) {
            stage_frame$feature_id <- stage_frame$feature_id
          } else {
            stage_frame$feature_id <- rownames(stage_frame)
          }
          stage_columns <- setdiff(names(stage_frame), c("txID", "geneID", "feature_id", "gene_id"))
          names(stage_frame)[match(stage_columns, names(stage_frame))] <- paste0("stageR_", stage_columns)
          significant_genes <- extract_ids(stageR::getSignificantGenes(stage))
          dtu_gene <- screen
          if ("geneID" %in% names(stage_frame) && "stageR_gene" %in% names(stage_frame)) {
            dtu_gene$stageR_gene <- stage_frame$stageR_gene[match(dtu_gene$gene_id, stage_frame$geneID)]
          } else {
            dtu_gene$stageR_gene <- NA_real_
          }
          dtu_gene$stageR_significant <- dtu_gene$gene_id %in% significant_genes
          dtu_gene$contrast_id <- cid
          dtu_gene <- dtu_gene[, c("contrast_id", "gene_id", "lr", "df", "pvalue", "adj_pvalue", "stageR_gene", "stageR_significant"), drop = FALSE]
          significant_tx <- stageR::getSignificantTx(stage)
          significant_tx_ids <- if (is.data.frame(significant_tx) || is.matrix(significant_tx)) {
            if ("txID" %in% colnames(significant_tx)) {
              extract_ids(significant_tx, "txID")
            } else if (nrow(significant_tx) > 0L) {
              rownames(significant_tx)
            } else {
              character()
            }
          } else {
            extract_ids(significant_tx)
          }
          stage_columns <- setdiff(names(stage_frame), "feature_id")
          dtu_tx <- merge(confirm, stage_frame[, c("feature_id", stage_columns), drop = FALSE], by = "feature_id", all.x = TRUE, sort = FALSE)
          feature_coefficients <- as.data.frame(DRIMSeq::coefficients(d, level = "feature"), check.names = FALSE)
          coefficient_name <- colnames(design)[coefficient]
          if (!"feature_id" %in% names(feature_coefficients) || !coefficient_name %in% names(feature_coefficients)) {
            stopf("Contrast %s: DRIMSeq did not return the expected transcript group coefficient", cid)
          }
          dtu_tx$logit_change <- feature_coefficients[[coefficient_name]][match(dtu_tx$feature_id, feature_coefficients$feature_id)]
          dtu_tx$contrast_id <- cid
          dtu_tx$stageR_significant <- dtu_tx$feature_id %in% significant_tx_ids
          dtu_tx <- dtu_tx[, c("contrast_id", "gene_id", "feature_id", "logit_change", "lr", "df", "pvalue", "adj_pvalue", stage_columns, "stageR_significant"), drop = FALSE]
        }
      }
    }
    dtu_gene_output[[length(dtu_gene_output) + 1L]] <- dtu_gene
    dtu_tx_output[[length(dtu_tx_output) + 1L]] <- dtu_tx
    contrast_summary[[length(contrast_summary) + 1L]] <- list(
      contrast_id = cid,
      test_group = test_group,
      reference_group = ref,
      sample_count = nrow(samples_sub),
      samples_by_group = list(reference = sum(samples_sub$group == ref), test = sum(samples_sub$group == test_group)),
      batch_levels = if ("batch" %in% names(samples_sub) && any(nzchar(as.character(samples_sub$batch)))) length(unique(samples_sub$batch)) else 0L,
      batch_adjusted = "batch" %in% names(samples_sub) && length(unique(samples_sub$batch)) > 1L,
      gene_de_features = nrow(gene_result),
      transcript_de_features = nrow(tx_result),
      dtu_genes_tested = nrow(dtu_gene),
      dtu_transcripts_tested = nrow(dtu_tx),
      dtu_genes_with_multiple_transcripts_before_filter = genes_before,
      low_power_warning = min(table(samples_sub$group)) < 3L
    )
  }
}

bind_or_empty <- function(items, columns) {
  nonempty <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, items)
  if (length(nonempty)) return(do.call(rbind, nonempty))
  as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns))
}
gene_table <- bind_or_empty(gene_output, empty_gene)
tx_table <- bind_or_empty(tx_output, empty_gene)
dtu_gene_table <- bind_or_empty(dtu_gene_output, empty_dtu_gene)
dtu_tx_columns <- unique(c(empty_dtu_tx, unlist(lapply(dtu_tx_output, function(x) setdiff(names(x), empty_dtu_tx)), use.names = FALSE)))
dtu_tx_table <- bind_or_empty(dtu_tx_output, dtu_tx_columns)
write_tsv(gene_table, file.path(opts[["output-dir"]], "gene_de.tsv"))
write_tsv(tx_table, file.path(opts[["output-dir"]], "transcript_de.tsv"))
write_tsv(dtu_gene_table, file.path(opts[["output-dir"]], "dtu_genes.tsv"))
write_tsv(dtu_tx_table, file.path(opts[["output-dir"]], "dtu_transcripts.tsv"))

summary <- list(
  stats_enabled = nrow(contrasts) > 0L,
  status = if (nrow(contrasts)) "complete" else "descriptive_only",
  method = list(gene_and_transcript_de = "edgeR robust quasi-likelihood negative-binomial GLM",
                dtu = "DRIMSeq Dirichlet-multinomial with stageR stage-wise adjustment"),
  alpha = 0.05,
  counts = list(gene_de_rows = nrow(gene_table), transcript_de_rows = nrow(tx_table),
                dtu_gene_rows = nrow(dtu_gene_table), dtu_transcript_rows = nrow(dtu_tx_table)),
  contrasts = contrast_summary,
  package_versions = setNames(lapply(c("edgeR", "DRIMSeq", "stageR"), function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "not available"
  }), c("edgeR", "DRIMSeq", "stageR"))
)
jsonlite::write_json(summary, file.path(opts[["output-dir"]], "summary.json"), auto_unbox = TRUE,
                     pretty = TRUE, null = "null", na = "null")
