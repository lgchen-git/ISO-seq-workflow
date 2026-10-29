#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script <- if (length(args)) args[[1L]] else "pipeline/scripts/differential.R"
script <- normalizePath(script, mustWork = TRUE)
required <- c("edgeR", "DRIMSeq", "stageR", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing statistical packages: ", paste(missing, collapse = ", "))

root <- tempfile("isoseq-differential-test-")
main <- function() {
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
samples <- expand.grid(replicate = seq_len(4L), group = c("Control", "DrugA", "DrugB"),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
samples$sample_id <- sprintf("S%02d", seq_len(nrow(samples)))
samples$batch <- rep(c("B1", "B2", "B1", "B2"), times = 3L)
samples <- samples[, c("sample_id", "group", "batch")]
samples$fastq <- paste0("inputs/", samples$sample_id, ".fastq.gz")
samples <- samples[, c("sample_id", "fastq", "group", "batch")]
contrasts <- data.frame(contrast_id = c("DrugA_vs_Control", "DrugB_vs_Control"),
                        test_group = c("DrugA", "DrugB"),
                        reference_group = "Control")
write.table(samples, file.path(root, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contrasts, file.path(root, "contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gene_ids <- sprintf("gene%02d", seq_len(60L))
tx_ids <- unlist(lapply(gene_ids, function(g) paste0(g, "_tx", 1:2)), use.names = FALSE)
gene_counts <- matrix(0L, nrow = length(gene_ids), ncol = nrow(samples),
                      dimnames = list(gene_ids, samples$sample_id))
tx_counts <- matrix(0L, nrow = length(tx_ids), ncol = nrow(samples),
                    dimnames = list(tx_ids, samples$sample_id))
for (g in seq_along(gene_ids)) {
  for (i in seq_len(nrow(samples))) {
    group <- samples$group[[i]]
    total <- 200L + 5L * g + if (group == "DrugA" && g <= 12L) 200L else 0L
    baseline <- if (g <= 24L && group == "DrugA") 0.25 else if (g <= 24L && group == "DrugB") 0.5 else 0.8
    proportion <- baseline + c(-0.06, -0.02, 0.02, 0.06)[(i - 1L) %% 4L + 1L]
    first <- as.integer(round(total * proportion))
    gene_counts[g, i] <- total
    tx_counts[2L * g - 1L, i] <- first
    tx_counts[2L * g, i] <- total - first
  }
}
gene_table <- data.frame(gene_id = gene_ids, gene_counts, check.names = FALSE)
tx_table <- data.frame(transcript_id = tx_ids, tx_counts, check.names = FALSE)
catalog <- data.frame(transcript_id = tx_ids,
                      gene_id = rep(gene_ids, each = 2L),
                      filter_status = "PASS")
write.table(gene_table, file.path(root, "gene_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tx_table, file.path(root, "transcript_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(catalog, file.path(root, "catalog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

out <- file.path(root, "out")
command_args <- c(script,
                  "--gene-counts", file.path(root, "gene_counts.tsv"),
                  "--transcript-counts", file.path(root, "transcript_counts.tsv"),
                  "--catalog", file.path(root, "catalog.tsv"),
                  "--samples", file.path(root, "samples.tsv"),
                  "--contrasts", file.path(root, "contrasts.tsv"),
                  "--output-dir", out)
status <- system2(file.path(R.home("bin"), "Rscript"),
                  args = vapply(command_args, shQuote, character(1)),
                  stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status")
if (!is.null(exit_status) && exit_status != 0L) stop(paste(status, collapse = "\n"))
summary <- jsonlite::read_json(file.path(out, "summary.json"), simplifyVector = TRUE)
n_contrasts <- if (is.data.frame(summary$contrasts)) nrow(summary$contrasts) else length(summary$contrasts)
stopifnot(isTRUE(summary$stats_enabled), n_contrasts == 2L)
gene_result <- read.delim(file.path(out, "gene_de.tsv"), check.names = FALSE)
tx_result <- read.delim(file.path(out, "transcript_de.tsv"), check.names = FALSE)
dtu_result <- read.delim(file.path(out, "dtu_transcripts.tsv"), check.names = FALSE)
stopifnot(all(c("Control", "DrugA", "DrugB") %in% samples$group))
stopifnot(setequal(unique(gene_result$contrast_id), contrasts$contrast_id))
stopifnot(setequal(unique(tx_result$contrast_id), contrasts$contrast_id))
stopifnot(all(c("logit_change", "stageR_gene", "stageR_transcript", "stageR_significant") %in% names(dtu_result)))
cat("Differential-expression and DTU multi-group/batch fixture passed\n")
}
main()
